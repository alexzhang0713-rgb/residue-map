"""Reproducible CLI comparison. Optional development dependency: pdb-tools==2.7.0.

python tools/compare_pdbtools.py --package-dir work/pdbtools --output-dir work/comparison
Only ATOM lines enter the comparison pipelines, matching Residue Map's export scope.
No application runtime dependency is added.
"""
import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from core import parse_pdb, export_crop


def atom(serial, number, insertion='', chain='A', alt='', name='ALA', x=1.0):
    return f'ATOM  {serial:5d}  CA {alt:1}{name:3} {chain:1}{number:4d}{insertion:1}   {x:8.3f}{2.:8.3f}{3.:8.3f}{1.:6.2f}{20.:6.2f}           C  '


def atoms(text):
    return [l.ljust(80)[:80] for l in text.splitlines() if l.startswith('ATOM  ')]


def atom_scope(text):
    return '\n'.join(atoms(text)) + '\n'


def labels(text):
    result = []
    for l in atoms(text):
        item = (l[21], int(l[22:26]), l[26].strip(), l[17:20].strip())
        if not result or result[-1] != item:
            result.append(item)
    return result


def identity_free(text):
    # All atom fields except atom serial and residue number/insertion code.
    return [(l[:6], l[11:22], l[27:]) for l in atoms(text)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package-dir', type=Path)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'work/comparison')
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    if args.package_dir:
        env['PYTHONPATH'] = str(args.package_dir.resolve())
    version = subprocess.run([sys.executable, '-c',
        "from importlib.metadata import version; print(version('pdb-tools'))"],
        env=env, text=True, capture_output=True, check=True).stdout.strip()
    commands = []

    def cli(module, options, text):
        command = [sys.executable, '-m', 'pdbtools.'+module, *options]
        proc = subprocess.run(command, input=text, text=True, capture_output=True, env=env)
        commands.append({'module': module, 'options': options, 'returncode': proc.returncode,
                         'stderr': proc.stderr, 'input_sha256': hashlib.sha256(text.encode()).hexdigest()})
        return proc

    def pipeline(text, steps):
        for module, options in steps:
            p = cli(module, options, text)
            if p.returncode:
                raise RuntimeError(f'{module}: {p.stderr}')
            text = p.stdout
        return text

    crn = (ROOT/'tests/fixtures/1CRN.pdb').read_text()
    igt = (ROOT/'tests/fixtures/1IGT.pdb').read_text()
    demo = (ROOT/'examples/numbering_demo.pdb').read_text()
    cases = [
        ('1crn_full', crn, 'A', '1', '46'),
        ('1crn_10_20', crn, 'A', '10', '20'),
        ('1igt_B_full', igt, 'B', '1', '474'),
        ('1igt_B_82_83', igt, 'B', '82', '83'),
        ('1igt_D_82_83', igt, 'D', '82', '83'),
        ('1igt_B_82A_83', igt, 'B', '82A', '83'),
        ('1igt_B_100H_100K', igt, 'B', '100H', '100K'),
        ('demo_125_150', demo, 'A', '125', '150'),
        ('demo_125A_150', demo, 'A', '125A', '150'),
        ('negative_chain', demo, 'B', '-1', '0'),
        ('alternate_atoms', atom(1,5,alt='A')+'\n'+atom(2,5,alt='B')+'\n'+atom(3,9), 'A', '5', '9'),
        ('nonmonotonic', '\n'.join([atom(1,10),atom(2,30),atom(3,20)]), 'A', '10', '20'),
    ]
    results = []
    for name, text, chain, first, last in cases:
        rs, _ = parse_pdb(text)
        selected = [r for r in rs if r.chain == chain]
        start = next(i for i,r in enumerate(selected) if r.label == first)
        end = next(i for i,r in enumerate(selected) if r.label == last)
        ours, table, mapping = export_crop(rs, chain, start, end)
        base = [('pdb_selchain', ['-'+chain])]
        # Natural numeric pipeline: records sharing an integer all get selected;
        # reres preserves insertion codes, so comparison is deliberately fieldwise.
        low = selected[start].number
        high = selected[end].number
        numeric = pipeline(atom_scope(text), base+[
            ('pdb_selres', [f'-{low}:{high}']), ('pdb_reres', ['-1']), ('pdb_reatom',['-1'])])
        # Equivalent recipe: clear insertion codes, assign observed ordinal indices,
        # select inclusive ordinal bounds, then restart residue/atom numbering at 1.
        equivalent = pipeline(atom_scope(text), base+[
            ('pdb_fixinsert', []), ('pdb_reres',['-1']),
            ('pdb_selres',[f'-{start+1}:{end+1}']),
            ('pdb_reres',['-1']), ('pdb_reatom',['-1'])])
        (out/(name+'.input.pdb')).write_text(text)
        (out/(name+'.residuemap.pdb')).write_text(ours)
        (out/(name+'.pdbtools-numeric.pdb')).write_text(numeric)
        (out/(name+'.pdbtools-equivalent.pdb')).write_text(equivalent)
        (out/(name+'.mapping.csv')).write_text(table)
        expected_atoms = [l for r in selected[start:end+1] for l in r.lines]
        expected_text = '\n'.join(expected_atoms)
        results.append({'case':name,'source_chain':chain,'first':first,'last':last,
            'observed_bounds':[start+1,end+1], 'residues':len(mapping),'atoms':len(atoms(ours)),
            'ours_preserves_selected_atom_fields':identity_free(ours)==identity_free(expected_text),
            'mapping_csv_roundtrip':[(int(r['original_resseq']),r['original_icode']) for r in csv.DictReader(io.StringIO(table))]==[(r.number,r.insertion) for r in selected[start:end+1]],
            'numeric_same_atoms':identity_free(ours)==identity_free(numeric),
            'numeric_same_residue_ids':labels(ours)==labels(numeric),
            'numeric_labels':labels(numeric),
            'equivalent_same_all_atom_columns':atoms(ours)==atoms(equivalent),
            'equivalent_same_residue_ids':labels(ours)==labels(equivalent)})
    # Behavior probes are observations, not expectations that filters validate input.
    probes = []
    for name,text in [
        ('duplicate_id', '\n'.join([atom(1,1),atom(2,2),atom(3,1)])),
        ('mixed_residue_names', atom(1,1)+'\n'+atom(2,1,name='GLY')),
        ('nan_coordinate', atom(1,1,x=float('nan'))),
        ('multiple_models', 'MODEL        1\n'+atom(1,1)+'\nENDMDL\nMODEL        2\n'+atom(1,1)+'\nENDMDL'),
    ]:
        try:
            parse_pdb(text); status='accepted'
        except ValueError as exc:
            status=str(exc)
        # Full input (including MODEL records) intentionally retained in probes.
        res = cli('pdb_reres',['-1'],text+'\n')
        probes.append({'probe':name,'residuemap':status,'pdb_reres_returncode':res.returncode,
                       'pdb_reres_output_atom_count':len(atoms(res.stdout)),'stderr':res.stderr})
        (out/(name+'.input.pdb')).write_text(text+'\n')
        (out/(name+'.pdbtools.pdb')).write_text(res.stdout)
    blank = atom(1,-1,chain=' ')+'\n'+atom(2,0,chain=' ')+'\n'
    blank_rs,_=parse_pdb(blank)
    blank_out,_,_=export_crop(blank_rs,' ',0,1)
    blank_cli=cli('pdb_selchain',['- '],blank)
    # API demonstrates this is a CLI argument limitation, not inability to handle a blank chain.
    api_probe=subprocess.run([sys.executable,'-c',
        "from pdbtools.pdb_selchain import run; import sys; sys.stdout.writelines(run(sys.stdin, {' '}))"],
        input=blank,text=True,capture_output=True,env=env,check=True)
    probes.append({'probe':'blank_chain','residuemap_labels':labels(blank_out),
                   'pdb_selchain_cli_returncode':blank_cli.returncode,
                   'pdb_selchain_api_preserves_input':atoms(api_probe.stdout)==atoms(blank),
                   'stderr':blank_cli.stderr.splitlines()[0]})
    insert_cli=cli('pdb_selres',['-82A:83'],atom_scope(igt))
    probes.append({'probe':'insertion_label_cli','command':'pdb_selres -82A:83',
                   'returncode':insert_cli.returncode,'stderr':insert_cli.stderr.strip()})
    # Exercise the actual HTTP export, including its sidecar report and rejection.
    import threading
    import urllib.request
    import urllib.error
    from http.server import ThreadingHTTPServer
    from server import Handler
    http_server = ThreadingHTTPServer(('127.0.0.1',0),Handler)
    thread = threading.Thread(target=http_server.serve_forever,daemon=True)
    thread.start()
    b_chain=[r for r in parse_pdb(igt)[0] if r.chain=='B']
    payload={'pdb':igt,'chain':'B','start':next(i for i,r in enumerate(b_chain) if r.label=='82'),
             'end':next(i for i,r in enumerate(b_chain) if r.label=='83'),'hotspots':'52A,82A'}
    def post(path, body):
        req=urllib.request.Request(f'http://127.0.0.1:{http_server.server_port}/api/'+path,
            data=json.dumps(body).encode(), headers={'Content-Type':'application/json'})
        with urllib.request.urlopen(req) as response:
            return json.load(response)
    try:
        checked=post('check',{'pdb':igt,'chain':'B','hotspots':'52A,82A,82D'})
        exported=post('export',payload)
        assert [{k:h[k] for k in ('original','new','status')} for h in exported['report']['hotspots']]==[
            {'original':'52A','new':None,'status':'outside crop'},
            {'original':'82A','new':2,'status':'retained'}]
        assert exported['report']['input_sha256']==hashlib.sha256(igt.encode()).hexdigest()
        for filename,content in [('1igt-B-hotspots.report.json',json.dumps(exported['report'],indent=2)),
                                 ('1igt-B-hotspots.mapping.csv',exported['csv']),
                                 ('1igt-B-hotspots.cropped.pdb',exported['pdb'])]:
            (out/filename).write_text(content)
        payload['hotspots']='82D'
        try:
            post('export',payload)
            raise AssertionError('missing hotspot must reject export')
        except urllib.error.HTTPError as exc:
            rejected_status=exc.code
            assert rejected_status==400
        probes.append({'probe':'hotspots_http','checks':checked,'mapping':exported['report']['hotspots'],
            'absent_hotspot_export_status':rejected_status,'input_sha256_verified':True})
    finally:
        http_server.shutdown(); http_server.server_close(); thread.join()
    # Raw-input metadata probe: don't conflate default retention with inability to exclude it.
    raw_pipeline=pipeline(igt,[('pdb_reres',['-1'])])
    no_het=pipeline(raw_pipeline,[('pdb_delhetatm',[])])
    hetero=lambda text:sum(l.startswith('HETATM') for l in text.splitlines())
    probes.append({'probe':'raw_input_metadata','pdbtools_full_structure_hetatm':hetero(raw_pipeline),
                   'after_pdb_delhetatm':hetero(no_het),
                   'residuemap_export_hetatm':hetero(exported['pdb']),
                   'residuemap_exclusion_warning':any('HETATM' in w for w in exported['report']['warnings']),
                   'note':'Raw full-structure probe is separate from the crop comparison; not an assertion of identical selections.'})
    success=all(r['ours_preserves_selected_atom_fields'] and r['mapping_csv_roundtrip'] and r['equivalent_same_all_atom_columns'] for r in results)
    report={'pdb_tools_version':version,'python':sys.version.split()[0],
        'scope':'ATOM-only normalized input for crop pipelines; compare 80 PDB ATOM columns, not metadata/TER/END/trailing whitespace.',
        'fixtures':{name:hashlib.sha256((ROOT/'tests/fixtures'/name).read_bytes()).hexdigest() for name in ['1CRN.pdb','1IGT.pdb']},
        'cases':results,'behavior_probes':probes,'cli_invocations':commands,
        'equivalent_cases_passed':sum(r['equivalent_same_all_atom_columns'] for r in results),
        'success':success}
    (out/'results.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('cases','cli_invocations','behavior_probes')},indent=2))
    for r in results:
        print(r['case'],r['residues'],r['atoms'],'numeric atoms/ids',r['numeric_same_atoms'],r['numeric_same_residue_ids'],'equivalent',r['equivalent_same_all_atom_columns'])
    for probe in probes: print(json.dumps(probe,ensure_ascii=False))
    return 0 if success else 1

if __name__=='__main__':
    sys.exit(main())

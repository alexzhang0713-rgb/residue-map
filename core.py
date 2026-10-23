"""PDB ATOM selection and explicit residue mapping. Python 3.10+, stdlib only."""
from dataclasses import dataclass, field
import csv
import io
import math
import re

AA = dict(zip('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(), 'ARNDCQEGHILKMFPSTWYV'))

@dataclass
class Residue:
    chain: str
    number: int
    insertion: str
    name: str
    lines: list[str] = field(default_factory=list)

    @property
    def label(self):
        return f'{self.number}{self.insertion}'

    @property
    def key(self):
        return self.chain, self.number, self.insertion


def parse_pdb(text):
    residues, by_key, warnings = [], {}, []
    model_count = 0
    active_key = None
    atom_count = 0
    for lineno, raw in enumerate(text.splitlines(), 1):
        record = raw[:6].strip()
        if record == 'MODEL':
            model_count += 1
            if model_count > 1:
                raise ValueError('Multiple MODEL records: extract one model before importing.')
        if record in ('TER', 'ENDMDL'):
            active_key = None
        if record == 'HETATM':
            if 'HETATM records excluded (including waters, ligands and modified residues).' not in warnings:
                warnings.append('HETATM records excluded (including waters, ligands and modified residues).')
        if record != 'ATOM':
            continue
        if len(raw) < 54:
            raise ValueError(f'Line {lineno}: truncated ATOM record.')
        line = raw.ljust(80)
        try:
            number = int(line[22:26])
            coords = [float(line[a:b]) for a, b in ((30,38),(38,46),(46,54))]
            if not all(math.isfinite(v) for v in coords):
                raise ValueError()
        except ValueError:
            raise ValueError(f'Line {lineno}: invalid residue number or coordinates.') from None
        chain, insertion, name = line[21], line[26].strip(), line[17:20].strip()
        key = (chain, number, insertion)
        if key != active_key:
            if key in by_key:
                raise ValueError(f'Line {lineno}: repeated non-contiguous residue {chain!r}:{number}{insertion}; ambiguous mapping.')
            residue = Residue(chain, number, insertion, name)
            residues.append(residue)
            by_key[key] = residue
            active_key = key
        residue = by_key[key]
        if residue.name != name:
            raise ValueError(f'Line {lineno}: multiple residue names at one position; resolve microheterogeneity first.')
        residue.lines.append(line)
        atom_count += 1
    if not residues:
        raise ValueError('No ATOM residues found. This version supports PDB, not mmCIF.')
    if atom_count > 99999:
        raise ValueError('Too many atoms for standard PDB serial numbering.')
    if any(r.name not in AA for r in residues):
        warnings.append('Non-standard ATOM residue names are shown as X; residue names are preserved.')
    if any(line[16].strip() for r in residues for line in r.lines):
        warnings.append('Alternate atom locations are preserved; no conformer is selected automatically.')
    for chain in dict.fromkeys(r.chain for r in residues):
        rs = [r for r in residues if r.chain == chain]
        if any(b.number-a.number != 1 or a.insertion or b.insertion for a,b in zip(rs,rs[1:])):
            warnings.append(f'Chain {chain.strip() or "(blank)"}: numbering has gaps, insertion codes or unusual order. Review mapping; this alone is not an error.')
    warnings.append('Export contains selected ATOM coordinates and TER/END only; headers, bonds, ANISOU and other records are not retained. No missing atoms or residues are repaired.')
    return residues, warnings


def check_hotspots(residues, chain, text):
    available = {r.label: r for r in residues if r.chain == chain}
    tokens = [t.strip() for t in text.split(',') if t.strip()]
    result = []
    for token in tokens:
        valid = bool(re.fullmatch(r'-?\d+[A-Za-z]?', token))
        result.append({'residue':token, 'exists':valid and token in available,
                       'message':'Found' if valid and token in available else 'Not present' if valid else 'Use numbers with optional insertion code, e.g. 125A'})
    return result


def select_ranges(residues, chain, ranges):
    """Return (chain residues, sorted unique observed indices). Never infer IDs."""
    chain_rs = [r for r in residues if r.chain == chain]
    if not isinstance(ranges, list) or not ranges:
        raise ValueError('Choose at least one crop range.')
    indices = set()
    for bounds in ranges:
        if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
            raise ValueError('Each crop range needs a start and an end.')
        start, end = bounds
        if (type(start) is not int or type(end) is not int
                or not 0 <= start <= end < len(chain_rs)):
            raise ValueError('Invalid selection: choose ordered endpoints in the same chain.')
        indices.update(range(start, end+1))
    return chain_rs, sorted(indices)


def export_selection(residues, chain, ranges, numbering='sequential', start_number=1):
    """Export the union of inclusive observed ranges in original file order.

    Preserve original identities, or assign consecutive numbers from start_number.
    A TER separates disjoint observed segments; it does not assert peptide bonds.
    """
    chain_rs, indices = select_ranges(residues, chain, ranges)
    if numbering not in ('sequential', 'original'):
        raise ValueError('Numbering must be sequential or original.')
    if numbering == 'sequential' and (type(start_number) is not int
            or not -999 <= start_number <= start_number+len(indices)-1 <= 9999):
        raise ValueError('New residue numbers must fit PDB columns (-999 to 9999).')
    lines, rows, serial = [], [], 1
    previous = None
    for offset, index in enumerate(indices):
        r = chain_rs[index]
        if previous is not None and index != previous+1:
            lines.append('TER')
        new_number = r.number if numbering == 'original' else start_number+offset
        new_insertion = r.insertion if numbering == 'original' else ''
        rows.append({'original_chain':r.chain, 'original_resseq':r.number,
                     'original_icode':r.insertion, 'resname':r.name,
                     'original_observed_index':index+1,
                     'new_chain':r.chain, 'new_resseq':new_number, 'new_icode':new_insertion})
        for line in r.lines:
            lines.append(f'{line[:6]}{serial:5d}{line[11:22]}{new_number:4d}{new_insertion or " "}{line[27:]}')
            serial += 1
        previous = index
    lines += ['TER', 'END']
    output = io.StringIO(newline='')
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return '\n'.join(lines)+'\n', output.getvalue(), rows


def export_crop(residues, chain, start, end):
    """Backward-compatible single-range export, numbered from 1."""
    return export_selection(residues, chain, [[start, end]])


def mapping_segments(rows):
    """Split mapping rows at omitted observed residues, not numerical gaps."""
    segments = []
    for row in rows:
        if not segments or row['original_observed_index'] != segments[-1][-1]['original_observed_index']+1:
            segments.append([])
        segments[-1].append(row)
    return segments


def export_fasta(rows):
    """One observed ATOM-derived sequence record per disjoint selection segment."""
    lines = []
    for i, segment in enumerate(mapping_segments(rows), 1):
        first, last = segment[0], segment[-1]
        chain = first['original_chain'].strip() or 'blank'
        first_label = f"{first['original_resseq']}{first['original_icode']}"
        last_label = f"{last['original_resseq']}{last['original_icode']}"
        lines.append(f'>residue_map_chain_{chain}_segment_{i} original={first_label}:{last_label} observed_ATOM_only')
        sequence = ''.join(AA.get(row['resname'], 'X') for row in segment)
        lines.extend(sequence[j:j+60] for j in range(0, len(sequence), 60))
    return '\n'.join(lines)+'\n'


def inspect(text):
    residues, warnings = parse_pdb(text)
    chains = []
    for chain in dict.fromkeys(r.chain for r in residues):
        rs = [r for r in residues if r.chain == chain]
        chains.append({'id':chain, 'sequence':''.join(AA.get(r.name,'X') for r in rs),
                       'residues':[{'label':r.label,'name':r.name} for r in rs]})
    return {'chains':chains, 'warnings':warnings}

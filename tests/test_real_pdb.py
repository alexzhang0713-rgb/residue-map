"""Offline regression fixtures; provenance and hashes: fixtures/SOURCES.md."""
import csv
import hashlib
import io
from pathlib import Path
import unittest
from core import parse_pdb, export_crop, check_hotspots

FIXTURES = Path(__file__).parent / 'fixtures'

class RealPDBTests(unittest.TestCase):
    def test_fixture_hashes(self):
        for name, expected in {
            '1CRN': '42199a30a0701864a2a5cc76cd7f35cc544cd0e65fbcf63e03c166543249b811',
            '1IGT': '0a71448e9cb985a3562a4a086f807f86deb235a4f04bdccabf483d5a860b0c95',
        }.items():
            self.assertEqual(hashlib.sha256((FIXTURES / (name+'.pdb')).read_bytes()).hexdigest(), expected)

    def assert_roundtrip(self, residues, chain, start, end):
        chosen = [r for r in residues if r.chain == chain][start:end+1]
        pdb, table, rows = export_crop(residues, chain, start, end)
        reread, _ = parse_pdb(pdb)
        self.assertEqual(len(reread), len(chosen))
        self.assertEqual(len(list(csv.DictReader(io.StringIO(table)))), len(rows))
        original = {r.key: r for r in residues}
        for i, (new, row) in enumerate(zip(reread, rows), 1):
            old = original[(row['original_chain'], row['original_resseq'], row['original_icode'])]
            self.assertIs(old, chosen[i-1])
            self.assertEqual(new.key, (chain, i, ''))
            self.assertEqual(row['original_observed_index'], start+i)
            self.assertEqual(new.name, old.name)
            # Every atom survives; only serial and residue numbering may change.
            self.assertEqual([(l[:6], l[11:22], l[27:]) for l in new.lines],
                             [(l[:6], l[11:22], l[27:]) for l in old.lines])
        return rows

    def test_crambin_full_and_internal_crop(self):
        residues, _ = parse_pdb((FIXTURES/'1CRN.pdb').read_text())
        self.assertEqual(len(residues), 46)
        self.assertEqual(sum(len(r.lines) for r in residues), 327)
        self.assert_roundtrip(residues, 'A', 0, 45)
        rows = self.assert_roundtrip(residues, 'A', 9, 19)
        self.assertEqual([r['original_resseq'] for r in rows], list(range(10,21)))

    def test_immunoglobulin_all_chains_and_insertions(self):
        residues, warnings = parse_pdb((FIXTURES/'1IGT.pdb').read_text())
        self.assertTrue(any('HETATM' in w for w in warnings))
        for chain, count in [('A',214), ('B',444), ('C',214), ('D',444)]:
            self.assertEqual(sum(r.chain == chain for r in residues), count)
            self.assert_roundtrip(residues, chain, 0, count-1)
        chain = [r for r in residues if r.chain == 'B']
        start = next(i for i,r in enumerate(chain) if r.label == '82')
        end = next(i for i,r in enumerate(chain) if r.label == '83')
        rows = self.assert_roundtrip(residues, 'B', start, end)
        self.assertEqual([(r['original_resseq'],r['original_icode'],r['new_resseq']) for r in rows],
                         [(82,'',1),(82,'A',2),(82,'B',3),(82,'C',4),(83,'',5)])
        self.assertEqual([c['exists'] for c in check_hotspots(residues,'B','82,82A,82B,82C,82D,83')],
                         [True,True,True,True,False,True])

import csv
import io
from pathlib import Path
import unittest
from core import parse_pdb, export_selection, export_fasta, mapping_segments
from test_core import atom

class SelectionTests(unittest.TestCase):
    def setUp(self):
        self.residues,_ = parse_pdb('\n'.join([atom(1,123,name='GLY'),atom(2,125),
            atom(3,125,'A',name='SER'),atom(4,150,name='LYS')]))

    def test_disjoint_segments_keep_observed_indices_and_atoms(self):
        pdb,table,rows=export_selection(self.residues,'A',[[2,3],[0,0]],start_number=10)
        self.assertEqual([r['original_observed_index'] for r in rows],[1,3,4])
        self.assertEqual([r['new_resseq'] for r in rows],[10,11,12])
        self.assertEqual(pdb.splitlines().count('TER'),2)
        reread,_=parse_pdb(pdb)
        self.assertEqual([r.number for r in reread],[10,11,12])
        for new,old in zip(reread,[self.residues[i] for i in [0,2,3]]):
            self.assertEqual([l[27:] for l in new.lines],[l[27:] for l in old.lines])
        self.assertEqual(len(list(csv.DictReader(io.StringIO(table)))),3)
        self.assertEqual([len(seg) for seg in mapping_segments(rows)],[1,2])
        fasta=export_fasta(rows)
        self.assertEqual([l for l in fasta.splitlines() if not l.startswith('>')],['G','SK'])

    def test_overlap_duplicates_and_adjacent_ranges_merge(self):
        pdb,_,rows=export_selection(self.residues,'A',[[2,3],[0,1],[1,2],[2,3]])
        self.assertEqual(len(rows),4)
        self.assertEqual(pdb.splitlines().count('TER'),1)
        self.assertEqual(len(mapping_segments(rows)),1)
        self.assertEqual(export_fasta(rows).splitlines()[1],'GASK')

    def test_original_ids_including_insertion_are_preserved(self):
        pdb,_,rows=export_selection(self.residues,'A',[[1,3]],numbering='original')
        reread,_=parse_pdb(pdb)
        self.assertEqual([r.key for r in reread],[r.key for r in self.residues[1:]])
        self.assertEqual([(r['new_resseq'],r['new_icode']) for r in rows],[(125,''),(125,'A'),(150,'')])

    def test_blank_chain_negative_and_zero_numbering(self):
        residues,_=parse_pdb(atom(1,-1,chain=' ')+'\n'+atom(2,0,chain=' '))
        for mode,start,expected in [('original',1,[-1,0]),('sequential',0,[0,1]),('sequential',-999,[-999,-998])]:
            with self.subTest(mode=mode,start=start):
                pdb,_,rows=export_selection(residues,' ',[[0,1]],mode,start)
                reread,_=parse_pdb(pdb)
                self.assertEqual([r.number for r in reread],expected)
                self.assertTrue(all(r.chain==' ' for r in reread))
                self.assertIn('chain_blank',export_fasta(rows))

    def test_ranges_reject_empty_reversed_outside_and_noninteger(self):
        for ranges in [[],None,{},[[2,1]],[[0,4]],[[0]],[[0,1,2]],[[0.0,1]],[[False,1]],[[0,'1']]]:
            with self.subTest(ranges=ranges),self.assertRaises(ValueError):
                export_selection(self.residues,'A',ranges)

    def test_numbering_bounds_and_mode(self):
        for mode,start in [('bad',1),('sequential',9999),('sequential',-1000),('sequential',1.5),('sequential',True),('sequential',None)]:
            with self.subTest(mode=mode,start=start),self.assertRaises(ValueError):
                export_selection(self.residues,'A',[[0,1]],mode,start)
        pdb,_,_=export_selection(self.residues,'A',[[0,0]],start_number=9999)
        self.assertEqual(parse_pdb(pdb)[0][0].number,9999)

    def test_fasta_unknown_and_wrapping(self):
        residues,_=parse_pdb('\n'.join(atom(i,i,name='UNK' if i==1 else 'ALA') for i in range(1,63)))
        _,_,rows=export_selection(residues,'A',[[0,61]])
        fasta=export_fasta(rows).splitlines()
        self.assertEqual(len(fasta[1]),60)
        self.assertEqual(fasta[1][0],'X')
        self.assertEqual(fasta[2],'AA')

    def test_real_igt_disjoint_original_mapping_and_fasta(self):
        residues,_=parse_pdb((Path(__file__).parent/'fixtures/1IGT.pdb').read_text())
        chain=[r for r in residues if r.chain=='B']
        indices={r.label:i for i,r in enumerate(chain)}
        ranges=[[indices['82A'],indices['83']],[indices['52A'],indices['52A']]]
        pdb,_,rows=export_selection(residues,'B',ranges,'original')
        selected=[chain[indices['52A']]]+chain[indices['82A']:indices['83']+1]
        reread,_=parse_pdb(pdb)
        self.assertEqual([r.key for r in reread],[r.key for r in selected])
        for new,old in zip(reread,selected):
            self.assertEqual([l[11:] for l in new.lines],[l[11:] for l in old.lines])
        self.assertEqual([l for l in export_fasta(rows).splitlines() if not l.startswith('>')],['N','SRLK'])

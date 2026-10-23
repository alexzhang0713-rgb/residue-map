import unittest
from core import parse_pdb, export_crop, check_hotspots, inspect


def atom(serial, number, insertion='', chain='A', name='ALA', alt='', x=1.0):
    return f'ATOM  {serial:5d}  CA {alt:1}{name:3} {chain:1}{number:4d}{insertion:1}   {x:8.3f}{2.:8.3f}{3.:8.3f}{1.:6.2f}{20.:6.2f}           C  '


class MappingTests(unittest.TestCase):
    def test_insertion_gap_roundtrip(self):
        text='\n'.join([atom(1,125),atom(2,125,'A'),atom(3,150)])
        residues,warnings=parse_pdb(text)
        pdb,csv,rows=export_crop(residues,'A',1,2)
        reread,_=parse_pdb(pdb)
        self.assertEqual([(r.number,r.insertion) for r in reread],[(1,''),(2,'')])
        self.assertEqual([(r['original_resseq'],r['original_icode']) for r in rows],[(125,'A'),(150,'')])
        self.assertEqual(residues[1].lines[0][30:],reread[0].lines[0][30:])
        self.assertIn('125,A',csv)
        self.assertTrue(any('gaps' in w for w in warnings))

    def test_chain_and_blank_and_negative(self):
        residues,_=parse_pdb('\n'.join([atom(1,-1,chain=' '),atom(2,-1,chain='B')]))
        self.assertTrue(check_hotspots(residues,' ','-1')[0]['exists'])
        pdb,_,_=export_crop(residues,' ',0,0)
        self.assertEqual(parse_pdb(pdb)[0][0].chain,' ')
        self.assertEqual(len(inspect('\n'.join(r.lines[0] for r in residues))['chains']),2)

    def test_alternates_preserved(self):
        residues,warnings=parse_pdb(atom(1,1,alt='A')+'\n'+atom(2,1,alt='B'))
        pdb,_,_=export_crop(residues,'A',0,0)
        reread,_=parse_pdb(pdb)
        self.assertEqual([line[16] for line in reread[0].lines],['A','B'])
        self.assertTrue(any('Alternate' in w for w in warnings))

    def test_hotspots_exact_identity(self):
        residues,_=parse_pdb(atom(1,125,'A'))
        checks=check_hotspots(residues,'A','125,125A,999,125-130')
        self.assertEqual([c['exists'] for c in checks],[False,True,False,False])

    def test_multimodel_rejected(self):
        with self.assertRaisesRegex(ValueError,'Multiple MODEL'):
            parse_pdb('MODEL        1\n'+atom(1,1)+'\nENDMDL\nMODEL        2\n'+atom(1,1))

    def test_noncontiguous_duplicate_rejected(self):
        with self.assertRaisesRegex(ValueError,'non-contiguous'):
            parse_pdb('\n'.join([atom(1,1),atom(2,2),atom(3,1)]))

    def test_truncated_invalid_and_empty(self):
        for text in ['ATOM      1',atom(1,1,x=float('nan')),'HEADER\nEND']:
            with self.assertRaises(ValueError):
                parse_pdb(text)

    def test_bad_crop_rejected(self):
        residues,_=parse_pdb(atom(1,1))
        for chain,start,end in [('A',1,0),('B',0,0),('A',-1,0)]:
            with self.assertRaises(ValueError):
                export_crop(residues,chain,start,end)

    def test_hetatm_excluded_and_reported(self):
        residues,warnings=parse_pdb(atom(1,1)+'\n'+atom(2,2).replace('ATOM  ','HETATM',1))
        self.assertEqual(len(residues),1)
        self.assertTrue(any('HETATM' in w for w in warnings))

    def test_single_model_allowed(self):
        residues,_=parse_pdb('MODEL        1\n'+atom(1,1)+'\nENDMDL')
        self.assertEqual(len(residues),1)

    def test_microheterogeneity_rejected(self):
        with self.assertRaisesRegex(ValueError,'multiple residue names'):
            parse_pdb(atom(1,1,name='ALA')+'\n'+atom(2,1,name='GLY'))

if __name__=='__main__':
    unittest.main()

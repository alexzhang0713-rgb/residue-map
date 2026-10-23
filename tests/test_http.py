import json
from pathlib import Path
import threading
import unittest
import urllib.request
import urllib.error
from http.server import ThreadingHTTPServer
from server import Handler

class HTTPTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
        cls.thread=threading.Thread(target=cls.server.serve_forever,daemon=True)
        cls.thread.start()
        cls.base=f'http://127.0.0.1:{cls.server.server_port}'
        cls.pdb=(Path(__file__).parent.parent/'examples/numbering_demo.pdb').read_text()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join()

    def post(self,path,payload,origin=None):
        headers={'Content-Type':'application/json'}
        if origin:headers['Origin']=origin
        req=urllib.request.Request(self.base+'/api/'+path,data=json.dumps(payload).encode(),headers=headers)
        with urllib.request.urlopen(req) as response:return json.load(response)

    def test_full_flow(self):
        self.assertEqual(self.post('inspect',{'pdb':self.pdb})['chains'][0]['sequence'],'GASK')
        payload={'pdb':self.pdb,'chain':'A','start':2,'end':3,'hotspots':'125A,150'}
        result=self.post('export',payload)
        self.assertEqual([m['new_resseq'] for m in result['mapping']],[1,2])
        self.assertEqual(result['report']['hotspots'][0]['new'],1)
        self.assertIn('125,A',result['csv'])
        payload['hotspots']='123,150'
        self.assertEqual(self.post('export',payload)['report']['hotspots'][0]['status'],'outside crop')
        payload['hotspots']='999'
        with self.assertRaises(urllib.error.HTTPError) as caught:
            self.post('export',payload)
        self.assertEqual(caught.exception.code,400)

    def test_multirange_original_ids_and_fasta(self):
        payload={'pdb':self.pdb,'chain':'A','ranges':[[2,3],[0,0]],
                 'numbering':'original','hotspots':'123,125,125A'}
        result=self.post('export',payload)
        self.assertEqual(result['report']['format_version'],2)
        self.assertEqual(result['report']['numbering'],{'mode':'original','start_number':None})
        self.assertEqual(len(result['report']['segments']),2)
        self.assertEqual([h['new_label'] for h in result['report']['hotspots']],['123',None,'125A'])
        self.assertEqual(result['report']['hotspots'][1]['status'],'outside crop')
        self.assertEqual(result['pdb'].splitlines().count('TER'),2)
        self.assertEqual([l for l in result['fasta'].splitlines() if not l.startswith('>')],['G','SK'])
        self.assertTrue(any('Disjoint' in w for w in result['report']['warnings']))

    def test_custom_numbering_zero_and_invalid_options(self):
        payload={'pdb':self.pdb,'chain':'A','ranges':[[2,3]],'hotspots':'125A',
                 'numbering':'sequential','start_number':0}
        result=self.post('export',payload)
        self.assertEqual(result['report']['hotspots'][0]['new_label'],'0')
        self.assertEqual(result['mapping'][1]['new_resseq'],1)
        for edit in [{'start_number':9999},{'start_number':1.5},{'start_number':None},
                     {'numbering':'bad'},{'ranges':[]},{'ranges':[[3,2]]}]:
            with self.subTest(edit=edit),self.assertRaises(urllib.error.HTTPError) as caught:
                self.post('export',dict(payload,**edit))
            self.assertEqual(caught.exception.code,400)

    def test_cross_origin_rejected(self):
        with self.assertRaises(urllib.error.HTTPError) as caught:
            self.post('inspect',{'pdb':self.pdb},'https://example.com')
        self.assertEqual(caught.exception.code,400)

    def test_static_routes(self):
        with urllib.request.urlopen(self.base) as r:
            self.assertIn(b'Residue Map',r.read())
        with urllib.request.urlopen(self.base+'/example') as r:
            self.assertIn(b'SYNTHETIC',r.read())

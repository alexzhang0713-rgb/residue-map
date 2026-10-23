"""Local-only UI server. No external services, uploads or dependencies."""
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import argparse
import hashlib
import json
import webbrowser
from core import inspect, parse_pdb, check_hotspots, export_selection, export_fasta, mapping_segments

ROOT = Path(__file__).parent
LIMIT = 10 * 1024 * 1024

class Handler(BaseHTTPRequestHandler):
    def log_message(self, *_):
        pass

    def do_GET(self):
        path = self.path.split('?')[0]
        if path == '/':
            data, mime = (ROOT/'index.html').read_bytes(), 'text/html; charset=utf-8'
        elif path == '/example':
            data, mime = (ROOT/'examples/numbering_demo.pdb').read_bytes(), 'text/plain'
        else:
            self.send_error(404)
            return
        self.send_response(200)
        self.send_header('Content-Type',mime)
        self.send_header('Content-Length',str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_POST(self):
        try:
            # Same-origin browser requests only; reject cross-origin requests to localhost.
            origin = self.headers.get('Origin')
            if origin and origin != f'http://{self.headers.get("Host")}':
                raise ValueError('Cross-origin requests are not allowed.')
            length = int(self.headers.get('Content-Length','0'))
            if not 0 < length <= LIMIT:
                raise ValueError('Request must be under 10 MB.')
            payload = json.loads(self.rfile.read(length))
            text = payload['pdb']
            if not isinstance(text,str):
                raise ValueError('PDB must be text.')
            if self.path == '/api/inspect':
                result = inspect(text)
            elif self.path == '/api/export':
                residues, warnings = parse_pdb(text)
                chain = payload['chain']
                checks = check_hotspots(residues,chain,payload.get('hotspots',''))
                if any(not c['exists'] for c in checks):
                    raise ValueError('Some hotspots do not exist. Correct them or clear the hotspot input.')
                ranges = payload.get('ranges', [[payload.get('start'), payload.get('end')]])
                numbering = payload.get('numbering', 'sequential')
                start_number = payload.get('start_number', 1)
                pdb, csv, rows = export_selection(residues,chain,ranges,numbering,start_number)
                segments = mapping_segments(rows)
                if len(segments) > 1:
                    warnings.append('Disjoint observed segments are separated by TER; FASTA has one record per segment. No peptide connectivity is inferred.')
                mapped = {f"{r['original_resseq']}{r['original_icode']}":r for r in rows}
                hotspot_mapping = []
                for check in checks:
                    row = mapped.get(check['residue'])
                    hotspot_mapping.append({'original':check['residue'],
                        'new':row['new_resseq'] if row else None,
                        'new_icode':row['new_icode'] if row else None,
                        'new_label':f"{row['new_resseq']}{row['new_icode']}" if row else None,
                        'status':'retained' if row else 'outside crop'})
                result = {'pdb':pdb,'csv':csv,'fasta':export_fasta(rows),'mapping':rows,
                          'report':{'format_version':2,'input_sha256':hashlib.sha256(text.encode()).hexdigest(),
                                    'selected_chain':chain,'warnings':warnings,
                                    'numbering':{'mode':numbering,'start_number':start_number if numbering == 'sequential' else None},
                                    'segments':[{'start_observed_index':seg[0]['original_observed_index'],
                                                 'end_observed_index':seg[-1]['original_observed_index']} for seg in segments],
                                    'fasta_scope':'One observed ATOM-only record per disjoint selected segment; unknown residue names are X.',
                                    'hotspots':hotspot_mapping,'mapping':rows}}
            elif self.path == '/api/check':
                residues,_ = parse_pdb(text)
                result = check_hotspots(residues,payload['chain'],payload.get('hotspots',''))
            else:
                self.send_error(404)
                return
            self.respond(result,200)
        except (ValueError, KeyError, TypeError) as exc:
            self.respond({'error':str(exc)},400)

    def respond(self,obj,status):
        data = json.dumps(obj).encode()
        self.send_response(status)
        self.send_header('Content-Type','application/json')
        self.send_header('Content-Length',str(len(data)))
        self.end_headers()
        self.wfile.write(data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--port',type=int,default=8765)
    parser.add_argument('--no-browser',action='store_true')
    args = parser.parse_args()
    server = ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    url = f'http://127.0.0.1:{args.port}'
    print(f'Residue Map running at {url}. Press Ctrl+C to stop.',flush=True)
    if not args.no_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        server.server_close()

// Pure JavaScript state regression; no browser or external packages required.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
class Element {
  constructor(){this.hidden=true;this.value='';this.children=[];this.listeners={};this.files=[];}
  replaceChildren(...nodes){this.children=nodes;}
  append(...nodes){this.children.push(...nodes);}
  add(node){this.children.push(node);if(this.children.length===1)this.value=node.value;}
  addEventListener(name, fn){this.listeners[name]=fn;}
  scrollIntoView(){}
}
const elements = new Map();
const get = id => {if(!elements.has(id))elements.set(id,new Element());return elements.get(id);};
const pending=[],requests=[];
const source=fs.readFileSync(require('node:path').join(__dirname,'../index.html'),'utf8').split('<script>')[1].split('</script>')[0];
const context=vm.createContext({document:{getElementById:get,createElement:()=>new Element()},
  Option:function(text,value){this.textContent=text;this.value=String(value);},
  fetch:(url,options)=>{requests.push({url,body:JSON.parse(options.body)});return new Promise(resolve=>pending.push(body=>resolve({ok:true,json:async()=>body})));},
  URL:{revokeObjectURL(){},createObjectURL(){return 'blob:test';}},Blob,console});
vm.runInContext(source,context);
const run = code => vm.runInContext(code,context);
const inspection={chains:[{id:'A',sequence:'AS',residues:[{label:'1',name:'ALA'},{label:'2',name:'SER'}]}],warnings:[]};
async function loadFixture(){const task=run("load('ATOM','fixture',beginLoad())");pending.shift()(inspection);await task;}
(async()=>{
 await loadFixture();
 get('numbering').value='sequential';get('newstart').value='1';
 const check=get('check').onclick();
 get('hotspots').value='2';get('hotspots').listeners.input();
 pending.shift()([{residue:'1',exists:true,message:'Found'}]);await check;
 assert.equal(get('checks').children.length,0,'old check must not reappear after edit');
 const exporting=get('export').onclick();get('start').listeners.input();
 pending.shift()({});await exporting;
 assert.equal(get('result').hidden,true,'old export must not reappear after edit');
 assert.equal(get('export').disabled,false);
 const first=run("load('first','first',beginLoad())");
 const second=run("load('second','second',beginLoad())");
 const replyFirst=pending.shift(),replySecond=pending.shift();
 replySecond(inspection);await second;replyFirst({});await first;
 assert.equal(get('filename').textContent,'  second','latest import wins');
 // New controls participate in request construction and async invalidation.
 get('add-range').onclick();
 assert.equal(run('extraRanges.length'),1);
 const range=run('extraRanges[0]');range.first.value='0';range.last.value='0';
 get('newstart').value='10';
 const rangeExport=get('export').onclick();
 assert.deepEqual(requests.at(-1).body.ranges,[[0,1],[0,0]]);
 assert.equal(requests.at(-1).body.start_number,10);
 range.first.listeners.input();pending.shift()({});await rangeExport;
 assert.equal(get('result').hidden,true);
 range.remove.onclick();assert.equal(run('extraRanges.length'),0);
 get('numbering').value='original';get('numbering').onchange();
 assert.equal(get('newstart').disabled,true);
 const originalExport=get('export').onclick();
 assert.equal(requests.at(-1).body.numbering,'original');
 pending.shift()({pdb:'ATOM',csv:'mapping',fasta:'>test\nA\n',mapping:[{
   original_chain:'A',original_resseq:125,original_icode:'A',resname:'ALA',new_resseq:125,new_icode:'A'}],
   report:{hotspots:[{original:'125A',new:125,new_icode:'A',new_label:'125A'}],
     segments:[{}],numbering:{mode:'original',start_number:null}}});
 await originalExport;
 assert.equal(get('downloads').children.length,4);
 assert.equal(get('mappedhotspots').children[0].textContent,'125A → 125A');
 get('newstart').listeners.input();assert.equal(get('result').hidden,true);
 get('add-range').onclick();run('chainChanged()');assert.equal(run('extraRanges.length'),0);
 get('result').hidden=false;get('downloads').append({});
 get('file').files=[{size:8*1024*1024+1}];await get('file').onchange({target:get('file')});
 assert.equal(get('workspace').hidden,true,'oversized upload clears old structure');
 assert.equal(get('downloads').children.length,0);
 assert.equal(get('error').hidden,false);
 console.log('UI state regressions passed: stale responses, multi-range controls, numbering, four downloads, oversized upload.');
})().catch(e=>{console.error(e);process.exitCode=1;});

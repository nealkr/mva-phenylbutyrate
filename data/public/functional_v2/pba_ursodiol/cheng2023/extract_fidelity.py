from pathlib import Path
import csv,json,hashlib,zipfile,xml.etree.ElementTree as ET
base=Path(__file__).resolve().parent
ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
plans=[('42255_2023_832_MOESM4_ESM.xlsx','Fig.2h','xl/worksheets/sheet7.xml',[('misalignment','ABC',4,7),('lagging_chromosomes','EFG',4,7),('micronuclei','IJK',4,6)]),('42255_2023_832_MOESM6_ESM.xlsx','Fig.3n','xl/worksheets/sheet9.xml',[('misalignment','ABCD',4,6),('lagging_chromosomes','FGHI',4,6),('micronuclei','KLMN',4,6)])]
obs=[]; cells=[]; inventory=[]
for fn,sheet,part,blocks in plans:
 b=(base/fn).read_bytes()
 with zipfile.ZipFile(base/fn) as z:
  strings=[''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si',ns)]
  root=ET.fromstring(z.read(part)); d={}
  for c in root.findall('.//m:sheetData/m:row/m:c',ns):
   v=c.find('m:v',ns)
   if v is None:continue
   val=strings[int(v.text)] if c.get('t')=='s' else v.text
   d[c.get('r')]=val
   cells.append(dict(file=fn,sheet=sheet,xml_part=part,cell=c.get('r'),value=val))
  for endpoint,cols,lo,hi in blocks:
   for col in cols:
    arm=d[col+'3']
    for row in range(lo,hi+1):
     addr=col+str(row)
     obs.append(dict(file=fn,sheet=sheet,cell=addr,endpoint=endpoint,arm=arm,source_row=row,within_arm_row_index=row-lo+1,value_percent=float(d[addr]),biological_unit='independent_experiment_per_figure_caption',original_experiment_id='',numerator='',denominator='',pairing_established=False))
  inventory.append(dict(file=fn,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),zip_test=z.testzip(),comments_present=any('comment' in p.lower() for p in z.namelist()),source_url='https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42255-023-00832-9/MediaObjects/'+fn,metadata_md5_match=None))
for name,data in [('fidelity_observations.csv',obs),('fidelity_source_cells.csv',cells)]:
 with (base/name).open('w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
(base/'source_workbook_inventory.json').write_text(json.dumps(inventory,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'observations':len(obs),'source_cells':len(cells),'workbooks':inventory},indent=2))
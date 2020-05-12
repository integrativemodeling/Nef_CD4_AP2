import ihm.dictionary
import urllib.request

d = ihm.dictionary.Dictionary()
for dn in 'mmcif_pdbx_v50.dic', 'mmcif_ihm.dic':
     url = 'http://mmcif.wwpdb.org/dictionaries/ascii/' + dn
     with urllib.request.urlopen(url) as fh:
         d += ihm.dictionary.read(fh)

with open('IM_Nef-CD4-AP2.cif') as fh:
     d.validate(fh)

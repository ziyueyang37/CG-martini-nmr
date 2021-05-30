import MDAnalysis as md
import os
#import subprocess
#from subprocess import Popen


CGDict = {}
CGDict['VAL'] = ['P5','AC2']
CGDict['ILE'] = ['P5','AC1']
CGDict['ASN'] = ['P5','P5']
CGDict['THR'] = ['P5','P1']
CGDict['PHE'] = ['P5','SC5','SC5','SC5']
CGDict['ASP'] = ['P5','Qa'] 
CGDict['GLY'] = ['P5']
CGDict['ALA'] = ['P4']
CGDict['TYR'] = ['P5','SC4','SC4','SP1']
CGDict['LEU'] = ['P5','AC1'] 
CGDict['GLN'] = ['P5','P4']
CGDict['HIS'] = ['P5','SC4','SP1','SP1']
CGDict['LYS'] = ['P5','C3','Qd']
CGDict['PRO'] = ['P4','C3']
CGDict['ASP'] = ['P5','Qa']
CGDict['SER'] = ['P5','P1']
CGDict['GLU'] = ['P5','Qa'] 
CGDict['TRP'] = ['P5','SC4','SNd','SC5','SC5']
CGDict['ARG'] = ['P5','N0','Qd']
CGDict['THR'] = ['P5','P1','P5','P1']
CGDict['MET'] = ['P5','C5']
CGDict['CYS'] = ['P5','C5']

def assign_beads(CGDict, path, file_name, new_path):
    CGpdb = path + '/' + file_name
    print(CGpdb)
    u = md.Universe(CGpdb)
    previous = u.atoms[-1]
    count = 0
    for atom in u.atoms:
        residue = atom.resname
        resid = atom.resid
        if residue == previous.resname and resid == previous.resid:
            count += 1
        else: count = 0
        #print(residue)
        print(CGDict[residue][count])
        atom.name = CGDict[residue][count]
        previous = atom
    protein = u.select_atoms("protein")
    protein.write(new_path + '/' + file_name)
    return 0
    

g = os.walk('../mtn')
new_path = '../CGpdb'
for path, dir_list, file_list in g:
    for file_name in file_list:
        try:
            assign_beads(CGDict, path, file_name, new_path)
        except:
            pass
        #print(file_name)
        #new_file_name = file_name[5:9]
        #MAX_ATOM = os.popen("python3 camshift.py %s" % path+'/'+file_name).read()
        #os.system("\cp /home/zyang43/nmr/test-dir/%s /home/zyang43/nmr/data/template.pdb" % file_name)
        #os.system("echo -e 'peptide: GROUP ATOMS=1-%d\nWHOLEMOLECULES ENTITY0=peptide\ncs: CS2BACKBONE ATOMS=peptide DATADIR=data/\nPRINT ARG=(cs\.ha-.*),(cs\.hn-.*) FILE=%s STRIDE=1' > plumed.dat" % (int(MAX_ATOM), new_path+new_file_name))        os.system("plumed driver --plumed plumed.dat --mf_pdb %s" % path+'/'+file_name)


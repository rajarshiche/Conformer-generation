# -*- coding: utf-8 -*-
"""
@author: rguha

"""

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
import time

start_time = time.time()
 

def confgen(input, output, prunermsthresh, numconf):

    
    mol = Chem.AddHs(Chem.MolFromSmiles(input), addCoords=True)
    param = rdDistGeom.ETKDGv2()
    param.pruneRmsThresh = prunermsthresh
    
    
    cids = AllChem.EmbedMultipleConfs(mol, numConfs = numconf, numThreads=0,useRandomCoords=True,
                           maxAttempts = 10000, pruneRmsThresh = prunermsthresh ) #can allow user to change numConfs
    
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant='MMFF94s')
    
    w = Chem.SDWriter(output)
    
    
    res = []
 
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid, e))
    sorted_res = sorted(res, key=lambda x:x[1])
    
    rdMolAlign.AlignMolConformers(mol)
    
    
    i=0
    for cid, e in sorted_res:
        
        print(e)
        mol.SetProp('CID', str(cid))
        mol.SetProp('Energy', str(e))
        
        if i == 0: ## Saving only the lowest energy conformer
            w.write(mol, confId=cid)
            w.close()
        i = i + 1    
    
 
if __name__=='__main__':
    
    
    smiles = 'CC1(C)[C@H]2CC[C@@]1(C)C(NC(=O)C(c1c(C(=O)O)[nH]c3cc(Cl)ccc13)N(C=O)Cc1ccc(OCc3ccccc3Br)cc1)C2'
    smiles = 'CC1(C)[C@H]2CC[C@@]1(C)C(NC(=O)C(c1c(C(=O)O)[nH]c3cc(Cl)ccc13)N(C=O)Cc1ccc(OCc3ccc(Br)cc3)cc1)C2'

    ### Generating the conformers (10 numbers), calculating the energies and saving the minimum energy conformer as .sdf
    confgen( smiles,  'gen_confs_MMFF94s.sdf', 0.5,  10) 
    
    
# ###----------------------------------Total Process Time-------------------------------------------------------------------------------------------------------
end_time = time.time()
print("The duration of run (in s): ", end_time - start_time)    
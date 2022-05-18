import os
from rdkit import Chem
#from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
#from rdkit.Chem import RDConfig
#from rdkit import rdBase
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolAlign
#from rdkit.Chem import TemplateAlign
#from rdkit.Chem import rdDepictor
#rdDepictor.SetPreferCoordGen(True)

try:
    from openbabel import openbabel
except:
    print(">>> Warning:\n"
          "            Could not find openbabel!!! SMILES2MOL2 is not available!\n"
          ">>> To install openbabel:\n"
          "            conda install -c conda-forge openbabel")

def smiles2sdf(smiles, labels=None, ref=None, mergesdf=False, onlysdf=True):
    if isinstance(smiles, str):
        mols = []
        mols.append(smiles)
        if labels:
            label = []
            label.append(labels)
        else:
            mlabel = []
        mergesdf=False
    elif isinstance(smiles, list):
        mols = smiles
        if labels:
            label = labels
        else:
            mlabel = []
        if mergesdf:
            w0 = Chem.SDWriter('all_mols.sdf')
    for i,smile in enumerate(mols):
        try:
            mol = Chem.MolFromSmiles(smile)
            if labels:
                mol.SetProp("_Name",label[i])
            mol = Chem.AddHs(mol)
            ps = AllChem.ETKDG()
            ps.randomSeed = 0xf00d
            Chem.AllChem.EmbedMolecule(mol,ps)
            if ref:
                suppl = Chem.SDMolSupplier(ref)
                refmol = suppl[0]
                refmol = Chem.AddHs(refmol)
                o3d = rdMolAlign.GetO3A(mol,refmol)
                try:
                    o3d.Align()
                except:
                    pass
            if labels:
                w = Chem.SDWriter('%s.sdf'%(label[i]))
            else:
                l = "mol-%d"%(i+1)
                mlabel.append(l)
                w = Chem.SDWriter('%s.sdf'%(l))
            if mergesdf:
                w0.write(mol)
            else:
                w.write(mol)
        except:
            if labels:
                print("Failed to prepare molecule : %s with SMILES = %s" % (label[i],smile))
            else:
                print("Failed to prepare molecule : mol-%d with SMILES = %s"%(i+1,smile))
    if not onlysdf:
        if labels:
            return mols, label
        else:
            return mols, mlabel

def smiles2mol2(smiles, labels=None, ref=None):
    mergesdf=False
    mols, label = smiles2sdf(smiles, labels=labels, ref=ref, mergesdf=mergesdf, onlysdf=False)
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "mol2")
    for l in label:
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, "%s.sdf"%(l))
        obConversion.WriteFile(mol, "%s.mol2"%(l))
        os.system("rm %s.sdf"%(l))
    return mols, label

def smiles2pdbqt(smiles, labels=None, ref=None):
    mols, label = smiles2mol2(smiles, labels=labels, ref=ref)
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "pdbqt")
    for l in label:
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, "%s.mol2"%(l))
        obConversion.WriteFile(mol, "%s.pdbqt"%(l))
        os.system("rm %s.mol2"%(l))
    return mols, label

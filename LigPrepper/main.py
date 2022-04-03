from rdkit import Chem
#from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
#from rdkit.Chem import RDConfig
#from rdkit import rdBase
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolAlign
#from rdkit.Chem import TemplateAlign
#from rdkit.Chem import rdDepictor

def smiles2sdf(smiles, labels=None, ref=None, mergesdf=False):
    if isinstance(smiles, str):
        mols = []
        mols.append(smiles)
        if labels:
            label = []
            label.append(labels)
        else:
            label=None
        mergesdf=False
    elif isinstance(smiles, list):
        mols = smiles
        if labels:
            label = labels
        if mergesdf:
            w0 = Chem.SDWriter('all_mols.sdf')
    for i,smile in enumerate(mols):
        mol = Chem.MolFromSmiles(smile)
        if label:
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
        if label:
            w = Chem.SDWriter('%s.sdf'%(label[i]))
        else:
            w = Chem.SDWriter('mol-%d.sdf'%(i+1))
        if mergesdf:
            w0.write(mol)
        else:
            w.write(mol)

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
          "            Could not find OpenBabel!!! SMILES2MOL2 and SMILES2PDBQT are not available!\n"
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

def sdf2pdbqt(sdfile):
    l=str(sdfile)[0:-4]
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, "%s.sdf"%(l))
    obConversion.WriteFile(mol, "%s.mol2"%(l))
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("mol2", "pdbqt")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, "%s.mol2"%(l))
    obConversion.WriteFile(mol, "%s.pdbqt"%(l))
    os.system("rm %s.mol2"%(l))
    return "%s.pdbqt"%(l)

def pdbqt2sdf(ligand):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt","sdf")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, "%s.pdbqt"%(str(ligand)[0:-6]))
    obConversion.WriteFile(mol, "%s.sdf"%(str(ligand)[0:-6]))
    return "%s.sdf"%(str(ligand)[0:-6])

def splitsdf(sdfile,outputdir=None):
    print("Splitting %s"%(sdfile), end=" : ")
    f = open(sdfile, 'r')
    count=0
    molsep="$$$$"
    for line in f.readlines():
        l = str(line).split()
        if molsep in l:
            count += 1
    f.close()
    sdfiles = []
    if count == 1:
        sdfiles.append(sdfile)
        return count, sdfiles
    elif count > 1:
        f = open(sdfile, 'r')
        molnameline=True
        for line in f.readlines():
            l = str(line).split()
            if not molsep in l:
                if molnameline:
                    molname=str(l[0])
                    if outputdir:
                        sdfilename=outputdir+'/'+molname+".sdf"
                    else:
                        sdfilename=molname+".sdf"
                    newf = open(sdfilename,'w')
                    newf.write(line)
                    molnameline=False
                    sdfiles.append(sdfilename)
                else:
                    newf.write(line)
            else:
                newf.write(line)
                newf.close()
                molnameline=True
        return count, sdfiles
    else:
        print(" *** Error: sdfile could not be processed! *** ")
        exit()
    print("Done!")

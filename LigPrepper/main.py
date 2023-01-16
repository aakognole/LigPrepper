import os
from rdkit import Chem
from rdkit.Chem import Draw
#from rdkit.Chem import DataStructs
from rdkit.Chem import AllChem
#from rdkit.Chem import RDConfig
#from rdkit import rdBase
#from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdMolAlign
#from rdkit.Chem import TemplateAlign
#from rdkit.Chem import rdDepictor
from rdkit.Chem import rdFMCS
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

def splitsdf(sdfile,outputdir=None,parts=0,molspersdf=1):
    if parts > 0 and molspersdf > 1:
        parts = 0
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
    elif count > 1 and parts == 0 and molspersdf == 1:
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
    elif count > 1 and parts > 0 and molspersdf == 1:
        f = open(sdfile, 'r')
        bunch=int(count/parts)
        newbunch = True
        molnum, part = 0, 0
        for line in f.readlines():
            l = str(line).split()
            if molsep in l:
                molnum += 1
                newf.write(line)
                if molnum >= bunch:
                    newf.close()
                    newbunch = True
                    molnum = 0
                else:
                    newbunch = False
            else:
                if newbunch:
                    part += 1
                    filename=str(sdfile)+'_part'+str(part)
                    if outputdir:
                        sdfilename=outputdir+'/'+filename+".sdf"
                    else:
                        sdfilename=filename+".sdf"
                    newf = open(sdfilename,'w')
                    sdfiles.append(sdfilename)
                    newbunch = False
                newf.write(line)
    elif count > 1 and parts == 0 and molspersdf > 1:
        f = open(sdfile, 'r')
        bunch=int(molspersdf)
        newbunch = True
        molnum, part = 0, 0
        for line in f.readlines():
            l = str(line).split()
            if molsep in l:
                molnum += 1
                newf.write(line)
                if molnum >= bunch:
                    newf.close()
                    newbunch = True
                    molnum = 0
                else:
                    newbunch = False
            else:
                if newbunch:
                    part += 1
                    filename=str(sdfile)+'_part'+str(part)
                    if outputdir:
                        sdfilename=outputdir+'/'+filename+".sdf"
                    else:
                        sdfilename=filename+".sdf"
                    newf = open(sdfilename,'w')
                    sdfiles.append(sdfilename)
                    newbunch = False
                newf.write(line)
    else:
        print(" *** Error: sdfile could not be processed! *** ")
        exit()
    print("Done!")
    return count, sdfiles

def smiles2png(smiles, labels=None, ref=None):
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
    align = False
    if ref:
        if str(ref)[-3::] == 'sdf':
            try:
                suppl = Chem.SDMolSupplier(ref)
                template = Chem.MolFromSmiles(Chem.MolToSmiles(suppl[0]))
                AllChem.Compute2DCoords(template)
                align=True
            except:
                print("WARNING: ref molecule could not be used!")
                Exception
        else:
            try:
                template = Chem.MolFromSmiles(ref)
                AllChem.Compute2DCoords(template)
                align=True
            except:
                print("WARNING: ref molecule could not be used!")
                Exception

    mols_to_draw = []
    for i,smile in enumerate(mols):
        mol = Chem.MolFromSmiles(smile)
        if labels:
            molname = label[i]
            mol.SetProp('_Name',molname)
        else:
            molname = "mol_"+str(i+1)
            mlabel.append(molname)
            mol.SetProp('_Name',molname)
        print(i+1,molname)
        if align:
            AllChem.Compute2DCoords(mol)
            mcs = rdFMCS.FindMCS([template, mol])
            patt = Chem.MolFromSmarts(mcs.smartsString)
            query_match = mol.GetSubstructMatch(patt)
            template_match = template.GetSubstructMatch(patt)
            try:
                rms = AllChem.AlignMol(mol, template, atomMap=list(zip(query_match,template_match)))
            except RuntimeError:
                continue
        Draw.MolToFile(mol, molname+'.png', legend=molname, size=(600,600))
        mols_to_draw.append(mol)
    try:
        img = Draw.MolsToGridImage(mols_to_draw, molsPerRow=4, legends=[x.GetProp('_Name') for x in mols_to_draw], subImgSize=(350,350))
        img.save("all_mols_grid_image.png")
    except:
        print("WARNING: all_mols_grid_image.png could not be created!")
        Exception

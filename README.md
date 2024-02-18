# LigPrepper

Simple program to quickly prepare 3D SDF/MOL2/PDBQT files from SMILES for ligands using rdkit & openbabel.

Also, split sdf and convert formats sdf<->pdbqt for Autodock Vina

[![PyPI version](https://badge.fury.io/py/LigPrepper.svg)](https://badge.fury.io/py/LigPrepper) [![Downloads](https://pepy.tech/badge/ligprepper)](https://pepy.tech/project/ligprepper)

## Installation:

```
pip install LigPrepper
```

## Usage:

#### For single SMILES

```
LigPrepper.smiles2sdf('c1ccncc1', labels='pyridine')

LigPrepper.smiles2mol2('c1ccncc1', labels='pyridine')

LigPrepper.smiles2pdbqt('c1ccncc1', labels='pyridine')
```

#### For a list of SMILES

```
LigPrepper.smiles2sdf(smiles_list, labels=labels_list, mergesdf=False)
```

#### Want to align molecule to reference molecule in SDF format

```
LigPrepper.smiles2sdf(smiles_list, labels=labels_list, ref='ref.sdf', mergesdf=False)
```

#### convert formats

```
LigPrepper.sdf2pdbqt("your-molecule.sdf")
LigPrepper.pdbqt2sdf("your-molecule.pdbqt")
```

#### split sdf

```
Ligprepper.splitsdf("your-molecules.sdf",outputdir=None,parts=0,molspersdf=1,firstpart=1)

Arguments:
    outputdir   - directory to save output sdf files
    parts       - number of parts to split
    molspersdf  - number of molecules to be chunked into (parts is neglected if this is used)
    firstpart   - naumber of the first part 
```

#### Draw 2D structures

```
LigPrepper.smiles2png(smiles_list, labels=labels_list)

LigPrepper.smiles2png(smiles_list, labels=labels_list, ref='ref smiles')

LigPrepper.smiles2png(smiles_list, labels=labels_list, ref='ref.sdf')
```
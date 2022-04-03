# LigPrepper

Simple program to quickly prepare 3D SDF files from SMILES for ligands using rdkit.

## Installation:

```
pip install LigPrepper
```

## Usage:

#### For single SMILES

```
LigPrepper.smiles2sdf('c1ccncc1', labels='pyridine')
```

#### For a list of SMILES

```
LigPrepper.smiles2sdf(smiles_list, labels=labels_list, mergesdf=False)
```

#### Have reference molecule to align to in sdf format

```
LigPrepper.smiles2sdf(smiles_list, labels=labels_list, ref='examples/ref.sdf', mergesdf=False)
```


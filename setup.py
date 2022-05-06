from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="LigPrepper",
    version="0.3.0",
    author="Abhishek Kognole",
    author_email="aakognole@gmail.com",
    description="Simple program to quickly prepare ligand 3d SDF/MOL2/PDBQT from smiles",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aakognole/ligprepper",
    project_urls={
        "Bug Tracker": "https://github.com/aakognole/ligprepper/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    install_requires=['rdkit-pypi'],
    python_requires=">=3.6",
    zip_safe=False
)

from setuptools import setup, find_packages 

setup(name='LigPrepper', 
      version='0.1.0', 
      description='program to quickly prepare ligand 3d SDF from smiles', 
      long_description="",
      url='https://github.com/aakognole/ligprepper',
      author='Abhishek Kognole',
      author_email='aakognole@gmail.com',
      license='MIT',
      packages=find_packages(), 
      include_package_data=False,
      install_requires=['rdkit-pypi'], 
      zip_safe=False)

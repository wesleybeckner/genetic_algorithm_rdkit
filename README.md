### Genetic Algorithm Can Accurately Navigate Structure Space for Ionic Liquids with Desirable Properties

Genetic algorithms that use QSPR models in their fitness test can be used to navigate molecular structure space and design new materials.

1. X versions - no RDKit integration
2. X versions - RDKit integration
    1. atomic mutations, alpha-num fitness test, fixed bonds, non-random starting structure
    2. atomic mutations, similarity map fitness test, fixed bonds, non-random starting structure
    3. atomic mutations, similarity map fitness test, bond mutations, non-random starting structure
    4. atomic mutations, similarity map fitness test, bond mutations, random starting structure
3. X versions - viscosity fitness tests
4. X versions - performance upgrades using clustering/similarity maps
5. X versions - different search structures

this package requires rdkit

# MAC OS install

[link](http://www.rdkit.org/docs/Install.html)
after creating a conda environment:
conda install postgresql
conda install rdkit

# Windows install

conda install -c rdkit rdkit

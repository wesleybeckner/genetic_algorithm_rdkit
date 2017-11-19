### Genetic Algorithm Can Accurately Navigate Structure Space for Ionic Liquids with Desirable Properties

Genetic algorithms that use QSPR models in their fitness test can be used to navigate molecular structure space and design new materials.

1. X versions - no RDKit integration
2. X versions - RDKit integration
    1. atomic mutations, alpha-num fitness test, fixed bonds, non-random starting structure
    2. atomic mutations, similarity map fitness test, fixed bonds, non-random starting structure
    3. atomic mutations, similarity map fitness test, bond mutations, non-random starting structure
    4. atomic mutations, similarity map fitness test, bond mutations, random starting structure
3. X versions - density/viscosity fitness tests
    1. get other 20k density data online
    2. follow viscosity protocol for density property
4. X versions - performance upgrades using clustering/similarity maps
    1. compare RDKit similarity map types
    2. cationic vs anionic biasing
    3. search structures (children/parent pool, chromosome crossover)
5. X versions - expanded to IL mixtures
6. X versions - GROMACS integration
    1. generate pdb systems
    2. suggest GROMACS commands

this package requires rdkit

# Linux OS install

conda create -n py36 python=3.6 anaconda

conda install -c rdkit rdkit

# MAC OS install

[link](http://www.rdkit.org/docs/Install.html)
after creating a conda environment:
conda install postgresql
conda install rdkit rdkit

# Windows install

conda install -c rdkit rdkit



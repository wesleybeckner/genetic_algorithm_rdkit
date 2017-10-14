import rdkit
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import ShowMol
import statistics
import time
import random
import sys
"""
This GA uses RDKit to make atomic mutations to a starting imidazole.
The starting structure is not random.
Fitness test uses RDKit FingerprintSimilarity.
Number of atoms in parent/children are fixed.
"""
class Benchmark:
    @staticmethod
    def run(function):
        timings = []
        stdout = sys.stdout
        for i in range(5):
            sys.stdout = None
            startTime = time.time()
            function()
            seconds = time.time() - startTime
            sys.stdout = stdout
            timings.append(seconds)
            mean = statistics.mean(timings)
            print("{} {:3.2f} {:3.2f}".format(
                1 + i, mean,
                statistics.stdev(timings, mean) if i > 1 else 0))
class Chromosome(Chem.rdchem.Mol):
    def __init__(self, genes, fitness):
        Chem.rdchem.Mol.__init__(self)
        self.Genes = genes
        self.Fitness = fitness
        self.Mol = Chem.MolFromSmiles(genes)
        self.RWMol = Chem.MolFromSmiles(genes)
        self.RWMol = Chem.RWMol(Chem.MolFromSmiles(genes))

def _generate_parent(length, geneSet, get_fitness):
    genes = "[NH+]1C=CN(C1C)C"
    fitness = get_fitness(genes)
    return Chromosome(genes, fitness)

def _mutate(parent, geneSet, get_fitness, fcat):
    def replace_atom(childGenes, geneSet, oldGene):
        if childGenes.RWMol.GetAtomWithIdx(oldGene).IsInRing() == True:
            genes = Chem.MolToSmiles(parent.Mol)
            return Chromosome(genes, 0)
        newGene = random.sample(geneSet, 1)[0]
        childGenes.RWMol.GetAtomWithIdx(oldGene).SetAtomicNum(newGene) 
        return childGenes  
    def add_atom(childGenes, geneSet, oldGene):
        newGeneNumber = childGenes.RWMol.GetNumAtoms()  
        newGene = random.sample(geneSet, 1)[0]
        childGenes.RWMol.AddAtom(Chem.Atom(newGene))
        childGenes.RWMol.AddBond(newGeneNumber,oldGene,Chem.BondType.SINGLE) 
        return childGenes
    def remove_atom(childGenes, geneSet, oldGene):
        if childGenes.RWMol.GetAtomWithIdx(oldGene).GetExplicitValence() != 1:
            genes = Chem.MolToSmiles(parent.Mol)
            return Chromosome(genes, 0)
        childGenes.RWMol.RemoveAtom(oldGene)
        return childGenes
    def add_rdkit_fragment(childGenes, geneSet, oldGene):
        try:
            newGene = Chromosome(Chem.MolToSmiles(fparams.GetFuncGroup(\
            random.sample(range(fparams.GetNumFuncGroups()), 1)[0])),0)
        except:
            return 0
        oldGene = oldGene + newGene.Mol.GetNumAtoms()
        combined = Chem.EditableMol(Chem.CombineMols(newGene.Mol,childGenes.Mol))
        combined.AddBond(1,oldGene,order=Chem.rdchem.BondType.SINGLE)
        combined.RemoveAtom(0)
        try:
            childGenes = Chromosome(Chem.MolToSmiles(childGenes),0)  
            return childGenes
        except:
            return 0
    def add_custom_fragment(childGenes, geneSet, oldGene):
        newGene = Chromosome(fcat.GetEntryDescription(\
        random.sample(range(fcat.GetNumEntries()), 1)[0]),0)
        oldGene = oldGene + newGene.Mol.GetNumAtoms()
        combined = Chem.EditableMol(Chem.CombineMols(newGene.Mol,childGenes.Mol))
        combined.AddBond(0,oldGene,order=Chem.rdchem.BondType.SINGLE)
        childGenes = combined.GetMol()   
        try:
            childGenes = Chromosome(Chem.MolToSmiles(childGenes),0)  
            return childGenes
        except:
            return 0
    childGenes = Chromosome(parent.Genes,0)
    oldGene = random.sample(range(childGenes.RWMol.GetNumAtoms()), 1)[0]
    mutate_operations = [add_rdkit_fragment, add_custom_fragment, remove_atom, \
            replace_atom, add_atom]
    i = random.choice(range(len(mutate_operations)))
    childGenes = mutate_operations[i](childGenes, geneSet, oldGene)
    try:
        childGenes.RWMol.UpdatePropertyCache(strict=True)
        Chem.SanitizeMol(childGenes.RWMol)
        genes = Chem.MolToSmiles(childGenes.RWMol)
        fitness = get_fitness(genes, target)
        return Chromosome(genes, fitness)
    except:
        return Chromosome(parent.Genes, 0)
        

def get_best(get_fitness, targetLen, optimalFitness, geneSet, display,\
        show_ion, fcat):
    random.seed()
    bestParent = _generate_parent(targetLen, geneSet, get_fitness)
    display(bestParent)
    if bestParent.Fitness >= optimalFitness:
        return bestParent
    while True:
        child = _mutate(bestParent, geneSet, get_fitness, fcat)
        if bestParent.Fitness >= child.Fitness:
            continue
        display(child)
        if child.Fitness >= optimalFitness:
            show_ion()
            return child
        bestParent = child

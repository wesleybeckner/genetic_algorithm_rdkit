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

def _generate_parent(length, geneSet, get_fitness):
    genes = "NCNC[NH+]1C=CN(C1C)C"
    fitness = get_fitness(genes)
    return Chromosome(genes, fitness)

def _mutate(parent, geneSet, get_fitness):
    childGenes = parent
    while True:
        oldGene, alternate = random.sample(range(0,childGenes.Mol.GetNumAtoms()), 2)
        if childGenes.Mol.GetAtomWithIdx(alternate).IsInRing() == False and\
                childGenes.Mol.GetAtomWithIdx(oldGene).IsInRing() == False:
            break
    newGene = random.sample(geneSet, 1)[0]
    try:
        childGenes.Mol.GetAtomWithIdx(alternate).SetAtomicNum(newGene) if\
              childGenes.Mol.GetAtomWithIdx(oldGene).GetFormalCharge() != 0 else\
              childGenes.Mol.GetAtomWithIdx(oldGene).SetAtomicNum(newGene)
        Chem.SanitizeMol(childGenes.Mol)
    except:
        print("Sanitization Failure")
    genes = Chem.MolToSmiles(childGenes.Mol)
    fitness = get_fitness(genes)
    return Chromosome(genes, fitness)

def get_best(get_fitness, targetLen, optimalFitness, geneSet, display,\
        show_ion):
    random.seed()
    bestParent = _generate_parent(targetLen, geneSet, get_fitness)
    display(bestParent)
    if bestParent.Fitness >= optimalFitness:
        return bestParent
    while True:
        child = _mutate(bestParent, geneSet, get_fitness)
        if bestParent.Fitness >= child.Fitness:
            continue
        display(child)
        if child.Fitness >= optimalFitness:
            show_ion()
            return child
        bestParent = child

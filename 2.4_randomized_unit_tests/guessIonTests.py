#! /usr/bin/python
import genetic #our genetic engine code
import rdkit
import os
from rdkit import RDConfig
from rdkit.Chem import FragmentCatalog
from rdkit import RDConfig
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Draw import ShowMol
import random
import unittest
import datetime

class GuessIonTests(unittest.TestCase):
    geneSet = genetic.generate_geneset()
    
    def test_1_butyl_2_3_dimethyl_1H_imidazolium(self):
        target = "CN1C=C[N+](=C1)C"
        self.guess_password(target)

    def test_benchmark(self):
        genetic.Benchmark.run(self.test_1_butyl_2_3_dimethyl_1H_imidazolium)

    def guess_password(self, target):
        startTime = datetime.datetime.now()

        def fnGetFitness(genes):
            return get_fitness(genes, target)

        def fnDisplay(candidate):
            display(candidate, startTime)

        def fnShowIon():
            show_ion(target)

        optimalFitness = get_fitness(target, target)
        best = genetic.get_best(fnGetFitness,\
			optimalFitness, self.geneSet, fnDisplay,\
                        fnShowIon, target)
    
def display(candidate, startTime):
    timeDiff = datetime.datetime.now() - startTime
    print("{}\t{}\t{}".format(
	candidate.Genes, candidate.Fitness, timeDiff))
    
def get_fitness(genes, target):
    ms = [Chem.MolFromSmiles(target), Chem.MolFromSmiles(genes)]
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    return DataStructs.FingerprintSimilarity(fps[0],fps[1])

def show_ion(target):
    mol = Chem.MolFromSmiles(target)
    print("{}\t{}".format("number of atoms: ", mol.GetNumAtoms()))
    for atom in mol.GetAtoms():
        print("{}\t{}\t{}".format("atom %s ring status and valence:"\
                %atom.GetSymbol(), atom.IsInRing(), atom.GetExplicitValence()))



if __name__ == '__main__':
    unittest.main()



#from rdkit import Chem
#import argparse
#
#p = argparse.ArgumentParser(description='process genes with RDKit.')
#p.add_argument('genes', metavar='Genes', nargs="*", type=str,\
#        help='the SMILES alphanumeric string of a molecular object')
#args = p.parse_args()
#temp = Chem.MolFromSmiles(args.genes[0]).UpdatePropertyCache(strict=True)
#Chem.SanitizeMol(temp)

from openfermion.chem import MolecularData
import sys
from vqemol import VQEMOL

#path_list = sys.argv
#path = path_list[1]

basis = "sto-3g" #basis set
multiplicity = 1 #spin multiplicity
charge = 0   #total charge for the molecule
geometry = [("H",(0.37,0,0)),("H",(-0.37,0,0))]
molecule = MolecularData(geometry, basis, multiplicity, charge)
mol = VQEMOL(molecule)

#ans = mol.runVQE('UCCSD',init_guess='MP2')
mol.ansatz ='UCCSD'
mol.t_num = 1
mol.opt_param = [0.22555712, 0.         , 0.        ]
print(mol.cal_NaturalOrb())
print(mol.count_elenum())


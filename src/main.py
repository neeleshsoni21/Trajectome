"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

import os
import logging
import time

from System_Class import System
from Simulation_Class import Simulation

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":
	system = System(L=1000)

	# Example: Adding proteins
	#system.add_protein("Protein1", center=[1,1,1], radius=30, mass=1000, diffcoff=0.01, color=0)
	#system.add_protein("Protein2", center=[2,2,2], radius=60, mass=2000, diffcoff=0.02, color=1)
	

	data_path = './data/'
	pdbfile1 = os.path.join(data_path,'pdb/7r5j_C8_CA_center_v1.cif')
	fastafile1 = os.path.join(data_path,'fasta/NPC_Whole.fasta')

	#pdbfile1 = os.path.join(data_path,'pdb/1v5w_C8.cif')
	#fastafile1 = os.path.join(data_path,'fasta/1v5w_A.fasta')
	
	pdbfile2 = os.path.join(data_path, 'pdb/Nup2.pdb')
	fastafile2 = os.path.join(data_path, 'fasta/Nup2.fasta')

	pdbfile3 = os.path.join(data_path, 'pdb/Nup85.pdb')
	fastafile3 = os.path.join(data_path, 'fasta/Nup85.fasta')

	#TODO: Add rigid body and resolution option. Add logging lines also
	#system.add_protein_from_structure("NPC", pdbfile1, fastafile1, diffcoff=0.0, color=0, centerize = False)

	for idx in range(0,3):
		system.add_protein_from_structure("NUP2."+str(idx), pdbfile2, fastafile2, diffcoff=0.001, color=1, centerize = False)
		system.add_protein_from_structure("NUP85."+str(idx), pdbfile3, fastafile3, diffcoff=0.001, color=2, centerize = False)
		
	#Get all ProteinStructure objects and obtain its hierarchy to append in the root hierarchy
	for prot in system.proteins:
		system.h_root.add_child(prot.hier)

	# Example: Adding interaction
	#protein1 = system.proteins[0]
	#protein2 = system.proteins[1]
	#system.add_interaction(protein1, protein2, "binding", 10.0)

	# Adding Excluded Volume Restraint and boundary conditions
	system.add_excluded_volume_restraint()
	system.apply_boundary_conditions()

	# time this step
	t1 = time.time()

	simulation = Simulation(system, output_dir='./output/', time_steps=500, temperature=300)
	simulation.run()

	dt2 = time.time() - t1
	print("Time taken for Sim:",dt2)

	print("Score: kcal/mol")

"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

import os
import logging

from System_Class import System
from Simulation_Class import Simulation

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":
	system = System()

	# Example: Adding proteins
	#system.add_protein("Protein1", center=[1,1,1], radius=30, mass=1000, diffcoff=0.01, color=0)
	#system.add_protein("Protein2", center=[2,2,2], radius=60, mass=2000, diffcoff=0.02, color=1)
	

	data_path = './data/'
	pdbfile1 = os.path.join(data_path,'pdb/7r5j_C8.cif')
	fastafile1 = os.path.join(data_path,'fasta/NPC_Whole.fasta')
	
	pdbfile2 = os.path.join(data_path, 'pdb/Nup2.pdb')
	fastafile2 = os.path.join(data_path, 'fasta/Nup2.fasta')

	#TODO: Add rigid body and resolution option. Add logging lines also
	system.add_protein_from_structure("NUP2", pdbfile1, fastafile1, diffcoff=0.01, color=0)
	system.add_protein_from_structure("NPC", pdbfile2, fastafile2, diffcoff=0.01, color=0)


	for prot in system.proteins:
		print(prot)
		system.h_root.add_child(prot.hier)

	# Example: Adding interaction
	#protein1 = system.proteins[0]
	#protein2 = system.proteins[1]
	#system.add_interaction(protein1, protein2, "binding", 10.0)

	# Adding Excluded Volume Restraint and boundary conditions
	system.add_excluded_volume_restraint()
	system.apply_boundary_conditions()

	simulation = Simulation(system, output_dir='./output/', time_steps=1000, temperature=300)
	simulation.run()

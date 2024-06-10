"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

import os
import logging
import time

from System_Class import System
from Simulation_Class import Simulation
from Interaction_Class import Interaction

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if __name__ == "__main__":
	system = System(L=2000)

	# Example: Adding proteins
	#system.add_protein("Protein1", center=[1,1,1], radius=30, mass=1000, diffcoff=0.01, color=0)
	#system.add_protein("Protein2", center=[2,2,2], radius=60, mass=2000, diffcoff=0.02, color=1)
	
	data_path = './data/'

	'''
	NPC_Subunits=[]
	#HUMAN: IMPORT each NPC subunit from different pdb
	for idx in range(0,8):
		pdbfile1 = os.path.join(data_path,
			'pdb/Human_Proteins/C8_Symmetry/7r5j_C8_SU_'+str(idx+1)+'.cif')

		npc_subunit_i = system.add_protein_from_structure("NPC-SU-"+str(idx+1), 
			pdbfile1, pdb_multimodel=False, resolution = 100, 
			diffcoff=0.0, color=1, centerize = False)

		NPC_Subunits.append(npc_subunit_i)
	'''
	
	NPC_Subunits=[]
	#YEAST:IMPORT entire NPC from single pdb
	for idx in range(0,8):
		pdbfile1 = os.path.join(data_path,
			'pdb/Yeast_Proteins/Subunits/FullNPC_SU_'+str(idx+1)+'.cif')

		npc_subunit_i = system.add_protein_from_structure("NPC-SU-"+str(idx+1), 
			pdbfile1, pdb_multimodel=False, resolution = 100, 
			diffcoff=0.0, color=1, centerize = False)

		NPC_Subunits.append(npc_subunit_i)
	
	
	NGH_Prots=[];
	NGH_Prots1=[];NGH_Prots2=[];NGH_Prots3=[]

	NGH_Proteins_Copies = 100
	#YEAST NPC interactors
	for idx in range(0,NGH_Proteins_Copies):

		pdbfile2 = os.path.join(data_path, 'pdb/Yeast_Proteins/Nup2.pdb')
		ngh_prot1 = system.add_protein_from_structure("NUP2."+str(idx), 
			pdbfile2, pdb_multimodel=False, resolution = 50, 
			diffcoff=0.01, color=2, centerize = False)
		NGH_Prots.append(ngh_prot1)
		NGH_Prots1.append(ngh_prot1)

		pdbfile3 = os.path.join(data_path, 'pdb/Yeast_Proteins/Nup85.pdb')
		ngh_prot2 = system.add_protein_from_structure("NUP85."+str(idx), 
			pdbfile3, pdb_multimodel=False, resolution = 50,
			diffcoff=0.01, color=3, centerize = False)
		NGH_Prots.append(ngh_prot2)
		NGH_Prots2.append(ngh_prot2)

		pdbfile4 = os.path.join(data_path, 'pdb/Yeast_Proteins/Nup85.pdb')
		ngh_prot3 = system.add_protein_from_structure("NUP85."+str(idx), 
			pdbfile4, pdb_multimodel=False, resolution = 50,
			diffcoff=0.01, color=4, centerize = False)
		NGH_Prots.append(ngh_prot3)
		NGH_Prots3.append(ngh_prot3)

	#Get all ProteinStructure objects and obtain its hierarchy to append in the root hierarchy
	for prot in system.proteins:
		system.h_root.add_child(prot.hier)
	
	
	#-------------------------
	#ADD DISTANCE RESTRAINTS
	#-------------------------
	#YEAST NPC interactors
	for idx in range(0,NGH_Proteins_Copies):
		prot1_tuple = ("NUP85."+str(idx), 0); prot2_tuple = ("NUP2."+str(idx), 0)
		Int_obj = Interaction(system)
		Int_obj.add_binding_restraint(prot1_tuple, prot2_tuple, "binding", 20, 0.00001)
		pass

	
	# Adding Excluded Volume Restraint and boundary conditions
	system.add_excluded_volume_restraint()
	system.apply_boundary_conditions()
	#system.Iterate_Hierarchy()
	
	###system.add_membrane_restraint(NGH_Prots)


	system.apply_nucleoplasm_boundary_conditions(NGH_Prots1)
	system.apply_nucleoplasm_boundary_conditions(NGH_Prots2)
	system.apply_cytoplasm_boundary_conditions(NGH_Prots3)
	
	bounding_box="nucleoplasm"
	system.shuffle_neighborhood_proteins(bounding_box, NGH_Prots1)
	bounding_box="nucleoplasm"
	system.shuffle_neighborhood_proteins(bounding_box, NGH_Prots2)
	bounding_box="cytoplasm"
	system.shuffle_neighborhood_proteins(bounding_box, NGH_Prots3)

	# time this step
	t1 = time.time()

	#simulation_time = 0.0001 #In Seconds
	simulation = Simulation(system, output_dir='./output/', simulation_time = 0.0001, temperature=300)
	simulation.run()

	dt2 = time.time() - t1
	print("Time taken for Sim:",dt2)

	print("Score: kcal/mol")

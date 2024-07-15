"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

import IMP
import IMP.atom
import IMP.pmi
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.topology
import IMP.pmi.tools

import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

from Protein_Class import Protein, ProteinStructure
from Interaction_Class import Interaction

class System:
	def __init__(self, L=500, R=1300, K_BB=10.0):
		self.model = IMP.Model()

		self.sys = IMP.pmi.topology.System(self.model)
		self.state = self.sys.create_state()
		self.h_root = self.state.build()

		#self.EV_RESOLUTION = 50
		self.proteins = []
		self.interactions = []
		self.restraints = []
		self.L = L
		#self.R = R
		self.K_BB = K_BB

		tor_th      = 45.0
		self.tor_R = 390.0 + 150.0
		self.tor_r = 150.0 - tor_th/2.0
		self.mem_thickness = 30

		self.create_bounding_box_and_pbc()
		#logger.info(f"Initialized System_Class with file: {file_path}, format: {file_format}")

	def create_bounding_box_and_pbc(self):
		# PBC cytoplasm bounding sphere:
		#self.pbc_sphere = IMP.algebra.Sphere3D([0, 0, 0], self.R)

		# Outer bounding box for simulation:
		self.bb = IMP.algebra.BoundingBox3D(
			IMP.algebra.Vector3D(-self.L/2, -self.L/2, -self.L/2),
			IMP.algebra.Vector3D(self.L/2, self.L/2, self.L/2)
		)

		membrane_start_z = self.mem_thickness + self.tor_r/2.0
		self.bb_nucleus = IMP.algebra.BoundingBox3D(
			IMP.algebra.Vector3D(-self.L/2, -self.L/2, -self.L/2),
			IMP.algebra.Vector3D(self.L/2, self.L/2, -1*membrane_start_z)
		)

		membrane_start_z = self.mem_thickness + self.tor_r/2.0
		self.bb_cytoplasm = IMP.algebra.BoundingBox3D(
			IMP.algebra.Vector3D(-self.L/2, -self.L/2, 1*membrane_start_z),
			IMP.algebra.Vector3D(self.L/2, self.L/2, self.L/2)
			
		)

		self.bb_harmonic = IMP.core.HarmonicUpperBound(0, self.K_BB)
		self.outer_bbss = IMP.core.BoundingBox3DSingletonScore(self.bb_harmonic, self.bb)

		self.bb_nuc_harmonic = IMP.core.HarmonicUpperBound(0, self.K_BB)
		self.bb_nuc_bbss = IMP.core.BoundingBox3DSingletonScore(self.bb_nuc_harmonic, self.bb_nucleus)

		self.bb_cyt_harmonic = IMP.core.HarmonicUpperBound(0, self.K_BB)
		self.bb_cyt_bbss = IMP.core.BoundingBox3DSingletonScore(self.bb_cyt_harmonic, self.bb_cytoplasm)

	def add_protein(self, name, center, radius, mass, diffcoff, color):
		protein = Protein(self.model, name, center, radius, mass, diffcoff, color)
		self.proteins.append(protein)

	def add_protein_from_structure(self, name, pdbfile, pdb_multimodel,  resolution, diffcoff, color, centerize):
		protein = ProteinStructure(self.model, self.state, self.h_root, name, pdbfile, pdb_multimodel, resolution, diffcoff, color, centerize)
		self.proteins.append(protein)
		return protein

	def add_interaction(self, protein1, protein2, interaction_type, strength):
		interaction = Interaction(protein1, protein2, interaction_type, strength)
		self.interactions.append(interaction)

	def add_restraint(self, restraint):
		self.restraints.append(restraint)

	def apply_boundary_conditions(self):
		rgbmembers = [protein.prb.get_particle().get_index() for protein in self.proteins]
		self.add_restraint(IMP.container.SingletonsRestraint(self.outer_bbss, IMP.container.ListSingletonContainer(self.model, rgbmembers)))

	def apply_nucleoplasm_boundary_conditions(self, proteins):
		prot_members = [protein.prb.get_particle().get_index() for protein in proteins]
		self.add_restraint(IMP.container.SingletonsRestraint(self.bb_nuc_bbss, IMP.container.ListSingletonContainer(self.model, prot_members)))

	def apply_cytoplasm_boundary_conditions(self, proteins):
		prot_members = [protein.prb.get_particle().get_index() for protein in proteins]
		self.add_restraint(IMP.container.SingletonsRestraint(self.bb_cyt_bbss, IMP.container.ListSingletonContainer(self.model, prot_members)))


	def add_excluded_volume_restraint(self):
		rgbmembers=[]
		#ADD ExcludedVolume by considering each rigidbody/protein as a sphere. Too coarse
		###rgbmembers = [protein.prb.get_particle().get_index() for protein in self.proteins]

		#ADD ExcludedVolume by considering each rigidbody/protein MEMBERS (residues) as a sphere. Too fine grained
		for protein in self.proteins:
			#print("p:",protein.prb.get_name(), protein.prb, protein.prb.get_particle())
			for residue_idxs in protein.prb.get_member_particle_indexes():
				rgbmembers.append(residue_idxs)

		#sys.exit()
		ev = IMP.core.ExcludedVolumeRestraint(IMP.container.ListSingletonContainer(self.model, rgbmembers), 100, 10)
		self.add_restraint(ev)

	def add_membrane_restraint2(self, NGH_proteins):

		#PROT1_P1 =IMP.atom.Selection(self.system.h_root,molecule=prot1).get_selected_particles()[p1_idx]
		rgbmembers=[]
		for protein in self.proteins:
			for residue_idxs in protein.prb.get_member_particle_indexes():
				rgbmembers.append(residue_idxs)

		import IMP.npc
		import IMP.pmi.restraints.npc
		#ev = IMP.npc.MembraneExclusionRestraint(IMP.container.ListSingletonContainer(self.model, rgbmembers), 100, 10)
		#ev = IMP.npc.MembraneExclusionRestraint(self.model, IMP.container.ListSingletonContainer(self.model, rgbmembers), 100, 10)
		#self.add_restraint(ev)

		r = IMP.pmi.restraints.npc.MembraneExclusionRestraint(
            #hier=root_hier, protein=(1, 1, "Nup84"), label='Test',
            hier=self.h_root, protein=(1, 1, "Nup2.0"), label='Test',
            tor_R=30., tor_r=10.)
		r.add_to_model()
		self.add_restraint(r)

		return

	def add_membrane_restraint(self, NGH_proteins):

		#PROT1_P1 =IMP.atom.Selection(self.system.h_root,molecule=prot1).get_selected_particles()[p1_idx]

		#prot = NGH_proteins[0]
		#print("Prot name:",prot.hier.get_name())

		import IMP.npc
		import IMP.pmi.restraints.npc
		from MembraneRestraint_Class import MembraneRestraint

		for prot in NGH_proteins:
			print(prot.mol.get_name())
			print(prot.hier.get_name())
			#mr = IMP.pmi.restraints.basic.MembraneRestraint2(self.h_root,
			mr_obj = MembraneRestraint(self.h_root,
				objects_inside=None,
				objects_above=None,
				#objects_below=[prot.hier.get_name()],
				objects_below=[prot.mol.get_name()],
				#objects_below=["NUP2.9_Chain-A_frag_83-132"],
				
				#objects_below=[(41, 80, "helix_2")],
				weight=10)

			#self.add_restraint(mr_obj.rs)
			#self.add_restraint(mr_obj.mr)

			mr_obj2 = IMP.pmi.restraints.npc.MembraneExclusionRestraint(
				hier=self.h_root, protein=prot.mol.get_name(), label='Test',
				tor_R=540., tor_r=127)

			self.add_restraint(mr_obj2.rs)
			#self.add_restraint(mr_obj2.mex)
			#mr.add_to_model()

		return

	def add_membrane_exclusion_restraint(self, NGH_proteins):

		#PROT1_P1 =IMP.atom.Selection(self.system.h_root,molecule=prot1).get_selected_particles()[p1_idx]

		#prot = NGH_proteins[0]
		#print("Prot name:",prot.hier.get_name())

		import IMP.npc
		import IMP.pmi.restraints.npc
		from MembraneExclusionRestraint import 	MembraneExclusionRestraint

		for prot in NGH_proteins:
			print(prot.mol.get_name())
			print(prot.hier.get_name())
			#mr = IMP.pmi.restraints.basic.MembraneRestraint2(self.h_root,
			#mr_obj = MembraneRestraint(self.h_root,
			#	objects_inside=None,
			#	objects_above=None,
			#	#objects_below=[prot.hier.get_name()],
			#	objects_below=[prot.mol.get_name()],
			#	#objects_below=["NUP2.9_Chain-A_frag_83-132"],
				
			#	#objects_below=[(41, 80, "helix_2")],
			#	weight=10)

			#self.add_restraint(mr_obj.rs)
			#self.add_restraint(mr_obj.mr)

			#mr_obj2 = IMP.pmi.restraints.npc.MembraneExclusionRestraint(
			mr_obj2 = MembraneExclusionRestraint(
				hier=self.h_root, protein=prot.mol.get_name(), label='Test',
				tor_R=540., tor_r=127)

			self.add_restraint(mr_obj2.rs)
			#self.add_restraint(mr_obj2.mex)
			#mr.add_to_model()

		return

	def update(self):
		for interaction in self.interactions:
			interaction.apply_interaction()

	def Iterate_Hierarchy(self):

		#This extra for dont required with prism rmfs
		for state in self.h_root.get_children(): 
			print("state:",state)
			for prot in state.get_children():
				#sel = IMP.atom.Selection(prot,resolution=1)
				print("prot:",prot)
				for c in prot.get_children():
					print("c:",c)
					name = c.get_name()
					leaf = IMP.atom.get_leaves(c)
					
					protname = prot.get_name()
					coordsR = IMP.core.XYZR(leaf[0])
					xyz = coordsR.get_coordinates()
					radius = coordsR.get_radius()
					print(name, leaf, protname, xyz, radius)
		return

	def shuffle_neighborhood_proteins(self, bounding_box,  NGH_proteins):
		bb_restraint = None
		if bounding_box=='nucleoplasm':
			bb_restraint = self.bb_nucleus
		elif bounding_box=='cytoplasm':
			bb_restraint = self.bb_cytoplasm
		else:
			bb_restraint = self.bb

		for protein in NGH_proteins:

			#To skip fixed particles
			if protein.diffcoff==0:
				continue

			print("BB:",protein.name)
			translation = IMP.algebra.get_random_vector_in(bb_restraint)
			rotation = IMP.algebra.get_random_rotation_3d()
			transformation = IMP.algebra.Transformation3D(rotation, translation)
			print(translation)

			IMP.core.transform(protein.prb, transformation)

		pass

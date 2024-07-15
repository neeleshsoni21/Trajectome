"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024

Attributes:
	logger (TYPE): Description
"""

import logging

import IMP
import IMP.atom
import IMP.core


logger = logging.getLogger(__name__)

class Interaction:
	"""
	A class to represent interactions within a protein simulation system.

	Attributes:
		system (System): The system to which the interactions belong.
	"""
	
	def __init__(self, system):
		"""
		Initializes the Interaction instance with the provided system.

		Args:
			system (System): The system to which the interactions belong.

		Returns:
			None
		"""
		self.system = system
		logger.info("Interaction instance created.")

	def add_binding_restraint(self, prot1_tuple, prot2_tuple, name, mean_dist, kappa):
		"""
		Adds a binding restraint between two proteins.

		Args:
			prot1_tuple (tuple): A tuple containing the first protein and its index.
			prot2_tuple (tuple): A tuple containing the second protein and its index.
			name (str): The name of the restraint.
			mean_dist (float): The mean distance for the restraint.
			kappa (float): The kappa value for the harmonic function.

		Returns:
			None
		"""
		prot1, p1_idx = prot1_tuple
		prot2, p2_idx = prot2_tuple

		PROT1_P1 =IMP.atom.Selection(self.system.h_root,molecule=prot1).get_selected_particles()[p1_idx]

		PROT2_P1 =IMP.atom.Selection(self.system.h_root,molecule=prot2).get_selected_particles()[p2_idx]

		dr = IMP.core.DistanceRestraint(self.system.model, IMP.core.Harmonic(mean_dist, kappa), PROT1_P1, PROT2_P1)

		dr.set_name(name)
		
		self.system.add_restraint(dr)

		logger.info(f"Added binding restraint between {prot1} and {prot2} with mean distance {mean_dist} and kappa {kappa}.")


	def add_distance_restraint(self, prot1, prot2, dist, k):
		"""
		Adds a distance restraint between two proteins.

		Args:
			prot1 (Protein): The first protein.
			prot2 (Protein): The second protein.
			dist (float): The target distance for the restraint.
			k (float): The force constant for the restraint.

		Returns:
			None
		"""
		cr = IMP.atom.create_distance_restraint(
			IMP.atom.Selection(prot1.mol), 
			IMP.atom.Selection(prot2.mol), dist, k)

		self.system.add_restraint(cr)

		logger.info(f"Added distance restraint between {prot1.name} and {prot2.name} with distance {dist} and kappa {k}.")




	'''
	
	def calculate_super_particle_forces(self, super_particles):
		forces = []
		for i, sp1 in enumerate(super_particles):
			force = np.zeros(3)
			for j, sp2 in enumerate(super_particles):
				if i != j:
					r = sp1.position - sp2.position
					distance = np.linalg.norm(r)
					if distance < self.cutoff:
						f = self.lj_force(distance) * (r / distance)
						force += f
			forces.append(force)
		return forces

	def lj_force(self, r):
		epsilon = 1.0
		sigma = 1.0
		r6 = (sigma / r) ** 6
		r12 = r6 * r6
		force = 24 * epsilon * (2 * r12 - r6) / (r * r)
		return force

	def get_patch_particles(h_root):
		ret_value= []
		for h in IMP.atom.get_leaves(h_root):
			p= h.get_particle()
			if re.search("patch", p.get_name()):
				ret_value.append(p)
		return ret_value
	
	patch_particles= get_patch_particles(h_root)
	cpc= IMP.container.ClosePairContainer(patch_particles,
										  RANGE_PATCHES_A, # cutoff
										  100.0 # slack (affects only speed)
										  )
	lips= IMP.npctransport.LinearInteractionPairScore(K_EXCLUDED,
													  RANGE_PATCHES_A,
													  K_PATCHES)
	pr= IMP.container.PairsRestraint(lips, cpc) # score over dynamic container
	rs.append(pr)

	# Scoring Function from restraints
	sf = IMP.core.RestraintsScoringFunction(rs, "SF")
	#print(h_root.get_children())
	'''
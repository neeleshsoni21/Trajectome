"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class Interaction:
	def __init__(self, protein1, protein2, interaction_type, strength):
		self.protein1 = protein1
		self.protein2 = protein2
		self.interaction_type = interaction_type
		self.strength = strength

	def apply_interaction(self):
		# Implement interaction logic
		pass

	'''
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
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
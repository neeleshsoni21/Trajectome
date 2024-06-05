"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

import IMP
import IMP.atom
import IMP.rmf
import RMF
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class Simulation:
	def __init__(self, system, output_dir, time_steps, temperature):
		
		self.system = system
		self.time_steps = time_steps
		self.temperature = temperature
		self.model = system.model
		self.output_dir=output_dir
		self.bd = None
		self.setup_brownian_dynamics()

	def setup_brownian_dynamics(self):

		self.bd = IMP.atom.BrownianDynamics(self.model)
		print(self.system)
		scoring_function = IMP.core.RestraintsScoringFunction(self.system.restraints)
		self.bd.set_scoring_function(scoring_function)
		self.bd.set_time_step(10000)
		self.bd.set_temperature(self.temperature)
		self.bd.set_maximum_move(30)

		# RMF filename for trajectory
		RMF_FILENAME = self.output_dir+'/rigid.rmf'
		rmf = RMF.create_rmf_file(RMF_FILENAME)
		#hierarchy = IMP.atom.get_by_type(self.model, IMP.atom.MOLECULE_TYPE)
		IMP.rmf.add_hierarchy(rmf, self.system.h_root)
		IMP.rmf.add_restraints(rmf, self.system.restraints)
		IMP.rmf.add_geometry(rmf, IMP.display.BoundingBoxGeometry(self.system.bb))
		IMP.rmf.add_geometry(rmf, IMP.display.SphereGeometry(self.system.pbc_sphere))

		# Pair RMF with model using an OptimizerState ("listener")
		sos = IMP.rmf.SaveOptimizerState(self.model, rmf)
		sos.set_log_level(IMP.SILENT)
		sos.set_simulator(self.bd)
		self.bd.add_optimizer_state(sos)

		# Dump initial frame to RMF
		sos.update_always("initial conformation")

	def run(self):
		self.system.apply_boundary_conditions()  # Ensure boundary conditions are applied before running
		
		#for step in range(self.time_steps):
		#    self.perform_time_step(step)
		#    #time.sleep(0.1)  # Simulate time step duration
		
		#OR
		self.system.update()
		self.bd.optimize(self.time_steps)
		

	def perform_time_step(self, step):
		print(f"Performing time step {step}")
		self.system.update()
		self.bd.optimize(1000)
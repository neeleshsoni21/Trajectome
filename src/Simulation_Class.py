"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024

Attributes:
	logger (TYPE): Description
"""

import IMP
import IMP.atom
import IMP.rmf
import RMF
import logging


logger = logging.getLogger(__name__)

class Simulation:
	"""
	A class to represent and run simulations on a protein system.

	Attributes:
		bd (IMP.atom.BrownianDynamics): Brownian dynamics simulator.
		model (IMP.Model): The IMP model.
		output_dir (str): Directory for simulation output files.
		rmf_dump_interval_frames (int): Interval for dumping RMF frames.
		RMF_DUMP_INTERVAL_NS (float): Interval for dumping RMF frames in nanoseconds.
		scoring_function (IMP.core.RestraintsScoringFunction): Scoring function for restraints.
		sim_time_frames (int): Simulation time in frames.
		sim_time_ns (float): Simulation time in nanoseconds.
		simulation_time (float): Total simulation time.
		system (System): The system to be simulated.
		temperature (float): Temperature for the simulation.
	"""
	
	def __init__(self, system, output_dir, simulation_time, temperature):
		"""
		Initializes the Simulation instance with specified parameters.

		Args:
			system (System): The system to be simulated.
			output_dir (str): Directory for simulation output files.
			simulation_time (float): Total simulation time.
			temperature (float): Temperature for the simulation.

		Returns:
			None
		"""
		self.system = system
		self.simulation_time = simulation_time
		self.temperature = temperature
		self.model = system.model
		self.output_dir=output_dir
		self.bd = None
		self.setup_brownian_dynamics()

	def setup_brownian_dynamics(self):
		"""
		Sets up the Brownian dynamics simulation parameters.

		Args:
			None

		Returns:
			None
		"""
		# III. Time parameters:
		BD_STEP_SIZE_SEC= 10E-8
		SIM_TIME_SEC= self.simulation_time #0.0001 #0.050 #In seconds
		bd_step_size_fs= BD_STEP_SIZE_SEC * 1E+15
		self.sim_time_ns= SIM_TIME_SEC * 1E+9
		self.RMF_DUMP_INTERVAL_NS= self.sim_time_ns / 1000.0

		self.sim_time_frames= self.convert_time_ns_to_frames(self.sim_time_ns, bd_step_size_fs)
		self.rmf_dump_interval_frames= self.convert_time_ns_to_frames(self.RMF_DUMP_INTERVAL_NS, bd_step_size_fs)
		

		print("Simulation time {:.1e} ns / {} frames; "
		  "RMF dump interval {:.1e} ns / {} frames".format(self.sim_time_ns,
														  self.sim_time_frames,
														   self.RMF_DUMP_INTERVAL_NS,
														  self.rmf_dump_interval_frames))

		self.bd = IMP.atom.BrownianDynamics(self.model)
		#print(self.system)
		#print("Restraints:",self.system.restraints)
		#print([type(rs) for rs in self.system.restraints])
		self.scoring_function = IMP.core.RestraintsScoringFunction(self.system.restraints)
		self.bd.set_scoring_function(self.scoring_function)
		self.bd.set_time_step(10000)
		self.bd.set_temperature(self.temperature)
		self.bd.set_maximum_move(30)
		##self.bd.set_maximum_move(5)
		self.bd.set_maximum_time_step(bd_step_size_fs) # in femtoseconds



		# RMF filename for trajectory
		RMF_FILENAME = self.output_dir+'/rigid.rmf'
		rmf = RMF.create_rmf_file(RMF_FILENAME)
		#hierarchy = IMP.atom.get_by_type(self.model, IMP.atom.MOLECULE_TYPE)
		IMP.rmf.add_hierarchy(rmf, self.system.h_root)
		IMP.rmf.add_restraints(rmf, self.system.restraints)
		IMP.rmf.add_geometry(rmf, IMP.display.BoundingBoxGeometry(self.system.bb))
		IMP.rmf.add_geometry(rmf, IMP.display.BoundingBoxGeometry(self.system.bb_nucleus))
		IMP.rmf.add_geometry(rmf, IMP.display.BoundingBoxGeometry(self.system.bb_cytoplasm))
		#IMP.rmf.add_geometry(rmf, IMP.display.BoundingBoxGeometry(self.system.bb))
		

		# Pair RMF with model using an OptimizerState ("listener")
		sos = IMP.rmf.SaveOptimizerState(self.model, rmf)
		sos.set_log_level(IMP.SILENT)
		sos.set_simulator(self.bd)
		sos.set_period(self.rmf_dump_interval_frames)
		self.bd.add_optimizer_state(sos)

		# Dump initial frame to RMF
		sos.update_always("initial conformation")

	def run(self):
		"""
		Runs the simulation.

		Args:
			None

		Returns:
			None
		"""
		self.system.apply_boundary_conditions()  # Ensure boundary conditions are applied before running
		
		self.system.update()

		print("Running simulation")
		print("Score before: {:f}".format(self.scoring_function.evaluate(True)))
		self.bd.optimize(self.sim_time_frames)
		print("Run finished successfully")
		print("Score after: {:f}".format(self.scoring_function.evaluate(True)))

		

	def convert_time_ns_to_frames(self, time_ns, step_size_fs):
		"""
		Converts simulation time from nanoseconds to frames.

		Args:
			time_ns (float): Time in nanoseconds.
			step_size_fs (float): Step size in femtoseconds.

		Returns:
			int: Number of frames.
		"""
		FS_PER_NS= 1E6
		time_fs= time_ns * FS_PER_NS
		n_frames_float= (time_fs+0.0) / step_size_fs
		n_frames= int(round(n_frames_float))

		return max(n_frames, 1)
	
	
	
	
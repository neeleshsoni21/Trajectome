
"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024

This module defines the `Protein` class used for protein simulations.

Attributes:
	logger (logging.Logger): Logger for the module.
"""


import IMP
import IMP.core
import IMP.atom
import IMP.display
import IMP.pmi
import IMP.pmi.dof
import IMP.pmi.tools
import IMP.pmi.topology
import numpy as np
from collections import OrderedDict, defaultdict

import logging
import os


logger = logging.getLogger(__name__)

class Protein:
	"""
	A class to represent a protein in a simulation.

	Attributes:
		model (IMP.Model): The IMP model.
		name (str): Name of the protein.
		center (list): Coordinates of the protein's center.
		radius (int): Radius of the protein.
		mass (int): Mass of the protein.
		diffcoff (float): Diffusion coefficient of the protein.
		color (int): Color index for display.
		protp (IMP.Particle): The main particle representing the protein.
		hier (IMP.atom.Hierarchy): Hierarchy setup for the protein.
		prb (IMP.core.RigidBody): Rigid body of the protein.
		mol (IMP.atom.Molecule): Molecule representation of the protein.
		dif (IMP.atom.Diffusion): Diffusion setup for the protein.
	"""

	def __init__(self, model, name, center=[1,1,1], radius=30, mass=1000, diffcoff=0.01, color=0):
		"""
		Constructs all the necessary attributes for the Protein object.

		Args:
			model (IMP.Model): The IMP model.
			name (str): Name of the protein.
			center (list, optional): Coordinates of the protein's center. Default is [1, 1, 1].
			radius (int, optional): Radius of the protein. Default is 30.
			mass (int, optional): Mass of the protein. Default is 1000.
			diffcoff (float, optional): Diffusion coefficient of the protein. Default is 0.01.
			color (int, optional): Color index for display. Default is 0.
		"""
		self.model = model
		self.name = name

		self.protp = IMP.Particle(self.model, name)
		self.hier = IMP.atom.Hierarchy.setup_particle(self.protp)
		
		self.center = center
		self.radius = radius
		self.mass = mass
		self.diffcoff = diffcoff
		self.color = color

		self.create_rigid_body_protein()
		self.hier.add_child(self.mol)
		logger.info(f"Initialized Protein: {name} at center {center} with radius {radius}, mass {mass}, diffusion coefficient {diffcoff}, and color {color}.")


	def create_rigid_body_protein(self):
		"""
		Creates a rigid body representation of the protein.

		This method sets up the protein structure by combining residues into a rigid
		body and configuring necessary attributes such as diffusion coefficient.
		
		Args:
			None

		Returns:
			None
		"""
		self.prb = IMP.core.RigidBody.setup_particle(IMP.Particle(self.model), IMP.algebra.ReferenceFrame3D())
		self.prb.set_coordinates_are_optimized(True)
		self.prb.set_name(self.name + " rb")
		
		self.mol = IMP.atom.Molecule.setup_particle(IMP.Particle(self.model))
		self.mol.set_name(self.name)

		self.d = IMP.core.XYZR.setup_particle(
			IMP.Particle(self.model),
			IMP.algebra.Sphere3D(IMP.algebra.Vector3D(*self.center), self.radius)
		)
		IMP.display.Colored.setup_particle(self.d, IMP.display.get_display_color(self.color))
		
		self.mol.add_child(IMP.atom.Fragment.setup_particle(self.d))
		IMP.atom.Mass.setup_particle(self.d, self.mass)
		self.prb.add_member(self.d)
		self.d.set_name(self.name)
		
		self.dif = IMP.atom.Diffusion.setup_particle(self.prb, self.diffcoff)

		logger.info(f"Created rigid body protein {self.name}.")


class ProteinStructure:

	"""
	A class to represent and manipulate the structure of a protein for simulations.

	Attributes:
		All_Residues (dict): Dictionary containing all residues in the protein.
		amino_acid_radii (defaultdict): Default dictionary containing radii for different amino acids.
		Coarse_Residues (list): List of coarse-grained residues.
		color (int): Color index for display purposes.
		dif (IMP.atom.Diffusion): Diffusion setup for the protein.
		diffcoff (float): Diffusion coefficient of the protein.
		h_root (IMP.atom.Hierarchy): Root hierarchy of the protein.
		hier (IMP.atom.Hierarchy): Hierarchy setup for the protein.
		mass (float): Mass of the protein.
		model (IMP.Model): The IMP model.
		mol (IMP.atom.Molecule): Molecule representation of the protein.
		multi_model (bool): Indicates if the PDB file contains multiple models.
		name (str): Name of the protein.
		pdbfile (str): Path to the PDB file of the protein.
		prb (IMP.core.RigidBody): Rigid body of the protein.
		protein (IMP.atom.Protein): Protein representation.
		Protein_Residues (list): List of residues in the protein.
		protp (IMP.Particle): Main particle representing the protein.
		resolution (int): Resolution for the protein structure.
		Rot_diff (float): Rotational diffusion coefficient.
		state (IMP.pmi.topology.State): State of the protein in the simulation.
		Tran_dif (float): Translational diffusion coefficient.
	"""
	
	def __init__(self, model, state, h_root, name, pdbfile, pdb_multimodel= False, resolution= 50, diffcoff=0.01, color=0, centerize=False):
		"""
		Constructs all the necessary attributes for the ProteinStructure object.

		Args:
			model (IMP.Model): The IMP model.
			state (IMP.pmi.topology.State): State of the protein in the simulation.
			h_root (IMP.atom.Hierarchy): Root hierarchy of the protein.
			name (str): Name of the protein.
			pdbfile (str): Path to the PDB file of the protein.
			pdb_multimodel (bool, optional): Indicates if the PDB file contains multiple models. Default is False.
			resolution (int, optional): Resolution for the protein structure. Default is 50.
			diffcoff (float, optional): Diffusion coefficient of the protein. Default is 0.01.
			color (int, optional): Color index for display purposes. Default is 0.
			centerize (bool, optional): Whether to centerize the protein. Default is False.

		Returns:
			None
		"""
		def default_radius():
			"""
			Returns the default radius for amino acids.

			This method provides the default radius value used for amino acids when no
			specific radius is defined.
			
			Args:
				None

			Returns:
				float: The default radius value.
			"""
			return 4.0
		self.amino_acid_radii = defaultdict(default_radius,{
			"GLY": 3.4,"G": 3.4,
			"ALA": 3.8,"A": 3.8,
			"VAL": 4.0,"V": 4.0,
			"LEU": 4.0,"L": 4.0,
			"ILE": 4.0,"I": 4.0,
			"MET": 4.1,"M": 4.1,
			"PHE": 4.4,"F": 4.4,
			"TYR": 4.6,"Y": 4.6,
			"TRP": 4.8,"W": 4.8,
			"SER": 3.6,"S": 3.6,
			"THR": 3.8,"T": 3.8,
			"CYS": 3.8,"C": 3.8,
			"PRO": 4.0,"P": 4.0,
			"ASN": 3.8,"N": 3.8,
			"GLN": 4.0,"Q": 4.0,
			"ASP": 3.8,"D": 3.8,
			"GLU": 3.8,"E": 3.8,
			"LYS": 4.3,"K": 4.3,
			"ARG": 4.3,"R": 4.3,
			"HIS": 4.1,"H": 4.1
		})

		self.model = model
		self.name = name
		self.state = state
		self.h_root = h_root

		self.protp = IMP.Particle(self.model, name)
		self.hier = IMP.atom.Hierarchy.setup_particle(self.protp)

		self.pdbfile = pdbfile
		self.multi_model=pdb_multimodel

		self.resolution = resolution

		self.mass = 1.0
		self.diffcoff = diffcoff
		self.color = color

		self.Get_Hierarchy_From_PDB()

		self.create_rigid_body_protein()
		self.hier.add_child(self.mol)

		if centerize==True:
			self.translate_protein_to_center()

		logger.info(f"Initialized ProteinStructure: {name} from PDB file {pdbfile} with resolution {resolution}, diffusion coefficient {diffcoff}, and color {color}.")
	

		return

	def Get_Hierarchy_From_PDB(self):
		"""
		Retrieves the hierarchy of the protein from a PDB file.

		This method reads a PDB file and extracts the hierarchy of the protein structure
		for further processing and simulation.
		
		Args:
			pdbfile (str): The path to the PDB file.

		Returns:
			IMP.atom.Hierarchy: The hierarchy of the protein structure.
		"""
		#protein_hierarchies=None
		protein_hierarchy=None
		self.All_Residues =OrderedDict()

		splittup = os.path.splitext(self.pdbfile)
		if splittup[1] == ".pdb":

			if self.multi_model==True:
				protein_hierarchies = IMP.atom.read_multimodel_pdb(self.pdbfile, self.model)
				for prot_hierarchy in protein_hierarchies:
					residues = IMP.atom.get_by_type(prot_hierarchy, IMP.atom.RESIDUE_TYPE)
					self.All_Residues.append(residues)
			else:
				protein_hierarchy = IMP.atom.read_pdb(self.pdbfile, self.model)
				residues = IMP.atom.get_by_type(protein_hierarchy, IMP.atom.RESIDUE_TYPE)
				#self.All_Residues.append(residues)
				proteins = IMP.atom.get_by_type(protein_hierarchy, IMP.atom.CHAIN_TYPE)
				for prot in proteins:
					self.All_Residues['-'.join(prot.get_name().split(' '))]=prot.get_children()
		
		elif splittup[1] == ".cif":
			
			#For multi model hierarchy needs to be iterated
			if self.multi_model==True:
				protein_hierarchies = IMP.atom.read_multimodel_mmcif(self.pdbfile, self.model)
				for prot_hierarchy in protein_hierarchies:
					residues = IMP.atom.get_by_type(prot_hierarchy, IMP.atom.RESIDUE_TYPE)
					self.All_Residues.append(residues)
			else:
				protein_hierarchy = IMP.atom.read_mmcif(self.pdbfile, self.model)
				residues = IMP.atom.get_by_type(protein_hierarchy, IMP.atom.RESIDUE_TYPE)
				
				proteins = IMP.atom.get_by_type(protein_hierarchy, IMP.atom.CHAIN_TYPE)
				for prot in proteins:
					self.All_Residues['-'.join(prot.get_name().split(' '))]=prot.get_children()
					

				#TODO: Try selection method as they are giving actual prot name, instead of chain_0
				#mols = IMP.atom.get_by_type(protein_hierarchy, IMP.atom.MOLECULE_TYPE)
				#print("Mols:",mols)
				#print("RES:",list(self.All_Residues.values())[0:10])

				#residues =IMP.atom.Selection().get_selected_particles()
				#print("RES2:",residues[0:10])
				#sys.exit()
			
		else:
			raise ValueError(f"Unsupported file format: {self.pdbfile}")
			logger.error(f"Error reading file {self.pdbfile}")

		return

	def create_rigid_body_protein(self):
		"""
		Creates a rigid body representation of the protein.

		This method sets up the protein structure by combining residues into a rigid
		body and configuring necessary attributes such as diffusion coefficient.
		
		Args:
			None

		Returns:
			None
		"""
		self.protein = IMP.Particle(self.model)

		#Define a molecule within a rigid body. In this case it is a protein
		self.mol = IMP.atom.Molecule.setup_particle(self.protein)
		self.mol.set_name(self.name)
		##d = IMP.core.XYZ.setup_particle(self.mol.get_particle())#,IMP.algebra.ReferenceFrame3D())

		self.Protein_Residues =[]
		#for mdl_id, residues in enumerate(self.All_Residues):
		for mol_id, residues in self.All_Residues.items():
			self.setup_residues(mol_id, residues)

		#Create a rigid body for the molecule. We can add more molecules in this rigid body
		self.prb = IMP.core.RigidBody.setup_particle(IMP.Particle(self.model), self.Protein_Residues)
		self.prb.set_coordinates_are_optimized(True)
		self.prb.set_name(self.name + "_rb")

		Mol_radius_of_gyration = IMP.atom.get_radius_of_gyration(self.mol)
		IMP.core.XYZR(self.prb.get_particle()).set_radius(Mol_radius_of_gyration)
		
		self.Tran_dif = IMP.atom.Diffusion.setup_particle(self.prb, self.diffcoff)

		logger.info(str(self.name)+ " Radius final:"+str(IMP.core.XYZR(self.prb.get_particle())))
		
		rot_diff_scale=10
		if self.diffcoff==0:
			rot_diff_scale=0
		#Rotational diffusion ceofficient is calculated automatically considering radius. However, radius is not correctly identified
		self.Rot_diff = IMP.atom.RigidBodyDiffusion.setup_particle(self.prb)
		self.Rot_diff.set_rotational_diffusion_coefficient(
		 self.Rot_diff.get_rotational_diffusion_coefficient() * rot_diff_scale)

		logger.info("\nRot Diff coef radians2/fs"+str(self.Rot_diff.get_rotational_diffusion_coefficient()))
		logger.info("Trans Diff coef A^2/fs"+str(self.Tran_dif.get_diffusion_coefficient())+"\n")
		
		return

	def setup_residues(self, mol_id, residues):
		"""
		Sets up residues for the protein.

		This method initializes the residues and configures them with necessary
		attributes such as mass, diffusion, and color.
		
		Args:
			residues (list): A list of residues to be set up.

		Returns:
			None
		"""
		self.Coarse_Residues = []
		Temp_Frag=[]

		#TODO:CHECK FOR RESIDUE CONTINUATION
		for i, residue in enumerate(residues):
			if i % self.resolution==0:
				if len(Temp_Frag)!=0:
					self.combine_residues_to_fragment(Temp_Frag)
				Temp_Frag=[]
			
			#All atoms in the residue
			leaf1 = IMP.atom.get_leaves(residue)
			coordsR = IMP.core.XYZR(leaf1[0]) #Add first atom, mostly this will be N for full atom pdb file
			xyz = coordsR.get_coordinates()
			#R_rad = coordsR.get_radius() #Individual atom radius
			##if residue.get_name() not in self.amino_acid_radii.keys():
			##	print(residue)
			##	continue
			radius = self.amino_acid_radii[residue.get_name()] #Entire residue radius
			mass = self.mass
			residx = IMP.pmi.tools.get_residue_indexes(residue)[0]

			Temp_Frag.append([xyz,radius,mass,residx])

		#For final remaining particles
		if len(Temp_Frag)!=0:
			self.combine_residues_to_fragment(Temp_Frag)

		#self.Protein_Residues =[]
		for i, residue in enumerate(self.Coarse_Residues):
			
			xyz = residue[0]
			radius = residue[1]
			mass = residue[2]
			residx = residue[3]

			#Create a particle for the residue
			particle = IMP.Particle(self.model)
			d = IMP.core.XYZR.setup_particle(
				particle,
				IMP.algebra.Sphere3D(xyz, radius)
			)
			IMP.display.Colored.setup_particle(d, IMP.display.get_display_color(self.color))
			self.mol.add_child(IMP.atom.Fragment.setup_particle(d))
			IMP.atom.Mass.setup_particle(d, mass)
			self.Protein_Residues.append(d)

			d.set_name(f"{self.name}_{mol_id}_frag_{residx}")


		return

	def combine_residues_to_fragment(self, Temp_Frag):
		"""
		Combines residues into a single fragment.

		Args:
			residues (list): A list of residues to be combined into a fragment.

		Returns:
			None
		"""

		Mean_xyz = [0,0,0]
		Total_mass = 0
		ALL_residx = ""
		New_radius = 0

		for frag in Temp_Frag:
		
			Mean_xyz[0]+=frag[0][0];Mean_xyz[1]+=frag[0][1];Mean_xyz[2]+=frag[0][2]

			New_radius+=frag[1]**3

			Total_mass+=frag[2]

		Mean_xyz[0]= Mean_xyz[0]/len(Temp_Frag)
		Mean_xyz[1]= Mean_xyz[1]/len(Temp_Frag)
		Mean_xyz[2]= Mean_xyz[2]/len(Temp_Frag)

		#D=M1/V1=Mc/Vc
		#sum(Mi)=n*M1=Mc; sum(4/3*pi*ri^3) = 4/3*pi*rc^3
		#sum(ri^3) = rc^3
		#packing density 0.74
		New_radius = int((New_radius**(1/3))/0.74)

		ALL_residx = str(Temp_Frag[0][3])+'-'+str(Temp_Frag[-1][3])

		self.Coarse_Residues.append([Mean_xyz, New_radius, Total_mass, ALL_residx])

		return


	def create_rigid_body_protein_singlemodel(self):
		"""
		Creates a rigid body protein for a single model PDB file.

		This method sets up the protein structure from a single model PDB file and
		configures the necessary components such as residues and rigid bodies.
		
		Args:
			None

		Returns:
			None
		"""
		splittup = os.path.splitext(self.pdbfile)
		if splittup[1] == ".pdb":
			protein_hierarchy = IMP.atom.read_pdb(self.pdbfile, self.model)
		elif splittup[1] == ".cif":
			protein_hierarchy = IMP.atom.read_mmcif(self.pdbfile, self.model)

			#protein_hierarchy = IMP.atom.read_multimodel_mmcif(self.pdbfile, self.model)
		else:
			raise ValueError(f"Unsupported file format: {self.pdbfile}")
			logger.error(f"Error reading file {self.pdbfile}")

		#protein_hierarchy = IMP.atom.read_pdb(self.pdbfile, self.model)
		residues = IMP.atom.get_by_type(protein_hierarchy, IMP.atom.RESIDUE_TYPE)

		self.prb = IMP.core.RigidBody.setup_particle(IMP.Particle(self.model), IMP.algebra.ReferenceFrame3D())
		self.prb.set_coordinates_are_optimized(True)
		self.prb.set_name(self.name + " rb")
		
		self.mol = IMP.atom.Molecule.setup_particle(IMP.Particle(self.model))
		self.mol.set_name(self.name)

		for i, residue in enumerate(residues):
			
			particle = IMP.Particle(self.model)
			leaf1 = IMP.atom.get_leaves(residue)
			coordsR = IMP.core.XYZR(leaf1[0]) #Add first atom, mostly this will be N for full atom pdb file
			
			residx = IMP.pmi.tools.get_residue_indexes(residue)[0]
			xyz = coordsR.get_coordinates()
			radius = self.amino_acid_radii[residue.get_name()]

			d = IMP.core.XYZR.setup_particle(
				particle,
				IMP.algebra.Sphere3D(xyz, radius)
			)
			IMP.display.Colored.setup_particle(d, IMP.display.get_display_color(self.color))
			
			self.mol.add_child(IMP.atom.Fragment.setup_particle(d))
			IMP.atom.Mass.setup_particle(d, self.mass)
			self.prb.add_member(d)
			d.set_name(f"{self.name}_residue_{residx}")

		
		self.dif = IMP.atom.Diffusion.setup_particle(self.prb, self.diffcoff)

		return

	def translate_protein_to_center(self, center = [0, 0, 0]):
		"""
		Translates the protein to the specified center coordinates.

		Args:
			center (list, optional): New center coordinates. Default is [0, 0, 0].

		Returns:
			None
		"""
		coords_xyz = []
		xyzs = set();rbs = set()

		for p in IMP.atom.get_leaves(self.hier):
			coords_xyz.append(IMP.core.XYZ(p).get_coordinates())
			#xyzs.add(p)

		coords_mean = np.mean(coords_xyz, axis=0)
		translation_vector = [val2-val1 for val1,val2 in zip(coords_mean,center)]
		transformation = IMP.algebra.Transformation3D(IMP.algebra.Vector3D(translation_vector))
		
		for p in IMP.atom.get_leaves(self.hier):
			if IMP.core.RigidBody.get_is_setup(p):
				rbs.add(IMP.core.RigidBody(p))
			elif IMP.core.RigidMember.get_is_setup(p):
				rb = IMP.core.RigidMember(p).get_rigid_body()
				rbs.add(rb)
			else:
				xyzs.add(p)
		for xyz in xyzs:
			IMP.core.transform(IMP.core.XYZ(xyz), transformation)

		for rb in rbs:
			IMP.core.transform(rb, transformation)








"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
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

import logging
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class Protein:
	def __init__(self, model, name, center=[1,1,1], radius=30, mass=1000, diffcoff=0.01, color=0):
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

	def create_rigid_body_protein(self):
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

class ProteinStructure:
	def __init__(self, model, state, name, pdbfile, pdb_multimodel= False, resolution= 50, diffcoff=0.01, color=0, centerize=False):
		
		self.amino_acid_radii = {
			".": 4.0,
			"UNK": 4.0,
			"GLY": 3.4,
			"ALA": 3.8,
			"VAL": 4.0,
			"LEU": 4.0,
			"ILE": 4.0,
			"MET": 4.1,
			"PHE": 4.4,
			"TYR": 4.6,
			"TRP": 4.8,
			"SER": 3.6,
			"THR": 3.8,
			"CYS": 3.8,
			"PRO": 4.0,
			"ASN": 3.8,
			"GLN": 4.0,
			"ASP": 3.8,
			"GLU": 3.8,
			"LYS": 4.3,
			"ARG": 4.3,
			"HIS": 4.1
		}

		self.model = model
		self.name = name

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

		return

	def Get_Hierarchy_From_PDB(self):

		
		#protein_hierarchies=None
		protein_hierarchy=None
		self.All_Residues =[]

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
				self.All_Residues.append(residues)
		
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
				self.All_Residues.append(residues)
			
		else:
			raise ValueError(f"Unsupported file format: {self.pdbfile}")
			logger.error(f"Error reading file {self.pdbfile}")

		return

	def create_rigid_body_protein(self):
		
		self.protein = IMP.Particle(self.model)

		#Define a molecule within a rigid body. In this case it is a protein
		self.mol = IMP.atom.Molecule.setup_particle(self.protein)
		self.mol.set_name(self.name)

		for mdl_id, residues in enumerate(self.All_Residues):
			self.setup_residues(residues)

		#Create a rigid body for the molecule. We can add more molecules in this rigid body
		self.prb = IMP.core.RigidBody.setup_particle(IMP.Particle(self.model), self.Protein_Residues)
		self.prb.set_coordinates_are_optimized(True)
		self.prb.set_name(self.name + "_rb")

		Mol_radius_of_gyration = IMP.atom.get_radius_of_gyration(self.mol)
		IMP.core.XYZR(self.prb.get_particle()).set_radius(Mol_radius_of_gyration)
		
		self.Tran_dif = IMP.atom.Diffusion.setup_particle(self.prb, self.diffcoff)

		print(self.name, "Radius final:", IMP.core.XYZR(self.prb.get_particle()))
		
		rot_diff_scale=10
		if self.diffcoff==0:
			rot_diff_scale=0
		#Rotational diffusion ceofficient is calculated automatically considering radius. However, radius is not correctly identified
		self.Rot_diff = IMP.atom.RigidBodyDiffusion.setup_particle(self.prb)
		self.Rot_diff.set_rotational_diffusion_coefficient(
         self.Rot_diff.get_rotational_diffusion_coefficient() * rot_diff_scale)

		print("\nRot Diff coef radians2/fs",self.Rot_diff.get_rotational_diffusion_coefficient())
		print("Trans Diff coef A^2/fs",self.Tran_dif.get_diffusion_coefficient())
		print("\n")
		
		return

	def setup_residues(self, residues):

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
			radius = self.amino_acid_radii[residue.get_name()] #Entire residue radius
			mass = self.mass
			residx = IMP.pmi.tools.get_residue_indexes(residue)[0]

			Temp_Frag.append([xyz,radius,mass,residx])

		#For final remaining particles
		if len(Temp_Frag)!=0:
			self.combine_residues_to_fragment(Temp_Frag)

		self.Protein_Residues =[]
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

			d.set_name(f"{self.name}_fragment_{residx}")


		return

	def combine_residues_to_fragment(self, Temp_Frag):
		
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







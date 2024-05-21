"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

import IMP
import IMP.atom
from Protein_Class import Protein
from Interaction_Class import Interaction


class System:
    def __init__(self, L=3000, R=1300, K_BB=0.1):
        self.model = IMP.Model()
        self.p_root= IMP.Particle(self.model, "root")
        self.h_root = IMP.atom.Hierarchy.setup_particle(self.p_root)

        self.proteins = []
        self.interactions = []
        self.restraints = []
        self.L = L
        self.R = R
        self.K_BB = K_BB

        self.create_bounding_box_and_pbc()

    def create_bounding_box_and_pbc(self):
        # PBC cytoplasm bounding sphere:
        self.pbc_sphere = IMP.algebra.Sphere3D([0, 0, 0], self.R)

        # Outer bounding box for simulation:
        self.bb = IMP.algebra.BoundingBox3D(
            IMP.algebra.Vector3D(-self.L/2, -self.L/2, -self.L/2),
            IMP.algebra.Vector3D(self.L/2, self.L/2, self.L/2)
        )

        # Add enclosing spheres for pbc and outer simulation box
        self.bb_harmonic = IMP.core.HarmonicUpperBound(0, self.K_BB)
        self.pbc_bsss = IMP.core.BoundingSphere3DSingletonScore(self.bb_harmonic, self.pbc_sphere)
        self.outer_bbss = IMP.core.BoundingBox3DSingletonScore(self.bb_harmonic, self.bb)

    def add_protein(self, name, center, radius, mass, diffcoff, color):
        protein = Protein(self.model, name, center, radius, mass, diffcoff, color)
        self.proteins.append(protein)

    def add_interaction(self, protein1, protein2, interaction_type, strength):
        interaction = Interaction(protein1, protein2, interaction_type, strength)
        self.interactions.append(interaction)

    def add_restraint(self, restraint):
        self.restraints.append(restraint)

    def apply_boundary_conditions(self):
        rgbmembers = [protein.d.get_particle() for protein in self.proteins]
        # Add boundary conditions restraints
        self.add_restraint(IMP.container.SingletonsRestraint(self.pbc_bsss, IMP.container.ListSingletonContainer(self.model, rgbmembers)))
        self.add_restraint(IMP.container.SingletonsRestraint(self.outer_bbss, IMP.container.ListSingletonContainer(self.model, rgbmembers)))

    def add_excluded_volume_restraint(self):
        rgbmembers = [protein.d.get_particle() for protein in self.proteins]
        ev = IMP.core.ExcludedVolumeRestraint(IMP.container.ListSingletonContainer(self.model, rgbmembers), 1, 3)
        self.add_restraint(ev)

    def update(self):
        for interaction in self.interactions:
            interaction.apply_interaction()

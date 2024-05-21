
"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

import IMP
import IMP.core
import IMP.atom
import IMP.display

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
        self.hier.add_child(self.ph)

    def create_rigid_body_protein(self):
        
        self.prb = IMP.core.RigidBody.setup_particle(IMP.Particle(self.model), IMP.algebra.ReferenceFrame3D())
        self.prb.set_coordinates_are_optimized(True)
        self.prb.set_name(self.name + " rb")
        
        self.ph = IMP.atom.Molecule.setup_particle(IMP.Particle(self.model))
        self.ph.set_name(self.name)

        self.d = IMP.core.XYZR.setup_particle(
            IMP.Particle(self.model),
            IMP.algebra.Sphere3D(IMP.algebra.Vector3D(*self.center), self.radius)
        )
        IMP.display.Colored.setup_particle(self.d, IMP.display.get_display_color(self.color))
        
        self.ph.add_child(IMP.atom.Fragment.setup_particle(self.d))
        IMP.atom.Mass.setup_particle(self.d, self.mass)
        self.prb.add_member(self.d)
        self.d.set_name(self.name)
        
        self.dif = IMP.atom.Diffusion.setup_particle(self.prb, self.diffcoff)
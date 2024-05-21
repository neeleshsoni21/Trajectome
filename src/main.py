"""
Author: Neelesh Soni, neelesh@salilab.org, neeleshsoni03@gmail.com
Date: April 5, 2024
"""

from System_Class import System
from Simulation_Class import Simulation


if __name__ == "__main__":

    system = System()
    
    # Example: Adding proteins
    system.add_protein("Protein1", [1,1,1], 30, 1000, 0.01, 0)
    system.add_protein("Protein2", [2,2,2], 30, 1000, 0.01, 1)
    
    for prot in system.proteins:
        system.h_root.add_child(prot.hier)

    # Example: Adding interaction
    protein1 = system.proteins[0]
    protein2 = system.proteins[1]
    system.add_interaction(protein1, protein2, "binding", 10.0)
    
    # Adding Excluded Volume Restraint and boundary conditions
    system.add_excluded_volume_restraint()
    system.apply_boundary_conditions()

    # Running the simulation
    simulation = Simulation(system, output_dir='./output/',time_steps=1000, temperature=300)
    simulation.run()

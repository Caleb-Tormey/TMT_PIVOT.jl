# Main program for solving the Two Molecule Theory
#function Solver_Loop(Simulation_Parameters, Molecular_Information, Two_Molecule_Information)
function Solver_Loop(simulation_parameters::Simulation_Parameters, Molecule_1::Molecule_Information,Molecule_2::Molecule_Information)        
    # I. Initialize the simulation information arrays etc.
    # II. Initialize tweaking vector (need to addd the functions to do this)
    # III. Run solvation potential loop (Converge on W(r))
    #  ^  A. Pivot (using C(k) and Omega(k)) --> Omega(k)_new
    #  |  B. Run C(r) loop (Converge on C(r))
    #  |  ^   1. W(r)_new
    #  |  |   2. Direct Sampling (using W(r)_new) --> g(r)_simulation
    #  |  |   3. hr_fixed using g(r)_simulation and tweaking vector
    #  |  |   4. C(k)_new using C(k)_old + dC(k)
    #  |  |___5. Check for convergece of C(r) using C(r)_error using C(r)_new and C(r)_old
    #  |___C. Check for convergence on W(r) using W(r)_new and W(r)_old
    # Iv. Report Results
    #--------------------------------------------------------------------------------------#

    # I. Initialize the simulation information arrays etc.
    #--------------------------------------------------------------------------------------#
    #1) Initialize the Sigma and Epsilon static matrices for inter and intra molecular interactions for both molecules
        #a) InterMolecular
    simulation_parameters.Simulation_Interaction_Arrays.Sigma_Intermolecular, simulation_parameters.Simulation_Interaction_Arrays.Epsilon_Intermolecular = initialze_InterMolecular(Molecule_1,Molecule_2)
        #b) IntraMolecular
             #i) Molecule 1
    simulation_parameters.Simulation_Interaction_Arrays.Sigma_Intramolecular1, simulation_parameters.Simulation_Interaction_Arrays.Epsilon_Intramolecular1 = initialze_IntraMolecular(Molecule_1)
             #ii) Molecule 2 
    simulation_parameters.Simulation_Interaction_Arrays.Sigma_Intramolecular2, simulation_parameters.Simulation_Interaction_Arrays.Epsilon_Intramolecular2 = initialze_IntraMolecular(Molecule_2)
end

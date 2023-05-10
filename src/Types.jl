# define the data types for the simulation
#
#
# Simulation_Parameters will include the following:
# 1. Mixing parameter
# 2. Temperature 
# 3. Boltzmann constant
# 4. Equation and simulation grid space information.
#    a. Number of grid points (ngrid)
#    b. Grid spacing in r space (dr)
#    c. Grid spacing in k space (dk)
#    d. Discrete values in r space (rvalues)
#    e. Discrete values in k space (kvalues)
# 5. Solvation factor (solvFactor)
# 6. Direct sampling information
#    a. Number of direct sampling measurements (NSamples)
#    b. Sampling start index for sampling (SampleStart)
#    c. Sampling stop index for sampling (Sa1mpleEnd)
#    d. Number of threads to use (number_of_threads)
# 7. Pivot information
#    a. Do pivot (Do_Pivot(Bool)) this is incase of rigid molecules no need to do pivot. 
#    b. Number of configurations to save (NConfig)
#    c. Number of steps between saves (config_save_step)
#    d. Number of excluded neighbors in the intramolecular energy (exclude_neighbor)
struct Grid_Space
    ngrid::Int64
    dr::Float64
    dk::Float64
    #rvalues::SVector{Float64,1}
    rvalues::StaticArray
    #kvalues::SVector{Float64,1}
    kvalues::StaticArray
end

function Grid_Space(;
        ngrid = 1024,
        dr = 0.1,
        dk = pi/(dr*ngrid),
        rvalues = SVector{ngrid,Float64}(collect(0:dr:(ngrid-1)*dr)),
        kvalues = SVector{ngrid,Float64}(collect(0:dk:(ngrid-1)*dk))
    )
    return Grid_Space(ngrid,dr,dk,rvalues,kvalues)
end

struct Direct_Sampling_Information
    NSamples::Int64
    SampleStart::Int64
    SampleEnd::Int64
    number_of_threads::Int64
    Tail_Start::Int64
end

function Direct_Sampling_Information(;
        NSamples = 256*128,
        SampleStart = 30,
        SampleEnd = 600,
        number_of_threads = 4
    )
    return Direct_Sampling_Information(NSamples,SampleStart,SampleEnd,number_of_threads)
end

struct Pivot_Information
    Do_Pivot::Bool
    Do_WarmUp::Bool
    WarmUP_Steps::Int64
    NConfig::Int64
    config_save_step::Int64
    exclude_neighbor::Int64
    angle_range::Float64 
    dihedral_range::Float64
end

function Pivot_Information(;
        Do_Pivot = true,
        Do_WarmUp = true,
        WarmUP_Steps = 20000, 
        NConfig = 5000,
        config_save_step = 400,
        exclude_neighbor = 4,
        angle_range::Float64 = 20.0*pi/180.0,#20 degrees in radians
        dihedral_range::Float64 = pi/180.0 #90 degrees in radians
    )
    return Pivot_Information(Do_Pivot,Do_WarmUp,WarmUP_Steps,NConfig,config_save_step,exclude_neighbor,angle_range,dihedral_range)
end

mutable struct Simulation_Interaction_Arrays
    Sigma_Intermolecular::StaticArray
    Sigma_Intramolecular1::StaticArray
    Sigma_Intramolecular2::StaticArray
    Epsilon_Intermolecular::StaticArray
    Epsilon_Intramolecular1::StaticArray
    Epsilon_Intramolecular2::StaticArray
    Coulombic_Intermolecular::StaticArray
    Coulombic_Intramolecular1::StaticArray
    Coulombic_Intramolecular2::StaticArray
end

function Simulation_Interaction_Arrays(;
        Sigma_Intermolecular = SMatrix{24,24,Float64}(zeros(24,24)), 
        Sigma_Intramolecular1 = SMatrix{24,24,Float64}(zeros(24,24)), 
        Sigma_Intramolecular2 = SMatrix{24,24,Float64}(zeros(24,24)), 
        Epsilon_Intermolecular = SMatrix{24,24,Float64}(zeros(24,24)), 
        Epsilon_Intramolecular1= SMatrix{24,24,Float64}(zeros(24,24)),  
        Epsilon_Intramolecular2= SMatrix{24,24,Float64}(zeros(24,24)),  
        Coulombic_Intermolecular = SMatrix{24,24,Float64}(zeros(24,24)), 
        Coulombic_Intramolecular1 = SMatrix{24,24,Float64}(zeros(24,24)),
        Coulombic_Intramolecular2 = SMatrix{24,24,Float64}(zeros(24,24)),
        )
    return Simulation_Interaction_Arrays(Sigma_Intermolecular,Sigma_Intramolecular1,Sigma_Intramolecular2,Epsilon_Intermolecular,Epsilon_Intramolecular1,Epsilon_Intramolecular2,Coulombic_Intermolecular,Coulombic_Intramolecular1,Coulombic_Intramolecular2)
    
end

struct Simulation_Parameters
    Mixing_Parameter_Wr::Float64
    Mixing_Parameter_Cr::Float64
    Tolerance_Wr::Float64
    Tolerance_Cr::Float64
    Temperature::Float64
    Density::Float64
    kBoltz::Float64
    grid_space::Grid_Space
    solvFactor::Float64
    direct_sampling_information::Direct_Sampling_Information
    pivot_information::Pivot_Information
    Simulation_Interaction_Arrays::Simulation_Interaction_Arrays
    LJ_Strengh::Float64
    end

function Simulation_Parameters(;
        Mixing_Parameter_Wr = 0.25,
        Mixing_Parameter_Cr = 0.25,
        Tolerance_Wr = 1e-6,
        Tolerance_Cr = 1e-6,
        Temperature = 405.0, #K
        Density = 0.03123, #Angstrom^-3
        kBoltz = 0.001985875, #kcal/mol/K
        grid_space = Grid_Space(), 
        solvFactor = 1.0,
        direct_sampling_information = Direct_Sampling_Information(),
        pivot_information = Pivot_Information(),
        Simulation_Interaction_Arrays = Simulation_Interaction_Arrays(),
        LJ_Strengh = 1.1125 #This is the LJ strength for the a completely repulsve LJ potential
    )
    return Simulation_Parameters(Mixing_Parameter_Wr,Mixing_Parameter_Cr,Tolerance_Wr,Tolerance_Cr,Temperature,Density,kBoltz,grid_space,solvFactor,direct_sampling_information,pivot_information,Simulation_Interaction_Arrays,LJ_Strengh)
end

#Molecular_Information will include the following 
# 1. Molecule
#     a. Number of sites (Number_of_Sites)
#     b. Number of site types (Number_of_Types)
#     c. Topology of the molecule (Topology)
#         i. Molecule(Vector of Site types in order of the topology)
#         ii. Bonds (Bonds)
#         iii. Angles (Angles)
#         iv. Dihedrals (Dihedrals)
#         v. Impropers (Impropers)
#     d. IntraMolecular_Energy_Data
#         i. Bond_Length  (Bond_Energy)
#         ii. Angle Information (Angle_Energy)
#         iii. Dihedral Information (Dihedral_Energy)
#         iv. Improper Information (Improper_Energy)    
#     d. Array of site types (Site_Types)
#         i. Site information
#             1. Site Name(Site_Name)
#             2. Site Index(Site_Index)
#             3. Lennard Jones Parameters (LJ_Parameters)
#                a. Sigma (Sigma)
#                b. Epsilon (Epsilon)
#             4. Coulombic On (Coulombic_On)
#             5. Charge (Charge)
# 2. 
struct Site
    Site_Name::String
    Site_Index::Int64
    Sigma::Float64
    Epsilon::Float64
    Coulombic_On::Bool
    Charge::Float64
end


function Site(;
        Site_Name = "PE_A",
        Site_Index = 1,
        Sigma = 3.93,#Angstrom
        Epsilon = 0.09234,#kcal/mol
        Coulombic_On = false,
        Charge = 0.0
    )
    return Site(Site_Name,Site_Index,Sigma,Epsilon,Coulombic_On,Charge)

end

struct Topology
    Molecule::Vector
    Bonds::Vector{SArray}
    Angles::Vector{SArray}
    Dihedrals::Vector{SArray}
    Impropers::Vector{SArray}
end


function Topology(;
        Molecule = PE24(),
        Bonds = [SA[i,i+1] for i in 1:24-1],
        Angles = [SA[i,i+1,i+2] for i in 1:24-2],
        Dihedrals = [SA[i,i+1,i+2,i+3] for i in 1:24-3],
        Impropers = [SA[0,0,0,0]]
    )
    return Topology(Molecule,Bonds,Angles,Dihedrals,Impropers)
end
#This function will create the PE =24 molecule
function PE24()
    temp = Vector{}()
    for i in 1:24
        if i%2 != 0
            Site_Name = "PE_A"
            Site_Index = 1
        else
            Site_Name = "PE_B"
            Site_Index = 2
        end
        Sigma = 3.93
        Epsilon = 0.09234
        Coulombic_On = false
        Charge = 0.0
        site_temp= Site(Site_Name,Site_Index,Sigma,Epsilon,Coulombic_On,Charge)
        push!(temp,site_temp)
    end
    return temp
end
struct IntraMolecular_Energy_Data
    bondLength::Float64
    K_Bend::Float64
    Bend_Angle::Float64
    alist_Dihedral::SArray
    improper::Float64
end

function IntraMolecular_Energy_Data(;
        bondLength = 1.54,#Angstrom
        K_Bend = 124.18,#kcal/mol/rad^2
        Bend_Angle = 114.00,#degrees
        alist_Dihedral = SA[2.007,4.012,0.271,-6.290],#kcal/mol
        improper = 0.0
    )
    return IntraMolecular_Energy_Data(bondLength,K_Bend,Bend_Angle,alist_Dihedral,improper)
end

struct Molecule_Information
    Number_of_Sites::Int64
    Number_of_Types::Int64
    Topology::Topology
    Current_Coordinates::Vector{SArray}
    Population_Coordinates::Vector{SArray}
    Site_Types::Vector{Int64}
end

function Molecule_Information(;
        Number_of_Sites = 24,
        Number_of_Types = 2,
        Topology = Topology(),
        Current_Coordinates = [SA[0.0,0.0,0.0,1.0] for i in 1:24],
        Population_Coordinates = [SA[0.0,0.0,0.0,1.0] for i in 1:24],
        Site_Types = [1,2] 
    )
    return Molecule_Information(Number_of_Sites,Number_of_Types,Topology,Current_Coordinates,Population_Coordinates,Site_Types)
end


#Two_Molecule_Theory_Information will include the following
struct Two_Molecule_Data
    bigW::Matrix{Float64}
    bigW_Old::Matrix{Float64}
    bigWk::Matrix{Float64}
    gr_sim::Matrix{Float64}
    hr_sim::Matrix{Float64}
    hk_sim::Matrix{Float64}
    hr_fixed::Matrix{Float64}



end
function Two_Molecule_Data(;
    bigW = zeros(Float64,2,2),
    bigW_Old = zeros(Float64,2,2),
    bigWk = zeros(Float64,2,2),
    gr_sim = zeros(Float64,2,2),
    hr_sim = zeros(Float64,2,2),
    hk_sim = zeros(Float64,2,2),
    hr_fixed = zeros(Float64,2,2)
    )
    return Two_Molecule_Data(bigW,bigW_Old,bigWk,gr_sim,hr_sim,hk_sim,hr_fixed)
end

struct RISM_Data
    Cr::Matrix{Float64}
    Cr_Old::Matrix{Float64}
    Ck::Matrix{Float64}
    Ck_Old::Matrix{Float64}
    delta::Matrix{Float64}
    dCk::Matrix{Float64}
    CSC::Matrix{Float64}
    hr_RISM::Matrix{Float64}
    hk_RISM::Matrix{Float64}

end
function RISM_Data(;
    Cr = zeros(Float64,2,2),
    Cr_Old = zeros(Float64,2,2),
    Ck = zeros(Float64,2,2),
    Ck_Old = zeros(Float64,2,2),
    delta = zeros(Float64,2,2),
    dCk = zeros(Float64,2,2),
    CSC = zeros(Float64,2,2),
    hr_RISM = zeros(Float64,2,2),
    hk_RISM = zeros(Float64,2,2),
    )
    return RISM_Data(Cr,Cr_Old,Ck,Ck_Old,delta,dCk,CSC,hr_RISM,hk_RISM)
end

struct Theory_Data
    Two_Molecule_Data::Two_Molecule_Data    
    RISM_Data::RISM_Data
    Number_of_Sites::Int64
end

function Theory_Data(;
        Two_Molecule_Data = Two_Molecule_Data(),
        RISM_Data = RISM_Data()
        Number_of_Sites = 2
    )
    return Theory_Data(Two_Molecule_Data,RISM_Data,Number_of_Sites)
end

struct Sum_Rules
    Rules::Vector{Vector{Vector{Float64}}}
    Order::Vector{Float64}
end

function Sum_Rules(;
    Rules =[[[1, -1] [-1, 1]],[[1, 0] [-1, 0]],[[1, -1] [0, 0]],[[1, 0] [0,-1]]], #These rules are for a 2 site system that is nonlinear. 
    Order = [2,0,0,0]
    )
    return Sum_Rules(Rules,Order)
end


struct Simulation_Indexes
    Sim_Distance::Float64 #This is how far out to consider simulation data to enforce sum rules and is multiplied by the sigma value. 
    Begin_Sim_Correction::Matrix
    End_Sim_Correction::Matrix
end

mutable struct Correcting_Data
    A_Matrix::Matrix{Float64} #This is the matrix that will be constructed to enforce sum rules
    A_Matrix_Lengths::Matrix{Int64} # This matrix will hold how long each portion of the row vector for the fititng points is
    RHS::Vector{Float64} #This is the right hand side of the equation that will be constructed to enforce sum rules given the simulation results
    thing::Vector{Float64} #This is the vector that will be added to the simulation data over the rigth range to enforce sum rules. 
    Sum_Rules::Sum_Rules
    Simulation_Indexes::Simulation_Indexes #This will give the range of points that will have their values modified to enforce the rules (via thing vector)
end 

#=function Correcting_Data(;
    A_Matrix = zeros(Float64,2,2),
    RHS = zeros(Float64,2),
    thing = zeros(Float64,2),
    Sum_Rules = Sum_Rules(),
    Simulation_Correcting = Simulation_Correcting(),
    )
    return Correcting_Data(A_Matrix,RHS,thing,Sum_Rules)
end
=#

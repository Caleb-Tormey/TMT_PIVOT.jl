# Initialize the simulation information arrays etc. 
function initialze_IntraMolecular(Molecule::Molecule_Information)
    Sigma_Temp = zeros(Float64,Molecule.Number_of_Sites,Molecule.Number_of_Sites)
    Epsilon_Temp = zeros(Float64,Molecule.Number_of_Sites,Molecule.Number_of_Sites)
    for i in 1:Molecule.Number_of_Sites
        for j in 1:Molecule.Number_of_Sites
            Sigma_Temp[i,j] = berthelotS(Molecule.Topology.Molecule[i].Sigma,Molecule.Topology.Molecule[j].Sigma)
            Epsilon_Temp[i,j] = berthelotE(Molecule.Topology.Molecule[i].Epsilon,Molecule.Topology.Molecule[j].Epsilon)
        end
    end
    Sigma = SMatrix{Molecule.Number_of_Sites,Molecule.Number_of_Sites,Float64}(Sigma_Temp)
    Epsilon = SMatrix{Molecule.Number_of_Sites,Molecule.Number_of_Sites,Float64}(Epsilon_Temp)
    return Sigma,Epsilon
end



function initialze_InterMolecular(Molecule_1::Molecule_Information,Molecule_2::Molecule_Information)
    Sigma_Temp = zeros(Float64,Molecule_1.Number_of_Sites,Molecule_2.Number_of_Sites)
    Epsilon_Temp = zeros(Float64,Molecule_1.Number_of_Sites,Molecule_2.Number_of_Sites)
    for i in 1:Molecule_1.Number_of_Sites
        for j in 1:Molecule_2.Number_of_Sites
            Sigma_Temp[i,j] = berthelotS(Molecule_1.Topology.Molecule[i].Sigma,Molecule_2.Topology.Molecule[j].Sigma)
            Epsilon_Temp[i,j] = berthelotE(Molecule_1.Topology.Molecule[i].Epsilon,Molecule_2.Topology.Molecule[j].Epsilon)
        end
    end
    Sigma = SMatrix{Molecule_1.Number_of_Sites,Molecule_2.Number_of_Sites,Float64}(Sigma_Temp)
    Epsilon = SMatrix{Molecule_1.Number_of_Sites,Molecule_2.Number_of_Sites,Float64}(Epsilon_Temp)
    return Sigma,Epsilon
end

function Simulation_Indexes(Sim_Distance::Float64, Site_Sigma_List::Vector, grid_space::Grid_Space, theory_data::Theory_Data)
    #This needs to be fixed in the future to be more general. This means keeping track of the total number of site types in the system
    #Rather than just the molecules. Likely that is the best way to do it is come up with a System type that includes all the molecules in the
    #in the system of type Molecular_Information. 
     nos::Int64 = theory_data.Number_of_Sites
     Begin_Sim_Correction = zeros(Int,nos,nos)
     End_Sim_Correction = zeros(Int,nos,nos)
     for i in 1:length(Site_Sigma_List)
         for j in 1:length(Site_Sigma_List)
             Begin_Sim_Correction[i,j] = round(Int,(Site_Sigma_List[i] + Site_Sigma_List[j])/(2*grid_space.dr))
             End_Sim_Correction[i,j] = round(Int,Sim_Distance*(Site_Sigma_List[i] + Site_Sigma_List[j])/(2*grid_space.dr))
         end
     end
     return Simulation_Indexes(Sim_Distance,Begin_Sim_Correction,End_Sim_Correction)
 end
# Initialize A Matrix used for solving hr of r sum rules 

function initialize_A_Matrix(grid_space::Grid_Space, sum_rules::Sum_Rules,simulation_indexes::Simulation_Indexes,theory_data::Theory_Data)
    nos::Int64 = theory_data.Number_of_Sites
    total_rlist = simulation_indexes.End_Sim_Correction - simulation_indexes.Begin_Sim_Correction  .+ 1
    #rvec_len = sum(total_rlist)
    rvec_values = zeros(4,ngrid)
    for g in 1:length(sumr_rules.Order)
        for i in 1:grid_space.ngrid
            n = sum_rules.Order[g]
            rvec_values[g,i] = (real(im^n)*dr*4*pi)/(n+1)*grid_space.rvalues[i]^(n+2)
        end
    end
    #println(rvec_values[begin_sim[1,1]:end_sim[1,1]])
    A_vec_len = 0
    for i in 1:nos
        for j in 1:nos
            A_vec_len +=end_sim[i,j] - begin_sim[i,j] + 1
        end
    end
    A_matrix = zeros(length(sum_rules.Order),A_vec_len)
    for g in 1:length(sum_rules.Order)
        rvector = Array{Float64}(undef,0)
        for i in 1:nos
            for j in 1:nos
                temp = rvec_values[g,simulation_indexes.Begin_Sim_Correction[i,j]:simulation_indexes.End_Sim_Correction[i,j]].*sum_rules.Rules[g][i,j]
                rvector = vcat(rvector,temp)
            end
        end
        A_matrix[g,:] = rvector
    end
    return A_matrix,total_rlist
end

function Correcting_Data()
    #----NEED TO CLEAN THIS UP IN HERE AND IN THE DATATYPES---
    A_Matrix, A_Matrix_Lengths = initialize_A_Matrix(grid_space, sum_rules,simulation_correcting,theory_data)
    #Initialize with blank RHS and thing as these will be calculated later
    RHS = zeros(Float64,length(sumorder))
    thing = zeros(Float64,length(sumorder))
    Sum_Rules()
    Simulation_Correcting(Sim_Distance, Site_Sigma_List, grid_space, theory_data)
    return Correcting_Data(A_Matrix,A_Matrix_Lengths,RHS,thing,Sum_rules,Simulation_Correcting)
end
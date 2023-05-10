#Function used in calculating the energy of the system
#
function berthelotE(a,b)
    return sqrt(a*b)
end

function berthelotS(a,b)
    return 0.5*(a+b)
end
function fullLJ(sigma::Float64,epsilon::Float64,dist::Float64)
   return 4*epsilon*((sigma/dist)^12-(sigma/dist)^6)
end
function LJ_repulsive(sigma::Float64,epsilon::Float64,dist::Float64)
    return 4*epsilon*((sigma/dist)^12-(sigma/dist)^6 + 0.25)
end
function fullLJ(sigma::Float64,epsilon::Float64,dist::Float64,r_cut::Float64)
    if dist>r_cut
        return 0.0
    end
    return 4*epsilon*((sigma/dist)^12-(sigma/dist)^6+0.25)
end
function fullLJ(sigma::Float64,epsilon::Float64,dist::Float64,r_cut::Float64,lj_Shift::Float64)
    if dist>r_cut
        return 0.0
    end
    return 4*epsilon*((sigma/dist)^12-(sigma/dist)^6+0.25) - lj_Shift
end
function calc_distance(site1::SVector{4,Float64},site2::SVector{4,Float64})
    #return LinearAlgebra.norm(site1-site2)
    return sqrt((site1[1]-site2[1])^2+(site1[2]-site2[2])^2+(site1[3]-site2[3])^2)
end

function system_energy!(molecule1::Vector{SVector{4,Float64}},molN1::Int64,energyMol1::Float64,molecule2::Vector{SVector{4,Float64}},molN2::Int64,energyMol2::Float64,Sigma::SMatrix,Epsilon::SMatrix,index_temp::MMatrix,r_cut::Float64,LJ_shift::Float64,dr::Float64,solv_pot::Vector,TopMol1::SVector,TopMol2::SVector,ngrid::Int64,solFactor::Float64)
    total_energy::Float64 = 0.0
    solE_temp::Float64 = 0.0
    LJE_temp::Float64 = 0.0
    r_max::Float64 = ngrid*dr
    for i in 1:molN1
        for j in 1:molN2
            dist::Float64 = calc_distance(molecule1[i],molecule2[j])
            dist_2::Float64 = dist/dr
            dist_index::Int64 = round(Int,dist_2)
            rIndex::Float64 = ceil(dist_2)
            rFactor::Float64 = (rIndex-dist_2)
            LJE_temp += fullLJ(Sigma[i,j],Epsilon[i,j],dist,r_cut,LJ_shift)
            #println("i = "*string(i)*", j = "*string(j)*", distance = "*string(dist)*", energy = "*string(LJE_temp))
            if dist_index >=1
                index_temp[i,j] = dist_index
            elseif dist_index == 0
                dist_index = 1
                index_temp[i,j] = 1
            end
            if rIndex == 0.0
                solE_temp += solFactor*solv_pot[dist_index][TopMol1[i],TopMol2[j]]

            elseif rIndex > r_max
                solE_temp+=0.0

            else
                solE_temp +=solFactor*(solv_pot[dist_index][TopMol1[i],TopMol2[j]]*rFactor + (1 - rFactor)*solv_pot[dist_index + 1][TopMol1[i],TopMol2[j]])
            end

        end
    end
    #index_tempS = SMatrix{2,2,Int64}(index_temp)
    total_energy = solE_temp + LJE_temp
    # + energyMol1 + energyMol2


    return total_energy

end

function calcAngleE(molecule::Vector,kbend::Float64,theta::Float64)
    angleE::Float64 = 0.0
    n::Int64 = length(molecule)
    for i in 1:n - 2
        angle = calcAngle(molecule[i:i+2])
        angleE += angleEnergy(angle,kbend,theta)
    end
    return angleE
end
function calcDihedralE(molecule::Vector,alist::SVector)
    dihedralE::Float64 = 0.0
    n::Int64 = length(molecule)
    for i in 1:n - 3
        dh = calcDihedral(molecule[i:i+3])
        dihedralE += dihedralEnergy(dh,alist)
    end
    return dihedralE
end
function angleEnergy(angle::Float64, kbend::Float64, theta::Float64)
    return 0.5*kbend*(angle-theta)^2
end
function dihedralEnergy(dh::Float64,alist::SVector)
    dhE::Float64 = 0.0
    for i in 1:4
        dhE += alist[i]*(dh)^(i-1)
    end
    return dhE
end

function calcInternalLJE(dr::Float64,ngrid::Int64,exclude::Int64,molecule::Vector{SVector{4,Float64}},TopMol::SVector,Sigma_internal::SMatrix,Epsilon_internal::SMatrix,solv_pot::Vector,solFactor::Float64,r_cut::Float64,LJ_shift::Float64)
    total_internal_energy::Float64 = 0.0
    solE_internal_temp::Float64 = 0.0
    LJE_internal_temp::Float64 = 0.0
    r_max::Float64 = ngrid*dr
    n::Int64 = length(molecule)
    #count = 0
    #distTemp = 0.0
    #for i in 1:n-1
    #    for j in i+1:n
    #        dist = calc_distance(chain[i,1:3],chain[j,1:3])
            #remove this  because we will only do it ever 400 steps so no point copying data back and forth every time.
            #distArray[i,j] = dist
        #end
    #end
    for i in 1:n-exclude
        for j in i+exclude:n
            dist::Float64 = calc_distance(molecule[i],molecule[j])
            dist_2::Float64 = dist/dr
            dist_index::Int64 = round(Int,dist_2)
            rIndex::Float64 = ceil(dist_2)
            rFactor::Float64 = (rIndex-dist_2)
            LJE_internal_temp += fullLJ(Sigma_internal[i,j],Epsilon_internal[i,j],dist,r_cut,LJ_shift)
            if dist_index == 0
                dist_index = 1
            end
            if rIndex == 0.0
                solE_internal_temp += solFactor*solv_pot[dist_index][TopMol[i],TopMol[j]]

            elseif rIndex > r_max
                solE_internal_temp+=0.0

            else
                solE_internal_temp +=solFactor*(solv_pot[dist_index][TopMol[i],TopMol[j]]*rFactor + (1 - rFactor)*solv_pot[dist_index + 1][TopMol[i],TopMol[j]])
            end
            #dist = calcDistance(chain[i,1:3],chain[j,1:3])
            #r = distArray[i,j]
            #rindex = floor(Integer,r)
            #if rindex == 0
            #    solvpot = Wr[1]
            #else
            #    solvpot = (1-(r-rindex)/dr)*Wr[rindex]+((r-rindex)/dr)*Wr[rindex+1]
            #end
            #ljE += (fullLJ(sigma,epsilon,r)+solvpot)
            #distArray[i,j] = dist
            #count +=1
            #distTemp = dist
        end
    end
    total_internal_energy = solE_internal_temp + LJE_internal_temp
    return total_internal_energy
end

using DriftDiffusionPoissonSystems
using Distributions
using PyPlot


#* Source of randomness in u,v and V: discreteness and randomness of dopants
#* Semiconductor material (nanowire): Silicon (Si)
#* impurity atoms (dopants): Boron
#* background medium (electrolyte): water
#* insulator: Silicon dioxide (SiO2)
#* thickness of the oxide layer: 8nm
#* thickness of the nanowire: 40nm
#* length of the nanowire: 60nm
#* BCs:  
#         (a) Dirichlet BCS: Ve=0V (electrode voltage applied on the top of the device)
#                                      VG=-1V (back-gate voltage applied on the bottom of the device)
#
#          (b)Neuman BCS: everywhere else zero
#* Relative permittivities: ASi=11.7, Aox=3.9, Aliq=78, Adop=4.2
#* number of dopants: corresponds to the used doping concentration. For example
#       Ndop=6    ---> Cdop=5*10^15 cm-3
#       Ndop=600 ---> Cdop=4*10^17 cm-3
#* randomness: charge profile of dopants is modeled by Gaussian distributions:
#
#        Cdop(x) := \Sigma_j\frac{cj}{(2pi\sigma^2)^3/2}\exp{-\frac{(x-xj)^2}{2\sigma^2}},
#
#        where \Sigma_j= sum over j, \sigma=0.3nm influence parameter, cj=charge of jth atom, and          xj=position of the jth atom.
# * R=0, recombination rate
#*  UT= 26mV thermal voltage in the room temperature
#*  q= 1.6 *10^-19 Coulombs    elementary charge 
#*  ni= 1.5 * 10^10 cm-3 intrinsic charge of the semiconductor
#* \mun=1400, \mup=450 cm2/(V.s)
#*BCs:
#          (a) Dirichlet BCS: uD=ni^-1 e^(-V1/UT)nD
#                                       vD=ni^-1 e^(V1/UT)pD
#                where   nD= 1/2( Cdop + \sqrt{Cdop^2 + 4ni^2})
#                             pD= 1/2( -Cdop + \sqrt{Cdop^2 + 4ni^2}),
#and
#                            V1= U + UT \ln (nD/ni)
#with
#                            U= external applied voltage.
#          (b)Neuman BCS:  zero.

## Cdop_omega evaluates the function at (x,y) for a given nx2 matrix of dopant locations omega 
## sigma should be supplied in units of [nm]
## returns Cdop in [C/cm^3]
function Cdop(x,y, omega::Array{Float64}, sigma::Float64, signs::Array{Int64,1})

	q = 1.6021766*10.0^(-19)
	n = size(omega,1)
	c = zeros(size(x))
	## the factor 10^21 converts [nm^-3] to [cm^-3]
	f(x,y,w) = (10.0^21) * sqrt(2*pi)^(3/2) * (1/sigma)^3 *exp(-((x-w[1]).^2 + (y-w[2]).^2) ./ (2*sigma^2)) 

	for i in 1:n
		wi = omega[i,:]
		c += signs[i]* f(x,y,wi) 
	end
 	c
end

function point_in_triangle(p::Array{Float64}, p1::Array{Float64}, p2::Array{Float64}, p3::Array{Float64})
    # Determine whether p is a point in the triangle defined by nodes
    # p1, p2, p3
    epsilon = 0.001
    a, b, c, d = p2[1]-p1[1], p3[1]-p1[1], p2[2]-p1[2], p3[2]-p1[2]
    (s, t) = 1/(a*d-b*c)*[d -b; -c a]*[p[1] - p1[1], p[2] - p1[2]]
    if (0-epsilon<=s<=1+epsilon && 0-epsilon<=t<=1+epsilon && s+t <=1+epsilon)
        true
    else
        false
    end
end

function calculate_transistor_random_dopants(mesh::Mesh, omega::Array{Float64,2},signs::Array{Int64,1})
		
	q = 1.6021766*10.0^(-19) # elementary charge in C
	U_T = 0.026 # thermal voltage in V
	n_i = 1.5* 10^10 # intrinsic density in cm^-3
	mu_n = (x,y) -> 1400 # electron mobility in cm^2/Vs
	mu_p = (x,y) -> 450  # hole mobility in cm^2 /Vs
	eps_0 = 8.854*10.0^(-14) # vacuum permittivity in F/cm

	sigma = 0.3 # influence parameter in nm
	Endpoints = [0.0 0.0; 60.0 0.0; 60.0 40.0; 0.0 40.0]
	Vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> -1 (x,y) -> 0 (x,y) -> 0 (x,y) -> 0] #voltage of -1V at bottom


	C(x,y) = Cdop(x,y,omega,sigma,signs)
	nD(x,y) = 0.5*(C(x,y) + sqrt(C(x,y).^2+ 4n_i^2))
	pD(x,y) = 0.5*(-C(x,y) + sqrt(C(x,y).^2+ 4n_i^2))

	ubddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> e^(-Vbddata[3,1](x,y)/U_T) (x,y) -> 0 (x,y) -> e^(-Vbddata[3,3](x,y)/U_T) (x,y) -> 0]
	vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> ((1/n_i^2)* e^(Vbddata[3,1](x,y)/U_T) *pD(x,y))[1] (x,y) -> 0 (x,y) -> ((1/n_i^2)* e^(Vbddata[3,3](x,y)/U_T) * pD(x,y))[1] (x,y) -> 0]

	Asi = 11.7
	Adop = 4.2

	nq = size(mesh.nodes,2)
	nme = size(mesh.elements,2)

	#corner point indices
	a1 = vec(mesh.elements[1,:])
	a2 = vec(mesh.elements[2,:])
	a3 = vec(mesh.elements[3,:])

	#calculate coordinates of triangle centers
	centroid_x = (mesh.nodes[1,a1]+mesh.nodes[1,a2]+mesh.nodes[1,a3])./3
	centroid_y = (mesh.nodes[2,a1]+mesh.nodes[2,a2]+mesh.nodes[2,a3])./3

	vertex1(x) = mesh.nodes[:,mesh.elements[1,x]]
	vertex2(x) = mesh.nodes[:,mesh.elements[2,x]]
	vertex3(x) = mesh.nodes[:,mesh.elements[3,x]]

	#find indices of triangles that contain dopants
	ind_dopants = []
	for i in 1:size(omega,1)
		push!(ind_dopants, find(x->point_in_triangle(omega[i,:],vertex1(x),vertex2(x),vertex3(x)),1:size(mesh.elements,2)))
	end
	ind_dopants = collect(ind_dopants...)
	
	# returns permittivity of dopant if the point (x,y) is inside a triangle containing a dopant atom,
	# else returns the permittivity of silicon
	function epsilon(x::Float64,y::Float64)
		eps = 0
		this_triangle = find(z->point_in_triangle([x;y],vertex1(z),vertex2(z),vertex3(z)),1:size(mesh.elements,2))
		   if length(findin(ind_dopants,this_triangle)) != 0
			   eps = Adop*eps_0
	           else
			   eps = Asi*eps_0
	           end
		return eps
	end

	rec_mode = :zero # R = 0
	tau_p = 0.0
	tau_n = 0.0
	tol = 10.0^(-12)

	I,V,u,v=calculate_current(mesh,rec_mode,Endpoints,Vbddata,ubddata,vbddata,epsilon,U_T,n_i,tau_p,tau_n,mu_p,mu_n,C,tol)
		

end


function calculate_transistor(mesh::Mesh)
		
	q = 1.6021766*10.0^(-19) # elementary charge in C
	U_T = 0.026 # thermal voltage in V
	n_i = 1.5* 10^10 # intrinsic density in cm^-3
	mu_n = (x,y) -> 1400 # electron mobility in cm^2/Vs
	mu_p = (x,y) -> 450  # hole mobility in cm^2 /Vs
	eps_0 = 8.854*10.0^(-14) # vacuum permittivity in F/cm

	sigma = 0.3 # influence parameter in nm
	Endpoints = [0.0 0.0; 60.0 0.0; 60.0 40.0; 0.0 40.0]
	Vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> -1 (x,y) -> 0 (x,y) -> 0 (x,y) -> 0] #voltage of -1V at bottom

	function C(x,y)
		c = zeros(length(x))
		for i in 1:length(x)
			if x[i] < 20
				c[i] = 5*10^15 #Values correspond to 6 dopants
			elseif 20 <= x[i] <= 40
			        c[i] = -5*10^15
			elseif 40 < x[i]
			        c[i] = 5*10^15
			end
		end
		c
	end
		    
	nD(x,y) = 0.5*(C(x,y) + sqrt(C(x,y).^2+ 4n_i^2))
	pD(x,y) = 0.5*(-C(x,y) + sqrt(C(x,y).^2+ 4n_i^2))

	ubddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> e^(-Vbddata[3,1](x,y)/U_T) (x,y) -> 0 (x,y) -> e^(-Vbddata[3,3](x,y)/U_T) (x,y) -> 0]
	vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> ((1/n_i^2)* e^(Vbddata[3,1](x,y)/U_T) *pD(x,y))[1] (x,y) -> 0 (x,y) -> ((1/n_i^2)* e^(Vbddata[3,3](x,y)/U_T) * pD(x,y))[1] (x,y) -> 0]

	Asi = 11.7
	epsilon = (x,y) -> Asi*eps_0 

	rec_mode = :zero # R = 0
	tau_p = 0.0
	tau_n = 0.0
	tol = 10.0^(-12)

	I,V,u,v=calculate_current(mesh,rec_mode,Endpoints,Vbddata,ubddata,vbddata,epsilon,U_T,n_i,tau_p,tau_n,mu_p,mu_n,C,tol)
		

end
	






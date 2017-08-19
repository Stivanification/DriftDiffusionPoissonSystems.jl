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
## returns C_dop_scaled
function Cdop(x,y, omega::Array{Float64}, sigma::Float64, signs::Array{Int64,1})

	q = 1.6021766*10.0^(-19)
	n = size(omega,1)
	f(x,y,w) = exp(-((x-w[1]).^2 + (y-w[2]).^2) ./ (2*sigma^2)) 
	c = 0

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

function calculate_transistor_random_dopants(omega::Array{Float64,2},signs::Array{Int64,1})
		
	L = 6.0 * 10.0^(-6) # length of device in cm 
	H = 4.0 * 10.0^(-6) # height of device in cm (we assume L > H)

	V_bottom = -1.0 #external applied voltage at the bottom in Volt
	V_top = 0.0 #external applied voltage at the top in Volt

	## Mesh Generation
	Endpoints = [0.0 0.0 ; 1.0 0.0 ; 1.0 (H/L) ; 0.0 (H/L)] #mesh is scaled according to x_scaled * L = x
	h = 0.05
	#meshname = "mesh_h005H40"
	#construct_mesh4(Endpoints, h, meshname)
	#run(`gmsh -2 $meshname.geo`)

	mesh = read_mesh("mesh_h005H40.msh")

	## Physical Parameters
	q = 1.6021766*10.0^(-19) # elementary charge in C
	U_T = 0.0259 # thermal voltage in V
	n_i = 1.5* 10^10 # intrinsic density in cm^-3
	mu_n = (x,y) -> 1400 # electron mobility in cm^2/Vs
	mu_p = (x,y) -> 450  # hole mobility in cm^2 /Vs
	eps_0 = 8.854*10.0^(-14) # vacuum permittivity in F/cm
	A_Si = 11.7 # relative permittivity of Silicon
	sigma = 0.3 # influence parameter in nm

	C_scaled(x,y) = Cdop(x,y,omega,sigma,signs)

	## Scaling Parameters
	C_tilda = (10.0^21) * (sqrt(2*pi)^(3/2))*(1/sigma)^3 # the factor 10.0^21 converts [nm^-3] to [cm^-3] 
	mu_tilda = 1000 # reference mobility in cm^2 /Vs (used for scaling)
	V_tilda = max(abs(V_bottom),abs(V_top),U_T) # reference voltage used for scaling, V = V_tilda * V_scaled.

	lambda_0 = sqrt((eps_0 * V_tilda)/(q*C_tilda*L^2)) # scaling parameter lambda_0. We have lambda(x,y) = lambda_0 * sqrt(epsilon(x,y)) 
	delta = sqrt(n_i / C_tilda) # scaling parameter delta
	mu_n_scaled = (x,y) -> mu_n(x,y)/mu_tilda
	mu_p_scaled = (x,y) -> mu_p(x,y)/mu_tilda

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
	ind_dopants = [ind_dopants[i][1] for i in 1:length(ind_dopants)]	

	# returns permittivity of dopant if the point (x,y) is inside a triangle containing a dopant atom,
	# else returns the permittivity of silicon
	function epsilon(x::Float64,y::Float64)
		eps = 0
		this_triangle = find(z->point_in_triangle([x;y],vertex1(z),vertex2(z),vertex3(z)),1:size(mesh.elements,2))
		   if length(findin(ind_dopants,this_triangle)) != 0
			   eps = Adop
	           else
			   eps = Asi
	           end
		return eps
	end

	lambda = (x,y) -> lambda_0 * sqrt(epsilon(x,y))

	## We have V_D = U + V_bias, where U is the external applied voltage, and V_bias is given by
	## ln((1/(2*delta^2))*C_scaled + sqrt(C_scaled^2 + 4*delta^4))
	U_bottom = V_bottom / V_tilda   # the potential is scaled in multiples of V_tilda  
	U_top =  V_top / V_tilda 
	V_bias = (x,y) -> log((1/(2*delta^2)) * (C_scaled(x,y) + sqrt(C_scaled(x,y)^2 + 4*delta^4))) 

	Vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> U_bottom + V_bias(x,y)  (x,y) -> 0 (x,y) -> U_top + V_bias(x,y) (x,y) -> 0]
	# We have u_D = exp(-U(x,y)) and v_D = exp(U(x,y))
	ubddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> exp(-U_bottom) (x,y) -> 0 (x,y) -> exp(-U_top) (x,y) -> 0]
	vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> exp(U_bottom) (x,y) -> 0 (x,y) -> exp(U_top) (x,y) -> 0]

	rec_mode = :zero # R = 0
	tau_p = 0.0
	tau_n = 0.0
	tol = 10.0^(-12)

	I,V,u,v=calculate_current(mesh, rec_mode, Endpoints, Vbddata, ubddata, vbddata, lambda, delta, tau_p, tau_n, mu_p_scaled, mu_n_scaled,C_scaled,tol)

	I
end


function calculate_transistor()
	
	L = 6.0 * 10.0^(-6) # length of device in cm 
	H = 4.0 * 10.0^(-6) # height of device in cm (we assume L > H)

	V_bottom = -2.0 #external applied voltage at the bottom in Volt
	V_top = 0.0 #external applied voltage at the top in Volt

	## Mesh Generation
	Endpoints = [0.0 0.0 ; 1.0 0.0 ; 1.0 (H/L) ; 0.0 (H/L)] #mesh is scaled according to x_scaled * L = x
	h = 0.05
	#meshname = "mesh_h005H40"
	#construct_mesh4(Endpoints, h, meshname)
	#run(`gmsh -2 $meshname.geo`)

	mesh = read_mesh("mesh_h005H40.msh")

	## Physical Parameters
	q = 1.6021766*10.0^(-19) # elementary charge in C
	U_T = 0.0259 # thermal voltage in V
	n_i = 1.5* 10^10 # intrinsic density in cm^-3
	mu_n = (x,y) -> 1400 # electron mobility in cm^2/Vs
	mu_p = (x,y) -> 450  # hole mobility in cm^2 /Vs
	eps_0 = 8.854*10.0^(-14) # vacuum permittivity in F/cm
	A_Si = 11.7 # relative permittivity of Silicon

	## Scaling Parameters
	C_tilda = 5*10^15 # maximum of absolute value of density in cm^-3, used for scaling purposes
	mu_tilda = 1000 # reference mobility in cm^2 /Vs (used for scaling)
	V_tilda = max(abs(V_bottom),abs(V_top),U_T) # reference voltage used for scaling, V = V_tilda * V_scaled.

	lambda = (x,y) -> sqrt((A_Si * eps_0 * V_tilda)/(q*C_tilda*L^2)) # scaling parameter lambda 
	delta = sqrt(n_i / C_tilda) # scaling parameter delta
	mu_n_scaled = (x,y) -> mu_n(x,y)/mu_tilda
	mu_p_scaled = (x,y) -> mu_p(x,y)/mu_tilda

	## This doping function corresponds to a pnp doping, where each doping section is 20nm long
	function C_scaled(x,y) # C = C_tilda * C_scaled
			if (1/3)*(H/L) <= y <= (2/3)*(H/L)
			        c = -1
			else 
				c = 1
			end
			c
	end

	## We have V_D = U + V_bias, where U is the external applied voltage, and V_bias is given by
	## ln((1/(2*delta^2))*C_scaled + sqrt(C_scaled^2 + 4*delta^4))
	U_bottom = V_bottom / V_tilda   # the potential is scaled in multiples of V_tilda  
	U_top =  V_top / V_tilda 
	V_bias = (x,y) -> log((1/(2*delta^2)) * (C_scaled(x,y) + sqrt(C_scaled(x,y)^2 + 4*delta^4))) 

	Vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> U_bottom + V_bias(x,y)  (x,y) -> 0 (x,y) -> U_top + V_bias(x,y) (x,y) -> 0]
	# We have u_D = exp(-U(x,y)) and v_D = exp(U(x,y))
	ubddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> exp(-U_bottom) (x,y) -> 0 (x,y) -> exp(-U_top) (x,y) -> 0]
	vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> exp(U_bottom) (x,y) -> 0 (x,y) -> exp(U_top) (x,y) -> 0]

	rec_mode = :zero # R = 0
	tau_p = 0.0 #because R = 0, tau_n and tau_p can be set to arbitrary values
	tau_n = 0.0
	tol = 10.0^(-12)

	I,V,u,v=calculate_current(mesh, rec_mode, Endpoints, Vbddata, ubddata, vbddata, lambda, delta, tau_p, tau_n, mu_p_scaled, mu_n_scaled,C_scaled,tol)

	I = (q*V_tilda * C_tilda * mu_tilda)* I #rescales I to A/cm
 
end
	






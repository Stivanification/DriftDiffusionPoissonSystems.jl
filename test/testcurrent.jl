using PyPlot
using DriftDiffusionPoissonSystems

    function demoshockley()
	Endpoints = [-1.0 -1.0 ; 1.0 -1.0 ; 1.0 1.0 ; -1.0 1.0]
	Dirpoints = [-0.5 0.5; -0.5 0.5]
	h = 0.05
	meshname = "mesh_h$(h)_square8"
	construct_mesh8(Endpoints, Dirpoints, h, meshname)
	run(`gmsh -2 $meshname.geo`)

	mesh = read_mesh(meshname*".msh")
	c = 1.0

	Vbddata = [1 2 3 4 5 6 7 8; 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D'; (x,y)-> c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y);]
	ubddata = [1 2 3 4 5 6 7 8; 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D'; (x,y) ->exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y));]
	vbddata = [1 2 3 4 5 6 7 8; 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D'; (x,y) ->exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y));]

	U_T = 1.0
	tau_p = 1.0
	tau_n = 1.0
	mu_p = (x,y) -> 1.0
	mu_n = (x,y) -> 1.0
	n_i = 1.0
	epsilon = (x,y) -> 1.6021766*10.0^(-19)
	C = (x,y) -> 0.0
	tol= 10.0^(-12)
	rec_mode = :shockley

	I=calculate_current(mesh,rec_mode,Endpoints,Vbddata,ubddata,vbddata,epsilon,U_T,n_i,tau_p,tau_n,mu_p,mu_n,C, tol)
    end

function demozero()
	Endpoints = [-1.0 -1.0 ; 1.0 -1.0 ; 1.0 1.0 ; -1.0 1.0]
	h = 0.05
	meshname = "mesh_h$(h)_square4"
	construct_mesh4(Endpoints, h, meshname)
	run(`gmsh -2 $meshname.geo`)

	mesh = read_mesh(meshname*".msh")
	c = 1.0

	q = 1.6021766*10.0^(-19) # elementary charge in C
	U_T = 0.0259 # thermal voltage in V
	n_i = 1.5* 10^10 # intrinsic density in cm^-3
	mu_n = (x,y) -> 1400 # electron mobility in cm^2/Vs
	mu_p = (x,y) -> 450  # hole mobility in cm^2 /Vs
	eps_0 = 8.854*10.0^(-14) # vacuum permittivity in F/cm
	L = 6.0 * 10.0^(-6) # length of device in cm 
	H = 4.0 * 10.0^(-6) # height of device in cm (we assume L > H)
	C_tilda = 5*10^15 # maximum of absolute value of density in cm^-3
	mu_tilda = 1000 # reference mobility in cm^2 /Vs (used for scaling)
	lambda = (x,y) -> sqrt((eps_0 * U_T)/(q*C_tilda*L^2)) # scaling parameter lambda 
	delta = sqrt(n_i / C_tilda) # scaling parameter delta

	mu_n_scaled = (x,y) -> mu_n(x,y)/mu_tilda
	mu_p_scaled = (x,y) -> mu_p(x,y)/mu_tilda


	Vbddata = [1 2 3 4; 'D' 'D' 'D' 'D'; (x,y)-> c*(x^2-y^2) (x,y)->c*(x^2-y^2) (x,y)->c*(x^2-y^2) (x,y)->c*(x^2-y^2);]
	ubddata = [1 2 3 4; 'D' 'D' 'D' 'D'; (x,y) ->exp(-c*(x^2-y^2)) (x,y) -> exp(-c*(x^2-y^2)) (x,y) -> exp(-c*(x^2-y^2)) (x,y) -> exp(-c*(x^2-y^2));]
	vbddata = [1 2 3 4; 'D' 'D' 'D' 'D'; (x,y) ->exp(c*(x^2-y^2)) (x,y) -> exp(c*(x^2-y^2)) (x,y) -> exp(c*(x^2-y^2)) (x,y) -> exp(c*(x^2-y^2));]
	C = (x,y) -> 0.0
	tol= 10.0^(-12)
	tau_p = 1.0
	tau_n = 1.0
	rec_mode = :zero

	I,V,u,v=calculate_current(mesh, rec_mode, Endpoints, Vbddata, ubddata, vbddata, lambda, delta, tau_p, tau_n, mu_p, mu_n,C,tol)
	
	xcoord=zeros(size(mesh.nodes,2),1)
	ycoord=zeros(size(mesh.nodes,2),1)
	for i=1:size(mesh.nodes,2)
		xcoord[i]=mesh.nodes[1,i]
		ycoord[i]=mesh.nodes[2,i]
	end

	figure(1)
	f = (x,y) -> exp(c*(x.^2-y.^2))
	surf(vec(xcoord),vec(ycoord),vec(v),cmap="jet")
	exacterror = maximum(abs(vec(v) - f(xcoord,ycoord)))
	title("v(x,y)")
	println("maximum error v: $exacterror")

	figure(2)
	g = (x,y) -> x.^2-y.^2
	surf(vec(xcoord),vec(ycoord),vec(V),cmap="jet")
	exacterror = maximum(abs(vec(V) - g(xcoord,ycoord)))
	title("V(x,y)")
	println("maximum error V: $exacterror")

	I = (q*U_T * C_tilda * mu_tilda)* I #rescales I to A/cm
 
    end

		
	

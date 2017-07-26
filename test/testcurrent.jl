using DriftDiffusionPoissonSystems
using PyPlot

    function demoshockley()
	Endpoints = [-1.0 -1.0 ; 1.0 -1.0 ; 1.0 1.0 ; -1.0 1.0]
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
	epsilon = 1.6021766*10.0^(-19)
	C = (x,y) -> zeros(length(x))
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

	Vbddata = [1 2 3 4; 'D' 'D' 'D' 'D'; (x,y)-> c*(x^2-y^2) (x,y)->c*(x^2-y^2) (x,y)->c*(x^2-y^2) (x,y)->c*(x^2-y^2);]
	ubddata = [1 2 3 4; 'D' 'D' 'D' 'D'; (x,y) ->exp(-c*(x^2-y^2)) (x,y) -> exp(-c*(x^2-y^2)) (x,y) -> exp(-c*(x^2-y^2)) (x,y) -> exp(-c*(x^2-y^2));]
	vbddata = [1 2 3 4; 'D' 'D' 'D' 'D'; (x,y) ->exp(c*(x^2-y^2)) (x,y) -> exp(c*(x^2-y^2)) (x,y) -> exp(c*(x^2-y^2)) (x,y) -> exp(c*(x^2-y^2));]

	U_T = 1.0
	tau_p = 1.0
	tau_n = 1.0
	mu_p = (x,y) -> 1.0
	mu_n = (x,y) -> 1.0
	n_i = 1.0
	epsilon = (x,y) -> 1.6021766*10.0^(-19) * x.*y
	C = (x,y) -> zeros(length(x))
	tol= 10.0^(-12)
	rec_mode = :zero

	I,V,u,v=calculate_current(mesh,rec_mode,Endpoints,Vbddata,ubddata,vbddata,epsilon,U_T,n_i,tau_p,tau_n,mu_p,mu_n,C, tol)

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

	I
 
    end

		
	

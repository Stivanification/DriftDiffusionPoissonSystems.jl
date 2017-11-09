using DriftDiffusionPoissonSystems
using PyPlot
using Base.Test

function test_meshcreate()
	Endpoints = [-1.0 -1.0 ; 1.0 -1.0 ; 1.0 1.0 ; -1.0 1.0]
	Dirpoints = [-0.1 0.1 ; -0.1 0.1]
	h = 0.05
	meshname = "testmesh"
	construct_mesh(Endpoints, Dirpoints, h, meshname)
	run(`gmsh -2 $meshname.geo`)
end

function test_ddpe1()
	mesh = read_mesh("testmesh.msh")

	c = 1.0

	Vbddata = [1 2 3 4 5 6 7 8; 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D'; (x,y)-> c*(x+y) (x,y)->c*(x+y)	(x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y) (x,y)->c*(x+y);]
	ubddata = [1 2 3 4 5 6 7 8; 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D'; (x,y) ->exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y)) (x,y) -> exp(-c*(x+y));]
	vbddata = [1 2 3 4 5 6 7 8; 'D' 'D' 'D' 'D' 'D' 'D' 'D' 'D'; (x,y) ->exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y)) (x,y) -> exp(c*(x+y));]

	lambda = (x,y) -> 1.0
	delta = 1.0 
	tau_p = 0.0
	tau_n = 0.0
	mu_p = (x,y) -> 1.0
	mu_n = (x,y) -> 1.0
	rec_mode = :zero
	C = (x,y) -> 0.0
	tol= 10.0^(-12)

	tic()
	V,u,v=solve_ddpe(mesh,rec_mode,Vbddata,ubddata,vbddata,lambda,delta,tau_p,tau_n,mu_p,mu_n,C, tol)
	toc()

	xcoord=zeros(size(mesh.nodes,2),1)
	ycoord=zeros(size(mesh.nodes,2),1)
	for i=1:size(mesh.nodes,2)
		xcoord[i]=mesh.nodes[1,i]
		ycoord[i]=mesh.nodes[2,i]
	end

	f = (x,y) -> exp(c*(x+y))
	surf(vec(xcoord),vec(ycoord),vec(v),cmap="jet")
	exacterror = maximum(abs(vec(v) - f(xcoord,ycoord)))
	println("maximum error: $exacterror")
end

function test_ddpe2()
	mesh = read_mesh("testmesh.msh")

	c = 1.0

	Vbddata = [1 2 3 4 5 6 7 8; 'D' 'N' 'D' 'N' 'D' 'N' 'D' 'N'; (x,y)-> c (x,y)->0	(x,y)->c*(x+y) (x,y)->0 (x,y)->c*(x+y) (x,y)->0 (x,y)->c*(x+y) (x,y)->0;]
	ubddata = [1 2 3 4 5 6 7 8; 'D' 'N' 'D' 'N' 'D' 'N' 'D' 'N'; (x,y) ->exp(-c) (x,y) -> 0 (x,y) -> exp(-c)  (x,y) -> 0 (x,y) -> exp(-c)  (x,y) -> 0 (x,y) -> exp(-c)  (x,y) -> 0;]
	vbddata = [1 2 3 4 5 6 7 8; 'D' 'N' 'D' 'N' 'D' 'N' 'D' 'N'; (x,y) ->exp(c) (x,y) -> 0 (x,y) -> exp(c) (x,y) -> 0 (x,y) -> exp(c) (x,y) -> 0 (x,y) -> exp(c) (x,y) -> 0;]

	lambda = (x,y) -> 1.0
	delta = 1.0 
	tau_p = 1.0
	tau_n = 1.0
	mu_p = (x,y) -> 1.0
	mu_n = (x,y) -> 1.0
	rec_mode = :shockley
	C = (x,y) -> 0.0
	tol= 10.0^(-12)


	tic()
	V,u,v=solve_ddpe(mesh,rec_mode,Vbddata,ubddata,vbddata,lambda,delta,tau_p,tau_n,mu_p,mu_n,C, tol)
	toc()

	xcoord=zeros(size(mesh.nodes,2),1)
	ycoord=zeros(size(mesh.nodes,2),1)
	for i=1:size(mesh.nodes,2)
		xcoord[i]=mesh.nodes[1,i]
		ycoord[i]=mesh.nodes[2,i]
	end

	f = (x,y) -> exp(c)
	g = (x,y) -> exp(-c)
	surf(vec(xcoord),vec(ycoord),vec(v),cmap="jet")
	surf(vec(xcoord),vec(ycoord),vec(u),cmap="jet")

	exacterrorv = maximum(abs(vec(v) - f(xcoord,ycoord)))
	exacterroru = maximum(abs(vec(u) - g(xcoord,ycoord)))

	println("maximum error v: $exacterrorv")
	println("maximum error u: $exacterroru")
end

function test_semlin_poisson()
	mesh = read_mesh("mesh_p05.msh")
			f = (x,y,u) -> 4*exp(-x.^2 -y.^2).*(1 - x.^2 - y.^2) + sinh(exp(-x.^2-y.^2)) - sinh(u)
			g = (x,y,u) -> -cosh(u)
	    A = [(x) -> 1;(x)->0;(x)->0;(x)->1]
	bddata=[1 2 3 4;'D' 'N' 'D' 'N';(x,y)->exp(-x.^2-y.^2)  (x,y)->-2x.*exp(-x.^2-y.^2) (x,y)-> exp(-x.^2-y.^2) (x,y)->2x.*exp(-x.^2-y.^2)]
	tic()
	u=solve_semlin_poisson(mesh,A,bddata,f,g)
        toc()
	xcoord=zeros(size(mesh.nodes,2),1)
	ycoord=zeros(size(mesh.nodes,2),1)
	for i=1:size(mesh.nodes,2)
		xcoord[i]=mesh.nodes[1,i]
		ycoord[i]=mesh.nodes[2,i]
	end
        exakt = vec(exp(-xcoord.^2 -ycoord.^2))
        versuch = maximum(abs(vec(u)-exakt))
	println("error maximum:$versuch")
	surf(vec(xcoord),vec(ycoord),vec(u),cmap="jet")
end

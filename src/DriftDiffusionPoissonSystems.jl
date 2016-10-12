module DriftDiffusionPoissonSystems

using PyPlot

export Mesh,read_mesh,construct_mesh,solve_ddpe,solve_semlin_poisson

type Mesh
    #2xnrnodes array with [x;y] coords in each column
    nodes::Array{Float64, 2}

    #4xnredges array with [node1; node2; phys_prop; geom_prop] in each column
    edges::Array{Int64, 2}

    #5xnrelements array with [node1; node2; node3; phys_prop; geom_prop]
    elements::Array{Int64, 2}
end

###############################################################################
# construct_mesh(...) is a simple script that generates a .geo file to be processed with gmsh
# input: rectangle endpoints, endpoints of the dirichlet boundaries on the left and right side, triangle width h, meshname
# output: meshname.geo
function construct_mesh(Endpoints::Array, Dirpoints::Array, h::Float64, meshname::AbstractString)
	outfile = open("$meshname.geo","w")
	write(outfile, "h = $h;\n \n")

	write(outfile, "Point(1)={$(Endpoints[1,1]), $(Endpoints[1,2]), 0, h};\n")
	write(outfile, "Point(2)={$(Endpoints[2,1]), $(Endpoints[2,2]), 0, h};\n")
	write(outfile, "Point(3)={$(Endpoints[2,1]), $(Dirpoints[2,1]), 0, h};\n")
	write(outfile, "Point(4)={$(Endpoints[2,1]), $(Dirpoints[2,2]), 0, h};\n")
	write(outfile, "Point(5)={$(Endpoints[3,1]), $(Endpoints[3,2]), 0, h};\n")
	write(outfile, "Point(6)={$(Endpoints[4,1]), $(Endpoints[4,2]), 0, h};\n")
	write(outfile, "Point(7)={$(Endpoints[4,1]), $(Dirpoints[1,2]), 0, h};\n")
	write(outfile, "Point(8)={$(Endpoints[4,1]), $(Dirpoints[1,1]), 0, h};\n")

	write(outfile, "Line(1) = {1, 2};\n")
	write(outfile, "Line(2) = {2, 3};\n")
	write(outfile, "Line(3) = {3, 4};\n")
	write(outfile, "Line(4) = {4, 5};\n")
	write(outfile, "Line(5) = {5, 6};\n")
	write(outfile, "Line(6) = {6, 7};\n")
	write(outfile, "Line(7) = {7, 8};\n")
	write(outfile, "Line(8) = {8, 1};\n")

	write(outfile, "Line Loop(1) = {1,2,3,4,5,6,7,8};\n")
	write(outfile, "Plane Surface(2) ={1};\n")
	write(outfile, "Physical Surface(2) = {2};\n")
	#write(outfile, "Transfinite Surface {2}; \n")
	write(outfile, "Physical Line(1) = {1, 2, 3, 4, 5, 6, 7, 8};\n")
	#write(outfile, "Periodic Line{1} = {-3};\n")
	#write(outfile, "Periodic Line{4} = {-2};\n")
	write(outfile, "Coherence;")


	close(outfile)
end

###############################################################################
# read_mesh(...) is a function that converts a gmsh file to Type Mesh
# input: path to file
# output: Mesh
function read_mesh(path)

    in = open(path, "r")
    lines = readdlm(in)
    i = 1
    @assert lines[i] == "\$MeshFormat"
    i += 1
    s = lines[i]
    @assert s[1] == 2.2
    i += 1
    @assert lines[i] == "\$EndMeshFormat"
    i += 1

    if lines[i] == "\$PhysicalNames"
        i += 1
        nr_names = lines[i]
        i += nr_names + 2 #Skip names
    end

    #Save each node as type Node and collect nodes into nodes array
    @assert lines[i] == "\$Nodes"
    i += 1
    nr_nodes = lines[i]
    nodes = zeros(2, nr_nodes)
    @assert lines[i+1:i+nr_nodes,4]==vec(zeros(nr_nodes,1))
    nodes[:,:]=lines[i+1:i+nr_nodes,2:3]'
    i += nr_nodes+1
    @assert lines[i] == "\$EndNodes"
    i += 1

    #Save each element and edge
    @assert lines[i] == "\$Elements"
    i += 1
    nr_all_elements = lines[i]
    pos_edges = findin(lines[i+1:i+nr_all_elements,2],1)
    nr_edges = length(pos_edges)
    edges = Array{Int64}(4,nr_edges)
    pos_elements = findin(lines[i+1:i+nr_all_elements,2],2)
    nr_elements = length(pos_elements)
    elements = Array{Int64}(5,nr_elements)
    @assert lines[i+pos_edges,3] == vec(ones(nr_edges,1)*2)
    edges[:,:] = lines[i+pos_edges,[6,7,4,5]]'
    @assert lines[i+pos_elements,3] == vec(ones(nr_elements,1)*2)
    elements[:,:] = lines[i+pos_elements,[6,7,8,4,5]]'

    if nr_all_elements != nr_edges+nr_elements
            error("Unknown type of elements used! Please only use edges and triangular elements!")
    end

    i+=nr_all_elements+1
    @assert lines[i] == "\$EndElements"
    i += 1

    ## Read next section if available.
    if i <= size(lines,1)
    warn("There are lines in the mesh file that are not needed and cannot be taken into account.")
    end
    close(in)
    Mesh(nodes, edges, elements)
end


#############################################################
# aquire_boundary is a helping function that converts the boundary information into the form needed later.
# input: mesh, bddata
# output: x-coordinates of neumann-nodes (neu_nodes1),y-coordinates of neumann-nodes (neu_nodes2), indices of dirichlet nodes (n), values at dirichlet nodes (gD)

function aquire_boundary(mesh::Mesh, bddata::Array)
		### find out which BC applies where
		neu_bd=div(findin(bddata,'N'),3)+1
		diri_bd=div(findin(bddata,'D'),3)+1

		##check if all edges are covered
		if length(neu_bd)+length(diri_bd)!=size(bddata,2)
			error("Could not understand boundary file. Not all parts of the boundary have BC assigned.")
		end

		### Begin Neumann Part ###
		if neu_bd!=[]
			neu_nodes1=[]
			neu_nodes2=[]
			for j=1:size(neu_bd,1)
				neu_edges=findin(mesh.edges[4,:],bddata[1,neu_bd[j]])
				append!(neu_nodes1,vec(mesh.edges[1,neu_edges]))
				append!(neu_nodes2,vec(mesh.edges[2,neu_edges]))
			end
		end

		### Begin Dirichlet Part ###
		diri_nodes=[]
		diri_values=[]
		for j=1:size(diri_bd,1)
			diri_edges=findin(mesh.edges[4,:],bddata[1,diri_bd[j]])
			append!(diri_nodes,unique(vec(mesh.edges[1:2,diri_edges])))
			append!(diri_values,map(node->bddata[3,diri_bd[j]](mesh.nodes[1,node],mesh.nodes[2,node]),diri_nodes[end-length(diri_edges):end]))
		end

		diri_data=unique([diri_nodes diri_values],1)
		gD=map(Float64,diri_data[:,2])

		#n = boundary_nodes(mesh)
		n=diri_data[:,1]

	if neu_bd == []
		neu_nodes1 = []
		neu_nodes2 = []
	end

	neu_nodes1,neu_nodes2,gD,n
end


###################################################################
# MV_assemble assembles the mass matrix for the equation
# -\epsilon \Delta V = q n_i (u\cdot exp(V/U_T) - v\cdot exp(-V/U_T)) + qC(x)
# V = V_D on \Omega_D, \frac{\partial V}{\partial \nu} = 0 on \Omega_n
#
# input: mesh,n  (indices of dirichlet nodes), gD (values at dirichlet nodes), epsilon
# output: M (mass matrix for the potential equation), bdM (matrix containing boundary information), ar (vector of triangle areas)

function MV_assemble(mesh::Mesh, n::Array, gD::Array, epsilon::Float64)
	#data acquiry
    nq = size(mesh.nodes,2)
    nme = size(mesh.elements,2)

		ph_prop_11 = vec(epsilon.*ones(Float64,1,nme))
		ph_prop_22 = vec(epsilon.*ones(Float64,1,nme))

    #uu = 2xnr_elems array of node2-node3 with order [x;y]
    #vv = 2xnr_elems array of node3-node1
    #ww = 2xnr_elems array node1-node2
    #ar = determinant of hcat(ones(3), nodes_x, nodes_y)

    #corner point indices
    a1 = vec(mesh.elements[1,:])
    a2 = vec(mesh.elements[2,:])
    a3 = vec(mesh.elements[3,:])

    # generate stiffness matrix according to optim_2 algorithm of the paper
    q1 = mesh.nodes[:,a1]
    q2 = mesh.nodes[:,a2]
    q3 = mesh.nodes[:,a3]
    uu = q2-q3
    vv = q3-q1
    ww = q1-q2
    ar = (uu[1,:].*vv[2,:]-uu[2,:].*vv[1,:])./2
		ar = ar'

		#stiffness matrix assembly
		ar4 = abs(ar.*4)
		Kg = zeros(9,nme)
		ph_proper = [ph_prop_22 ph_prop_11;]'
		Kg[1,:] = sum(uu.*ph_proper.*uu,1)./ar4
		Kg[2,:] = sum(vv.*ph_proper.*uu,1)./ar4
		Kg[3,:] = sum(ww.*ph_proper.*uu,1)./ar4
		Kg[5,:] = sum(vv.*ph_proper.*vv,1)./ar4
		Kg[6,:] = sum(ww.*ph_proper.*vv,1)./ar4
		Kg[9,:] = sum(ww.*ph_proper.*ww,1)./ar4
		Kg[[4,7,8],:] = Kg[[2,3,6],:]
		Ig = vec(mesh.elements[[1, 2, 3, 1, 2, 3, 1, 2, 3],:])
		Jg = vec(mesh.elements[[1, 1, 1, 2, 2, 2, 3, 3, 3],:])
		Kg = vec(Kg[:])

	### Start changing stiffness matrix entries
    rows = findin(Ig, n)
    cols = findin(Jg, n)
    ### bdM gives the matrix necessary for subtraction from RHS
    bdM=sparse(vcat(Ig[rows],Ig[cols]),vcat(Jg[rows],Jg[cols]),vcat(Kg[rows],Kg[cols]),nq,nq)
    ### Change stiffness matrix entries
    Kg[rows] = 0
    Kg[cols] = 0
    append!(Ig, n)
    append!(Jg, n)
    append!(Kg, ones(Int64,length(n)))

    ### create matrix
    M = sparse(Ig[:], Jg[:], Kg[:], nq, nq)

	M,bdM,ar
end

###################################################################
# Vsolve solves the potential equation for fixed u0, v0 by employing Newton's method

function Vsolve(mesh::Mesh, bddata::Array, u::Array, v::Array, U_T::Float64, n_i::Float64, M::SparseMatrixCSC{Float64,Int64}, bdM::SparseMatrixCSC{Float64,Int64}, ar::Array, neu_nodes1::Array, neu_nodes2::Array, n::Array, gD::Array, C::Function, tol::Float64)

	#corner point indices
    a1 = vec(mesh.elements[1,:])
    a2 = vec(mesh.elements[2,:])
    a3 = vec(mesh.elements[3,:])

	nq = size(mesh.nodes,2)
    nme = size(mesh.elements,2)

	##	calculate positions of the triangle medians
    halfa_x=(mesh.nodes[1,a1]+mesh.nodes[1,a2])./2
    halfa_y=(mesh.nodes[2,a1]+mesh.nodes[2,a2])./2
    halfb_x=(mesh.nodes[1,a2]+mesh.nodes[1,a3])./2
    halfb_y=(mesh.nodes[2,a2]+mesh.nodes[2,a3])./2
    halfc_x=(mesh.nodes[1,a3]+mesh.nodes[1,a1])./2
    halfc_y=(mesh.nodes[2,a3]+mesh.nodes[2,a1])./2

	#initial guess V = 0
	V = transpose(zeros(Float64, 1, nq))

    V1 = reshape(V[a1],1,length(V[a1]))
    V2 = reshape(V[a2],1,length(V[a2]))
    V3 = reshape(V[a3],1,length(V[a3]))

	#Interpolate the solution values of the triangle nodes
    Sa = (V1 + V2)./(2*U_T)
    Sb = (V2 + V3)./(2*U_T)
    Sc = (V3 + V1)./(2*U_T)

	#Partition u0, v0 into parts belonging to the triangle nodes
    u1 = reshape(u[a1],1,length(u[a1]))
    u2 = reshape(u[a2],1,length(u[a2]))
    u3 = reshape(u[a3],1,length(u[a3]))

    v1 = reshape(v[a1],1,length(v[a1]))
    v2 = reshape(v[a2],1,length(v[a2]))
    v3 = reshape(v[a3],1,length(v[a3]))

	#interpolation of u0 and v0
	ua = (u1 + u2)./2
	ub = (u2 + u3)./2
	uc = (u3 + u1)./2

	va = (v1 + v2)./2
	vb = (v2 + v3)./2
	vc = (v3 + v1)./2

	#Right hand side vectors corresponding to a1, a2, a3
	fa=n_i*(exp(-Sa).*va - exp(Sa).*ua) + C(halfa_x, halfa_y)
	fb=n_i*(exp(-Sb).*vb - exp(Sb).*ub) + C(halfb_x, halfb_y)
	fc=n_i*(exp(-Sc).*vc - exp(Sc).*uc) + C(halfc_x, halfc_y)
    rhs=[fa fb fc]

    Irhs=vec([a1 a2 a3])
    Jrhs=vec(ones(Int64,length(Irhs)))
    rhsint=rhs.*abs([ar ar ar])./3 #triangle quadrature using medians
    Krhs=vec(rhsint)

	## Start boundary conditions
	append!(Irhs,vec([neu_nodes1;neu_nodes2]))
   	Jrhs=vec(ones(Int64,length(Irhs)))

	#Homogenous Neumann conditions
	neumanndata = zeros(Float64,1,size(neu_nodes1,1))
   	append!(Krhs,repmat(vec(neumanndata),2))

	## Subtract the columns from the right hand side.
    b=full(sparse(Irhs,Jrhs,Krhs,nq,1))
    b -= bdM[:, n] * gD
    b[n] = gD

    ## begin newton method on function (\varphi symbolises a vectorised function consisting of the \varphi_j)
    ## F(V) = M*V - \int f(x,y,V)* \varphi(x,y) dx dy

	while maximum(abs(M*V-b)) > tol

    	## The Jacobi matrix of the right hand side consists of the integrals \int g(x,y,V)*\varphi_i(x,y) *\varphi_j(x,y)
		## where g is the partial derivative of the right hand side with respect to V.
		## These integrals can be calculated using the algorithm to construct
    	## the weighted mass matrix WM according to the paper

			W1 = -n_i*(u1.*exp(V1) + v1.*exp(-V1)).*ar/30
			W2 = -n_i*(u2.*exp(V2) + v2.*exp(-V2)).*ar/30
			W3 = -n_i*(u3.*exp(V3) + v3.*exp(-V3)).*ar/30

    		WKg = zeros(Float64,9,nme)
    		WKg[1,:] = 3.*W1 + W2 + W3
    		WKg[2,:] = W1 + W2 + W3./2
    		WKg[3,:] = W1 + W2./2 + W3
    		WKg[5,:] = W1 + 3.*W2 + W3
    		WKg[6,:] = W1./2 + W2 + W3
    		WKg[9,:] = W1 + W2 + 3.*W3
    		WKg[[4,7,8],:] = WKg[[2,3,6],:]
			WIg = vec(mesh.elements[[1, 2, 3, 1, 2, 3, 1, 2, 3],:])
    		WJg = vec(mesh.elements[[1, 1, 1, 2, 2, 2, 3, 3, 3],:])
    		WKg = vec(WKg[:])

    		J0 = full(sparse(WIg[:], WJg[:], WKg[:], nq, nq))

			## newton iteration of V
 			DeltaV = -(M-J0)\(M*V-b)
    		V = DeltaV + V

			V1 = reshape(V[a1],1,length(V[a1]))
    		V2 = reshape(V[a2],1,length(V[a2]))
    		V3 = reshape(V[a3],1,length(V[a3]))

		    Sa = (V1 + V2)./(2*U_T)
			Sb = (V2 + V3)./(2*U_T)
			Sc = (V3 + V1)./(2*U_T)

			fa=n_i*(exp(-Sa).*va - exp(Sa).*ua) + C(halfa_x, halfa_y)
			fb=n_i*(exp(-Sb).*vb - exp(Sb).*ub) + C(halfb_x, halfb_y)
			fc=n_i*(exp(-Sc).*vc - exp(Sc).*uc) + C(halfc_x, halfc_y)
   			rhs=[fa fb fc]

    		rhsint=(rhs.*abs([ar ar ar]))./3 #triangle quadrature using medians
    		Krhs=vec(rhsint)

			if neu_nodes1 != []
				append!(Krhs,repmat(vec(neumanndata),2))
			end

    		## Subtract the columns from the right side.
    		b=full(sparse(Irhs,Jrhs,Krhs,nq,1))
    		b -= bdM[:, n] * gD
    		b[n] = gD

		control = maximum(abs(M*V - b))
   	end
    V
end

###################################################################
#   uvsolve solves the linear PDE
#	U_T * n_i div(\mu_n *exp(V/U_T) *grad(u)) = n_i * \frac{u*v0 -1}{\tau_p * (exp(V/U_T)*u + 1) + \tau_n * (exp(-V/U_T)*v0 + 1)}
#   for u, given v0.

function uvsolve(mesh::Mesh, V::Array, bddata::Array, U_T::Float64, tau_p::Float64, tau_n::Float64, mu::Function, neu_nodes1::Array, neu_nodes2::Array, n::Array, u0::Array, v0::Array)
	#data acquiry
	nq = size(mesh.nodes,2)
	nme = size(mesh.elements,2)

	#uu = 2xnr_elems array of node2-node3 with order [x;y]
	#vv = 2xnr_elems array of node3-node1
	#ww = 2xnr_elems array node1-node2
	#ar = determinant of hcat(ones(3), nodes_x, nodes_y)

	#corner point indices
	a1 = vec(mesh.elements[1,:])
	a2 = vec(mesh.elements[2,:])
	a3 = vec(mesh.elements[3,:])

	# generate stiffness matrix according to optim_2 algorithm of the paper
	q1 = mesh.nodes[:,a1]
	q2 = mesh.nodes[:,a2]
	q3 = mesh.nodes[:,a3]
	uu = q2-q3
	vv = q3-q1
	ww = q1-q2
	ar = (uu[1,:].*vv[2,:]-uu[2,:].*vv[1,:])./2
	ar = ar'

	V = V/U_T
	V1 = reshape(V[a1],1,length(V[a1]))
  V2 = reshape(V[a2],1,length(V[a2]))
  V3 = reshape(V[a3],1,length(V[a3]))

	#calculate coordinates of triangle centers
	centroid_x = (mesh.nodes[1,a1]+mesh.nodes[1,a2]+mesh.nodes[1,a3])./3
	centroid_y = (mesh.nodes[2,a1]+mesh.nodes[2,a2]+mesh.nodes[2,a3])./3

	#mu-values at the triangle centers
	mucentroid = vec(map(mu,centroid_x,centroid_y))
	# estimate permittivity on triangle by taking the permittivity value at the triangle center
	ph_prop = mucentroid.*vec(exp((V1+V2+V3)/3))

	ar4 = abs(ar.*4)
	Kg = zeros(9,nme)
	ph_proper = [ph_prop ph_prop;]'
	# println("uu: $(size(uu))")
	# println("ph_roper: $(size(ph_proper))")
	Kg[1,:] = sum(uu.*ph_proper.*uu,1)./ar4
	Kg[2,:] = sum(vv.*ph_proper.*uu,1)./ar4
	Kg[3,:] = sum(ww.*ph_proper.*uu,1)./ar4
	Kg[5,:] = sum(vv.*ph_proper.*vv,1)./ar4
	Kg[6,:] = sum(ww.*ph_proper.*vv,1)./ar4
	Kg[9,:] = sum(ww.*ph_proper.*ww,1)./ar4
	Kg[[4,7,8],:] = Kg[[2,3,6],:]
	Ig = vec(mesh.elements[[1, 2, 3, 1, 2, 3, 1, 2, 3],:])
	Jg = vec(mesh.elements[[1, 1, 1, 2, 2, 2, 3, 3, 3],:])
	Kg = vec(Kg[:])

    ## initial guess
    u = transpose(zeros(Float64, 1, nq))

    u1 = reshape(u[a1],1,length(u[a1]))
    u2 = reshape(u[a2],1,length(u[a2]))
    u3 = reshape(u[a3],1,length(u[a3]))

    v01 = reshape(v0[a1],1,length(v0[a1]))
    v02 = reshape(v0[a2],1,length(v0[a2]))
    v03 = reshape(v0[a3],1,length(v0[a3]))

    u01 = reshape(u0[a1],1,length(u0[a1]))
    u02 = reshape(u0[a2],1,length(u0[a2]))
    u03 = reshape(u0[a3],1,length(u0[a3]))

    ua = (u1 + u2)./2
    ub = (u2 + u3)./2
    uc = (u3 + u1)./2

	Va = (V1 + V2)./2
	Vb = (V2 + V3)./2
	Vc = (V3 + V1)./2

	v0a = (v01 + v02)./2
	v0b = (v02 + v03)./2
	v0c = (v03 + v01)./2

	u0a = (u01 + u02)./2
	u0b = (u02 + u03)./2
	u0c = (u03 + u01)./2

    fa=1./(tau_p*(exp(Va).*u0a + 1) + tau_n*(exp(-Va).*v0a + 1))
    fb=1./(tau_p*(exp(Vb).*u0b + 1) + tau_n*(exp(-Vb).*v0b + 1))
	fc=1./(tau_p*(exp(Vc).*u0c + 1) + tau_n*(exp(-Vc).*v0c + 1))
    rhs=[fa fb fc]

    Irhs=vec([a1 a2 a3])
    Jrhs=vec(ones(Int64,length(Irhs)))
    rhsint=rhs.*abs([ar ar ar])./3 #triangle quadrature using medians
    Krhs=vec(rhsint)

	### find out which BC applies where
		diri_bd=div(findin(bddata,'D'),3)+1

		### Begin Dirichlet Part ###
		diri_nodes=[]
		diri_values=[]
		for j=1:size(diri_bd,1)
			diri_edges=findin(mesh.edges[4,:],bddata[1,diri_bd[j]])
			append!(diri_nodes,unique(vec(mesh.edges[1:2,diri_edges])))
			append!(diri_values,map(node->bddata[3,diri_bd[j]](mesh.nodes[1,node],mesh.nodes[2,node]),diri_nodes[end-length(diri_edges):end]))
		end
		diri_data=unique([diri_nodes diri_values],1)
		gD=map(Float64,diri_data[:,2])

	### Start changing stiffness matrix entries
    rows = findin(Ig, n)
    cols = findin(Jg, n)
    ### bdM gives the matrix necessary for subtraction from RHS
    bdM=sparse(vcat(Ig[rows],Ig[cols]),vcat(Jg[rows],Jg[cols]),vcat(Kg[rows],Kg[cols]),nq,nq)
    ### Change stiffness matrix entries
    Kg[rows] = 0
    Kg[cols] = 0
    append!(Ig, n)
    append!(Jg, n)
    append!(Kg, ones(Int64,length(n)))

    ### create matrix
    M = sparse(Ig[:], Jg[:], Kg[:], nq, nq)

    ## Subtract the columns from the right side.
    b=full(sparse(Irhs,Jrhs,Krhs,nq,1))
    b -= bdM[:, n] * gD
    b[n] = gD

	# weighted mass matrix assembly
	W1 = -v01./(tau_p*(exp(V1).*u01+ 1) + tau_n*(exp(-V1).*v01 + 1)).*ar/30
	W2 = -v02./(tau_p*(exp(V2).*u02+ 1) + tau_n*(exp(-V2).*v02 + 1)).*ar/30
	W3 = -v03./(tau_p*(exp(V3).*u03+ 1) + tau_n*(exp(-V3).*v03 + 1)).*ar/30

	WKg = zeros(Float64,9,nme)
	WKg[1,:] = 3.*W1 + W2 + W3
	WKg[2,:] = W1 + W2 + W3./2
	WKg[3,:] = W1 + W2./2 + W3
	WKg[5,:] = W1 + 3.*W2 + W3
	WKg[6,:] = W1./2 + W2 + W3
	WKg[9,:] = W1 + W2 + 3.*W3
	WKg[[4,7,8],:] = WKg[[2,3,6],:]
	WIg = vec(mesh.elements[[1, 2, 3, 1, 2, 3, 1, 2, 3],:])
	WJg = vec(mesh.elements[[1, 1, 1, 2, 2, 2, 3, 3, 3],:])
	WKg = vec(WKg[:])

	J0 = full(sparse(WIg[:], WJg[:], WKg[:], nq, nq))

	u = (M-J0)\b
end


##################################################################
# G is the step function in the fixed point iteration used to solve the drift-diffusion poisson equations,
# a fixed point of G solves the system.
# recusively called helping function for solve_ddpe(...)

function G(mesh::Mesh, u0::Array, v0::Array, U_T::Float64, n_i::Float64, tau_p::Float64, tau_n::Float64, mu_p::Function, mu_n::Function, neu_nodes1::Array, neu_nodes2::Array, gD::Array, n::Array, MV::SparseMatrixCSC{Float64,Int64}, bdMV::SparseMatrixCSC{Float64,Int64}, ar::Array, Vbddata::Array, ubddata::Array, vbddata::Array, C::Function, tol::Float64)
	V = Vsolve(mesh, Vbddata, u0, v0, U_T, n_i, MV, bdMV, ar, neu_nodes1, neu_nodes2, n, gD, C, tol)
	u = uvsolve(mesh, V, ubddata, U_T, tau_p, tau_n, mu_p, neu_nodes1, neu_nodes2, n, u0, v0)
	v = uvsolve(mesh, -V, vbddata, U_T, tau_n, tau_p, mu_n, neu_nodes1, neu_nodes2, n, v0, u0)
 	u,v,V
end


###################################################################

#	solve_ddpe(...) solves a system of drift-diffusion poisson equations in slotboom variables
#
# 	\epsilon *\Delta V = q*n_i (exp(V/U_T)*u - exp(-V/U_T)*v) - q*C(x,y)
#	U_T * n_i div(\mu_n *exp(V/U_T) *grad(u)) = R
#	U_T * n_i div(\mu_p *exp(-V/U_T)*grad(v)) = R
#
#	given general Dirichlet boundary conditions mixed with homogenous Neumann
#	boundary conditions for (V,u,v).
#
#   R is the Shockley Read Hall recombination term
#	n_i * \frac{u*v -1}{\tau_p * (exp(V/U_T)*u + 1) + \tau_n * (exp(-V/U_T)*v + 1)}
#
#	the concentrations n and p can be recovered via
#	n = n_i * exp(V/U_T)*u, p = n_i * exp(-V/U_T)*v
#
#   solve_ddpe() solves the drift-diffusion poisson equations by fixed-point iteration,
#	which converges for small values of the potential.
#
#	input: * Mesh
#		   * constants epsilon, n_i, U_T, tau_p, tau_n (all Float64)
#		   * functions mu_p, mu_n, C
#		   * Arrays Vbddata, ubddata, vbddata (see readme)
#   output: potential V, slotboom variables u,v

function solve_ddpe(mesh::Mesh, Vbddata::Array, ubddata::Array, vbddata::Array, epsilon::Float64, U_T::Float64, n_i::Float64, tau_p::Float64, tau_n::Float64, mu_p::Function, mu_n::Function, C::Function, tol=10.0^(-9))

	nq = size(mesh.nodes,2)

	#initial guess u0 = v0 = 0
	u = zeros(Float64, nq,1)
	v = zeros(Float64, nq,1)
	V = zeros(Float64, nq,1)

	# elemental charge in Couloumb
	q = 1.6021766*10.0^(-19)

	#process boundary data for the potential equation
	neu_nodes1,neu_nodes2,gD,n = aquire_boundary(mesh,Vbddata)

	#assemble mass matrix for the potential equation
	MV,bdMV,ar = MV_assemble(mesh,n,gD,epsilon/q)

	control = 1
	count_dracula = 0

	#fixed point iteration of G
	while control > tol
		u0,v0,V0 = u,v,V
		u,v,V = G(mesh, u, v, U_T, n_i, tau_p, tau_n, mu_p, mu_n, neu_nodes1, neu_nodes2, gD, n, MV, bdMV, ar, Vbddata, ubddata, vbddata, C, tol)

		control = maximum([abs(u-u0) abs(v-v0) abs(V-V0)])
		count_dracula += 1
		println("iteration $count_dracula, tolerance: $control")
	end
    V,u,v
end

###################################################################
#
#	solve_semlin_poisson solves the semilinear Poisson equation
# 	- \nabla \cdot (A \nabla u) = f(x,y,u) with boundary conditions consisting of both Dirichlet and Neumann parts
# 	g := df/du, that is g is the partial derivative of the right hand side with respect to u

function solve_semlin_poisson(mesh, A::Array, bddata::Array, f::Function, g::Function, tol=10.0^(-14))
    #data acquiry
    nq = size(mesh.nodes,2)
    nme = size(mesh.elements,2)

    ph_prop_11 = zeros(Float64,1,nme)
    ph_prop_22 = zeros(Float64,1,nme)
    ph_prop_11=broadcast(A[1],float(mesh.elements[5,:]))
    ph_prop_22=broadcast(A[4],float(mesh.elements[5,:]))

    #uu = 2xnr_elems array of node2-node3 with order [x;y]
    #vv = 2xnr_elems array of node3-node1
    #ww = 2xnr_elems array node1-node2
    #ar = determinant of hcat(ones(3), nodes_x, nodes_y)

    #corner point indices
    a1 = vec(mesh.elements[1,:])
    a2 = vec(mesh.elements[2,:])
    a3 = vec(mesh.elements[3,:])

    # generate stiffness matrix according to optim_2 algorithm of the paper
    q1 = mesh.nodes[:,a1]
    q2 = mesh.nodes[:,a2]
    q3 = mesh.nodes[:,a3]
    uu = q2-q3
    vv = q3-q1
    ww = q1-q2
    ar = (uu[1,:].*vv[2,:]-uu[2,:].*vv[1,:])./2
		ar = ar'

		#stiffness matrix assembly
		ar4 = abs(ar.*4)
		Kg = zeros(9,nme)
		ph_proper = [ph_prop_22 ph_prop_11]'
		Kg[1,:] = sum(uu.*ph_proper.*uu,1)./ar4
		Kg[2,:] = sum(vv.*ph_proper.*uu,1)./ar4
		Kg[3,:] = sum(ww.*ph_proper.*uu,1)./ar4
		Kg[5,:] = sum(vv.*ph_proper.*vv,1)./ar4
		Kg[6,:] = sum(ww.*ph_proper.*vv,1)./ar4
		Kg[9,:] = sum(ww.*ph_proper.*ww,1)./ar4
		Kg[[4,7,8],:] = Kg[[2,3,6],:]
		Ig = vec(mesh.elements[[1, 2, 3, 1, 2, 3, 1, 2, 3],:])
		Jg = vec(mesh.elements[[1, 1, 1, 2, 2, 2, 3, 3, 3],:])
		Kg = vec(Kg[:])

		## calculate positions of the triangle medians

		halfa_x=(mesh.nodes[1,a1]+mesh.nodes[1,a2])./2
		halfa_y=(mesh.nodes[2,a1]+mesh.nodes[2,a2])./2
		halfb_x=(mesh.nodes[1,a2]+mesh.nodes[1,a3])./2
		halfb_y=(mesh.nodes[2,a2]+mesh.nodes[2,a3])./2
		halfc_x=(mesh.nodes[1,a3]+mesh.nodes[1,a1])./2
		halfc_y=(mesh.nodes[2,a3]+mesh.nodes[2,a1])./2

    ## initial guess V = 0
    V = transpose(zeros(Float64, 1, nq))

    V1 = reshape(V[a1],1,length(V[a1]))
    V2 = reshape(V[a2],1,length(V[a2]))
    V3 = reshape(V[a3],1,length(V[a3]))

		## start with intial guess for the right hand side integrals
		fa=transpose(vec(f(vec(halfa_x),vec(halfa_y), 0)))
		fb=transpose(vec(f(vec(halfb_x),vec(halfb_y), 0)))
		fc=transpose(vec(f(vec(halfc_x),vec(halfc_y), 0)))
	  rhs=[fa fb fc]

    Irhs=vec([a1 a2 a3])
    Jrhs=vec(ones(Int64,length(Irhs)))
    rhsint=rhs.*abs([ar ar ar])./3 #triangle quadrature using medians
    Krhs=vec(rhsint)


	###start boundary conditions

	### find out which BC applies where
	neu_bd=div(findin(bddata,'N'),3)+1
	diri_bd=div(findin(bddata,'D'),3)+1
	##check if all edges are covered
	if length(neu_bd)+length(diri_bd)!=size(bddata,2)
		error("Could not understand boundary file. Not all parts of the boundary have BC assigned.")
	end

	### Begin Neumann Part ###
	if neu_bd!=[]
   		neu_nodes1=[]
   		neu_nodes2=[]
   		neu_fct=[]
   		for j=1:size(neu_bd,1)
   			neu_edges=findin(mesh.edges[4,:],bddata[1,neu_bd[j]])
   			append!(neu_nodes1,vec(mesh.edges[1,neu_edges]))
   			append!(neu_nodes2,vec(mesh.edges[2,neu_edges]))
   			append!(neu_fct,vec(broadcast(bddata[3,neu_bd[j]],(mesh.nodes[1,vec(mesh.edges[1,neu_edges])]+mesh.nodes[1,vec(mesh.edges[2,neu_edges])])/2,(mesh.nodes[2,vec(mesh.edges[1,neu_edges])]+mesh.nodes[2,vec(mesh.edges[2,neu_edges])])/2)))
   		end
   		### The flow is given by the surface integral over the respective edges.
   		### We use the formula int_dE gn*eta ds ~ |E|/2*g(center_x,center_y)

   		##Calculate the flows at the edges.
   		neumanndata=sqrt((mesh.nodes[1,neu_nodes1]-mesh.nodes[1,neu_nodes2]).^2+(mesh.nodes[2,neu_nodes1]-mesh.nodes[2,neu_nodes2]).^2).*neu_fct/2
			neumanndata = vec(repmat(neumanndata,2))

    	#The vector [neu_nodes1;neu_nodes2] together with the repeated neumanndata adds all the Neumann data to the correct nodes.
  		append!(Irhs,vec([neu_nodes1;neu_nodes2]))
   		Jrhs=vec(ones(Int64,length(Irhs)))
   		append!(Krhs,neumanndata)

	end

	### Begin Dirichlet Part ###
	diri_nodes=[]
	diri_values=[]
	for j=1:size(diri_bd,1)
		diri_edges=findin(mesh.edges[4,:],bddata[1,diri_bd[j]])
		append!(diri_nodes,unique(vec(mesh.edges[1:2,diri_edges])))
       	append!(diri_values,map(node->bddata[3,diri_bd[j]](mesh.nodes[1,node],mesh.nodes[2,node]),diri_nodes[end-length(diri_edges):end]))
    end
    diri_data=unique([diri_nodes diri_values],1)
    gD=map(Float64,diri_data[:,2])


    #n = boundary_nodes(mesh)
    n=diri_data[:,1]

		## Start changing stiffness matrix entries
    rows = findin(Ig, n)
    cols = findin(Jg, n)
    ### bdM gives the matrix necessary for subtraction from RHS
    bdM=sparse(vcat(Ig[rows],Ig[cols]),vcat(Jg[rows],Jg[cols]),vcat(Kg[rows],Kg[cols]),nq,nq)
    ### Change stiffness matrix entries
    Kg[rows] = 0
    Kg[cols] = 0
    append!(Ig, n)
    append!(Jg, n)
    append!(Kg, ones(Int64,length(n)))

    ### create matrix
    M = sparse(Ig[:], Jg[:], Kg[:], nq, nq)

    ## Subtract the columns from the right side.
    b=full(sparse(Irhs,Jrhs,Krhs,nq,1))
    b -= bdM[:, n] * gD
    b[n] = gD

    ## begin newton method on function (\varphi symbolises a vectorised function consisting of the \varphi_j)
    ## F(V) = M*V - \int f(x,y,V)* \varphi(x,y) dx dy

    count = 0

	while maximum(abs(M*V-b)) > tol

    	## The Jacobi matrix of the right hand side consists of the integrals \int g(x,y,V)*\varphi_i(x,y) *\varphi_j(x,y).
    	## They can be calculated using the algorithm to construct
    	## the weighted mass matrix WM according to the paper

  		W1 = g(mesh.nodes[1,a1], mesh.nodes[2,a1], V1).*ar/30
  		W2 = g(mesh.nodes[1,a2], mesh.nodes[2,a2], V2).*ar/30
  		W3 = g(mesh.nodes[1,a3], mesh.nodes[2,a3], V3).*ar/30

			WKg = zeros(Float64,9,nme)
  		WKg[1,:] = 3.*W1 + W2 + W3
  		WKg[2,:] = W1 + W2 + W3./2
  		WKg[3,:] = W1 + W2./2 + W3
  		WKg[5,:] = W1 + 3.*W2 + W3
  		WKg[6,:] = W1./2 + W2 + W3
  		WKg[9,:] = W1 + W2 + 3.*W3
  		WKg[[4,7,8],:] = WKg[[2,3,6],:]
			WIg = vec(mesh.elements[[1, 2, 3, 1, 2, 3, 1, 2, 3],:])
  		WJg = vec(mesh.elements[[1, 1, 1, 2, 2, 2, 3, 3, 3],:])
  		WKg = vec(WKg[:])

  		J0 = full(sparse(WIg[:], WJg[:], WKg[:], nq, nq))
 			DeltaV = -(M-J0)\(M*V-b)
  		V += DeltaV

			V1 = reshape(V[a1],1,length(V[a1]))
  		V2 = reshape(V[a2],1,length(V[a2]))
  		V3 = reshape(V[a3],1,length(V[a3]))

			Ca = (V1 + V2)./2
			Cb = (V2 + V3)./2
			Cc = (V3 + V1)./2

			fa=transpose(vec(f(vec(halfa_x),vec(halfa_y), vec(Ca))))
			fb=transpose(vec(f(vec(halfb_x),vec(halfb_y), vec(Cb))))
			fc=transpose(vec(f(vec(halfc_x),vec(halfc_y), vec(Cc))))
		  rhs=[fa fb fc]

  		rhsint=(rhs.*abs([ar ar ar]))./3 #triangle quadrature using medians
  		Krhs=vec(rhsint)

			if neu_bd != []
				append!(Krhs,neumanndata)
	    end

  		## Subtract the columns from the right side.
  		b=full(sparse(Irhs,Jrhs,Krhs,nq,1))
  		b -= bdM[:, n] * gD
  		b[n] = gD

			count += 1
			control = maximum(abs(M*V - b))
			println("iteration $count, tolerance: $control")
 	end
  V
end

end#module

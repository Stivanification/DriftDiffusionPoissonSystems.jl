
mesh = read_mesh("mesh_p05.msh")
f = (x,y,u) -> 4*exp(-x.^2 -y.^2).*(1 - x.^2 - y.^2) + sinh(exp(-x.^2-y.^2)) - sinh(u)
g = (x,y,u) -> -cosh(u)
A = [(x) -> 1;(x)->0;(x)->0;(x)->1]
bddata=[1 2 3 4;'D' 'N' 'D' 'N';(x,y)->exp(-x.^2-y.^2)  (x,y)->-2x.*exp(-x.^2-y.^2) (x,y)-> exp(-x.^2-y.^2) (x,y)->2x.*exp(-x.^2-y.^2)]

#u=solve_semlin_poisson(mesh,A,bddata,f,g)
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
	println("$(typeof(ar))")
	#stiffness matrix assembly
	ar4 = abs.(ar.*4)
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
	fa=transpose(vec(f.(vec(halfa_x),vec(halfa_y), 0)))
	fb=transpose(vec(f.(vec(halfb_x),vec(halfb_y), 0)))
	fc=transpose(vec(f.(vec(halfc_x),vec(halfc_y), 0)))
	rhs=[fa fb fc]

	Irhs=vec([a1 a2 a3])
	Jrhs=vec(ones(Int64,length(Irhs)))
	rhsint=rhs.*abs.([ar ar ar])./3 #triangle quadrature using medians
	println("$(typeof(rhsint))")
	println(rhsint)
	Krhs=vec(rhsint)
	println("$(typeof(Krhs))")

	###start boundary conditions

	### find out which BC applies where
	neu_bd=div.(findin(bddata,'N'),3)+1
	diri_bd=div.(findin(bddata,'D'),3)+1
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
   		neumanndata=sqrt.((mesh.nodes[1,neu_nodes1]-mesh.nodes[1,neu_nodes2]).^2+(mesh.nodes[2,neu_nodes1]-mesh.nodes[2,neu_nodes2]).^2).*neu_fct/2
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
		diri_edges_ind=findin(mesh.edges[4,:],bddata[1,diri_bd[j]])
		append!(diri_nodes,unique(vec(mesh.edges[1:2,diri_edges_ind])))
       	append!(diri_values,map(node->bddata[3,diri_bd[j]](mesh.nodes[1,node],mesh.nodes[2,node]),diri_nodes[end-length(diri_edges_ind):end]))
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
u = V

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





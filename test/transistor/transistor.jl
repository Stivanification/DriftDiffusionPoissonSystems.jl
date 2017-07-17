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
#        where \Sigma_j= sum over j, \sigma=0.3 influence parameter, cj=charge of jth atom, and          xj=position of the jth atom.
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
function Cdop_omega(x::Float64, y::Float64, omega::Array{Float64}, sigma::Float64, charge::Array{Float64})

	q = 1.6021766*10.0^(-19)
	sigmavec = vec([sigma sigma])
	n = size(omega,1)
	c = 0.0

	for i in 1:n
		dist = MvNormal(vec(omega[i,:]), sigmavec)
		c += q*charge[i]*pdf(dist, vec([x y]))
		end
	end
 c
end

function calculate_transistor(omega::Array{Float64})
		
	q = 1.6021766*10.0^(-19)
	sigma = 0.3
	Endpoints = [0.0 0.0; 60.0 0.0; 60.0 40.0; 0.0 40.0]
	Vbddata = [1 2 3 4; 'D' 'N' 'D' 'N'; (x,y) -> -1 (x,y) -> 0 (x,y) -> 0 (x,y) -> 0]
end


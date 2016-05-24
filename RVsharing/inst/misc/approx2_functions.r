# Utility functions for method 2 to approximate sharing probabilities in presence of unknown relationships

# By Alexandre Bureau

# Version 0.1
# 2013/05/16

infer.nalleles = function(phi,nf)
# Returns the most likely number of distinct alleles among nf founders based on mean estimated kinship phi
{
a = nf:(2*nf)
term1 = (2*nf-a)*(2*nf-a-1)/(nf*(nf-1))
term2 = 2*(a-nf)*2*(2*nf-a)/(2*(a-nf)*2*(2*nf-a) + (a-nf)*(2*(a-nf)-1) + 2*(2*nf-a)*(2*nf-a-1))
phi.vec = 0.5*term1 +  0.25*term2
phi.diff = (phi.vec - phi)
a[phi.diff>0&c(phi.diff[-1],0)<0]
#phi.vec
}

compute.phi.vec = function(nf,amin=2*nf-2)
# Compute the vector of expected phi_a for nf founders for numbers of distinct alleles a from amin to 2*nf-1
{
a = amin:(2*nf-1)
term1 = (2*nf-a)*(2*nf-a-1)/(nf*(nf-1))
term2 = 2*(a-nf)*2*(2*nf-a)/(2*(a-nf)*2*(2*nf-a) + (a-nf)*(2*(a-nf)-1) + 2*(2*nf-a)*(2*nf-a-1))
(0.5*term1 +  0.25*term2)/(nf-1)
}

infer.theta = function(phi,phi.vec)
# Solve the parameter theta when the distribution of a is limited to values from 2*nf-2 to 2*nf
# phi is the mean estimated kinship between founders
# phi.vec contains phi_a for a = 2*nf-2 and 2*nf-1
{
phi.diff = (phi - phi.vec)
(-phi.diff[2] - sqrt(phi.diff[2]^2 - 2*phi.diff[1]*phi))/phi.diff[1]
}
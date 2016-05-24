# Utility functions for method 2 to approximate sharing probabilities in presence of unknown relationships

# By Alexandre Bureau

# 2013/06/05

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
# Sequential probability computation
{
a = amin:(2*nf-1)
term1 = (2*nf-a)*(2*nf-a-1)/(nf*(nf-1))
term2 = (a-nf)*(2*nf-a)/(nf*(nf-1)) + 2*(a-nf)*(2*nf-a)/(nf*(2*nf-1))
(0.5*term1 +  0.25*term2)/(nf-1)
}

infer.theta = function(phi,phi.vec)
# Solve the parameter theta for polynomial approximation of the distribution of the number of distinct alleles
# This is a general function for polynomials of order 2 to 5. The previous version of this function was doing
# the computation for a quadratic polynomial only

# Arguments:
# phi is the mean estimated kinship between founders
# phi.vec contains phi_a for a = 2*nf-ord to 2*nf-1, where ord must be between 2 and 5

# Value:
# Real roots of the polynomial approximation
{
ord = length(phi.vec)
phi.diff = (phi - phi.vec)
coef.vec = c(1,1/2,1/6,1/24,1/120)[1:ord]
racines = polyroot(c(phi,phi.diff[ord:1]*coef.vec))
# Return only the real roots
Re(racines)[abs(Im(racines))<1e-10]
}

# Special cases of infer.theta which are no longer needed
infer.theta.quadratic = function(phi,phi.vec)
# Solve the parameter theta when the distribution of a is limited to values from 2*nf-4 to 2*nf
# phi is the mean estimated kinship between founders
# phi.vec contains phi_a for a = 2*nf-2 and 2*nf-1
{
phi.diff = (phi - phi.vec)
(-phi.diff[2] - sqrt(phi.diff[2]^2 - 2*phi.diff[1]*phi))/phi.diff[1]
}

infer.theta.cubic = function(phi,phi.vec)
# Solve the parameter theta when the distribution of a is limited to values from 2*nf-3 to 2*nf
# phi is the mean estimated kinship between founders
# phi.vec contains phi_a for a = 2*nf-3 to 2*nf-1
{
phi.diff = (phi - phi.vec)
racines = polyroot(c(phi,phi.diff[3:1]*c(1,1/2,1/6)))
# Return only the real roots
Re(racines)[abs(Im(racines))<1e-10]
}

infer.theta.order4 = function(phi,phi.vec)
# Solve the parameter theta when the distribution of a is limited to values from 2*nf-3 to 2*nf
# phi is the mean estimated kinship between founders
# phi.vec contains phi_a for a = 2*nf-4 to 2*nf-1
{
phi.diff = (phi - phi.vec)
racines = polyroot(c(phi,phi.diff[4:1]*c(1,1/2,1/6,1/24)))
# Return only the real roots
Re(racines)[abs(Im(racines))<1e-10]
}

infer.theta.order5 = function(phi,phi.vec)
# Solve the parameter theta when the distribution of a is limited to values from 2*nf-3 to 2*nf
# phi is the mean estimated kinship between founders
# phi.vec contains phi_a for a = 2*nf-5 to 2*nf-1
{
phi.diff = (phi - phi.vec)
racines = polyroot(c(phi,phi.diff[5:1]*c(1,1/2,1/6,1/24,1/120)))
# Return only the real roots
Re(racines)[abs(Im(racines))<1e-10]
}


get.LODallshare <- function(vec,pshare)
{
if (any(pshare$ped.tocompute.vec%in%vec)) sum(pshare$mlog10pshare[pshare$ped.tocompute.vec%in%vec])
else NA
}

# Wrappers for pedigree object
# Returns only pshare
RVsharing.ped.pshare = function(ped)
{
id = ped$id
dad.id = mom.id = numeric(length(id))
dad.id[ped$findex>0] = ped$id[ped$findex]
mom.id[ped$mindex>0] = ped$id[ped$mindex]
RVsharing(id,dad.id,mom.id)$pshare
} 

# Returns object
RVsharing.ped = function(ped)
{
id = ped$id
dad.id = mom.id = numeric(length(id))
dad.id[ped$findex>0] = ped$id[ped$findex]
mom.id[ped$mindex>0] = ped$id[ped$mindex]
RVsharing(id,dad.id,mom.id)
}

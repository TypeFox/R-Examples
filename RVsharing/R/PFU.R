# Computation of P[FjU]

PFU.direct = function(nf,theta,ord=2)
{
# Computation of P[FjU] using equation 21 of Bureau et al.

# Arguments:
# nf: number of founders of the pedigree
# theta : value of the parameter of the polynomial distribution
# ord : order of the polynomial approximation to the distribution of the number 
#       of distinct alleles in the founders (noted d in Bureau et al.) Must be <= 5.

# Value: P[FjU] (scalar)

a = (2*nf):(2*nf-ord)
distri = c(1,theta,theta^2/2,theta^3/6,theta^4/24,theta^5/120)[1:(ord+1)]
weighted.mean(2/nf - 2/a,distri)
}

# Special cases of PFU.direct which are no longer needed
PFU.direct.cubic = function(nf,theta)
{
a = (2*nf-3):(2*nf)
weighted.mean(2/nf - 2/a,c(theta^3/6,theta^2/2,theta,1))
}

PFU.direct.order4 = function(nf,theta)
{
a = (2*nf-4):(2*nf)
weighted.mean(2/nf - 2/a,c(theta^4/24,theta^3/6,theta^2/2,theta,1))
}

## PFU.fromphi = function(phi.vec,theta)
## {
## # 2*nf - 2 alleles introduced
## PF2 = (2 - (Nf-1)*2*phi.vec[1])/(2*Nf-2)
## # 2*nf - 1 alleles introduced
## PF1 = (2 - (Nf-1)*2*phi.vec[2])/(2*Nf-1)

## PFU = (PF2 * theta^2/2 + PF1 * theta + 1/Nf)/(1 + theta + 1/2*theta^2)
## }


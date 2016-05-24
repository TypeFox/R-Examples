# Computation of P[FjU]
# 2*nf - 2 alleles introduced

PFU.fromphi = function(phi.vec,theta)
{
PF2 = (2 - (Nf-1)*2*phi.vec[1])/(2*Nf-2)
# 2*nf - 1 alleles introduced
PF1 = (2 - (Nf-1)*2*phi.vec[2])/(2*Nf-1)

PFU = (PF2 * theta^2/2 + PF1 * theta + 1/Nf)/(1 + theta + 1/2*theta^2)
}

# The function below needs testing
PFU.direct(nf,theta)
{
a = (2*nf-2):(2*nf)
weighted.mean((2/a - 42*(2*nf/a - 1)/nf),c(theta^2/2,theta,1))
}
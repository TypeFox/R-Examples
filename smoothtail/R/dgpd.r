dgpd <- function(x, gam, sigma = 1){
#
# Density function of the generalized Pareto
# distribution (GPD).
#
# Input:
# - x     : places where the density function shall be evaluated. May
#           be a vector. Note that x E (0,-1/gam) if gam<0, x E R if gam>0
# - gam   : tail index
# - sigma : scale parameter
#
# Kaspar Rufibach, 2006
#
n <- length(x)
d <- 1:n*NA
if (gam<0){
    II <- (1:n)*(x>=0)*(x < -sigma / gam)
    d[II] <- (1+gam*x[II]/sigma)^(-(1+1/gam))}

II <- (1:n)*(x>=0)

if (gam==0){d[II] <- exp(-x[II]/sigma)}

if (gam>0){d[II] <- (1+gam*x[II]/sigma)^(-(1+1/gam))}

return(d / sigma)
}

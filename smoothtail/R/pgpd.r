pgpd <- function(q, gam, sigma = 1){
#
# Distribution function of the generalized Pareto
# distribution (GPD).
#
# Input:
# - q     : places where the distribution function shall be evaluated. May
#           be a vector. Note that q E (0,-1/gam) if gam<0, q E R if gam>0
# - gam   : tail index
# - sigma : scale parameter
#
# Kaspar Rufibach, 2006
#
n <- length(q)
p <- 1:n*NA
if (gam<0){
    II <- (1:n)*(q>=0)*(q< -sigma/gam)
    p[II] <- 1-(1+gam*q[II]/sigma)^(-1/gam)}

II <- (1:n)*(q>=0)    

if (gam==0){p[II] <- 1-exp(-q[II]/sigma)}
    
if (gam>0){p[II] <- 1-(1+gam*q[II]/sigma)^(-1/gam)}

return(p)}


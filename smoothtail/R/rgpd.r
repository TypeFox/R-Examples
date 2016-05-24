rgpd <- function(n, gam, sigma = 1){
#
# Random number generation for the generalized Pareto
# distribution (GPD).
#
# Input:
# - n     : Number of random numbers to be generated.
# - gam   : tail index
# - sigma : scale parameter
#
# Kaspar Rufibach, 2006
#
u <- runif(n,0,1)
r <- qgpd(u, gam, sigma)

err <- 0
if (n<=0){
    err <- 1
    cat('Please give a strict positive n!')}

if (err==0){return(sort(r))}
}


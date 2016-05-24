qgpd <- function(p, gam, sigma = 1){
#
# Quantile function of the generalized Pareto
# distribution (GPD).
#
# Input:
# - p     : Vector of probabilities.
#
# - gam   : tail index
# - sigma : scale parameter
#
# Kaspar Rufibach, 2006
#
n <- length(p)
q <- 1:n*NA

for (i in 1:n){
    if (p[i] <= 0){q[i] <- -Inf}
    if (p[i] == 1){q[i] <- max(p)}
    if (p[i] >  1){q[i] <- Inf}
    
    if (is.na(q[i]) == 1){
        if (gam!=0){q[i] <- ((1-p[i])^(-gam)-1)/gam}
        if (gam==0){q[i] <- -log(1-p[i])}}
}
return(q * sigma)
}


# Calculate estimators for Hs and Ht
#
# This function calculates Nei and Chesser's estimators for Hs and Ht
# 
# Note, the individuals functions in the diff_stat family (listed below)
# return these estimates (becuase they use this function internally). This
# function is not exported (if someone has a use for it I can change this).

HsHt <- function(x){
    #start by dropping missing for this locus
    x <- x[complete.cases(x@tab)]
    n_by_pop <-  table(pop(x))
    if(length(n_by_pop) < 2){
        warning("Need at least two population to calculate differentiation")
        return(c(NA, NA, NA))
    }
    afreqs <- apply(x@tab, 2, "/", ploidy(x))
    n <- nPop(x)
    harmN <- harmonic_mean(n_by_pop[n_by_pop >0])
    a <- apply(afreqs, 2, function(A) tapply(A, pop(x), mean))
    HpS <- mean(1 - rowSums(a^2))
    Hs_est <- (2*harmN/(2*harmN-1))*HpS
    HpT <- 1 - sum(colMeans(a)^2)
    Ht_est <- HpT + Hs_est/(2*harmN*n)
    return(c(Ht_est, Hs_est, n))
}

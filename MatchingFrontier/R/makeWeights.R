# See j.mp/CEMweights for an explanation of weights

makeWeights <-
function(dataset, treatment){
    strata <- getMahalStrata(rownames(dataset), dataset$matched.to)
    num.treated <- sum(dataset[[treatment]] == 1)
    num.control <- sum(dataset[[treatment]] == 0)
    w <- rep(NA, nrow(dataset))
    w[dataset[[treatment]] == 1] <- 1
    for(s in unique(strata)){
        num.treated.strata <- sum(strata == s & dataset[[treatment]] == 1)
        num.control.strata <- sum(strata == s & dataset[[treatment]] == 0)
        w[strata == s & dataset[[treatment]] == 0] <- (num.control / num.treated) * (num.treated.strata / num.control.strata)
    }
    return(w)
}

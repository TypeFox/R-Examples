compare <- function(...){
    ## function to allow quick comparisons between various designs
    liste <- as.list(match.call())[-1]
    if (length(liste)<2) stop("at least two designs are needed")
    check.designs <- sapply(liste, function(obj) design.info(eval(obj))$type)
    if (!(all(check.designs=="lhs") | all(check.designs=="Dopt")))
       stop("invalid combination of designs")
    di <- design.info(eval(liste[[1]]))
    aus <- cbind(c(nruns=di$nruns, nfactors=di$nfactors, unlist(di$optimality.criteria)))
    for (i in 2:length(liste)){
            di <- design.info(eval(liste[[i]]))
            aus <- cbind(aus, c(nruns=di$nruns, nfactors=di$nfactors, unlist(di$optimality.criteria)))
            }
    colnames(aus) <- as.character(unlist(liste))
    aus
}
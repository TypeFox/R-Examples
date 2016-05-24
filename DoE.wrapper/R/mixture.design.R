mixture.design <- function(nfactors=NULL, nlevels, randomize=TRUE, seed=NULL, replications=1,
     repeat.only=FALSE, factor.names=if (!is.null(nfactors)) paste("X",1:nfactors,sep="") else NULL, ...){
    #    aus <- gen.mixture(nlevels, factor.names)
        ## incorporate randomization and replication
        ## could be run with large nlevels
        ## and subsequent conditioning
        ## as entry point to debarred mixtures
    cat("This does not work yet.\n")
}
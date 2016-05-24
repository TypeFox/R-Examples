#Procrustes rotation for each bootstrap step, uses procrustes() function from MCMCpack
qpcrustes <- function(loa, target, nfactors) {
    prox <- as.matrix(loa)
    #prores <- procrustes(prox, target)
    # loarot <- as.data.frame(prores[1])
    warning("The procrustes rotation is not working currently due to an issue in package dependency from 'MCMCpack' and 'graph'. If you'd like to try this: 1) see the code for the function 'qpcrustes' and 2) uncomment lines 4 and 5 and comment line 7.")
    loarot <- as.data.frame(prox)
    colnames(loarot) <- paste("loarot_f", 1:nfactors, sep="")
    return(loarot)
}
print.population = function(x,...){
    cat("Substructured population\n")
    cat("------------------------\n")
    cat(paste("Number of members:", x$nProfiles, "\n"))
    cat(paste("Number of subpopulations", x$nSubpops, "\n"))
    cat(paste("Number of loci:", x$nLoci, "\n"))
    cat(paste("Expected level of coancestry (theta):", x$theta, "\n"))
}

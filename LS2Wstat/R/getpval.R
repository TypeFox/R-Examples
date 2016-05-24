getpval <- function (statvec, verbose = TRUE){

statlen <- length(statvec)

if (verbose){
        cat("Observed bootstrap is ", round(statvec[1],3), "\n")
}

p <- (1+sum(statvec[1] <= statvec[2:statlen]))/statlen

if (verbose) {
        cat("p-value is ", p, "\n")
}

return(p)
}

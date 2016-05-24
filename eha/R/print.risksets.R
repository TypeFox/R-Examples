print.risksets <- function(x, ...){
    if (class(x) != "risksets") stop("Only for class 'risksets'")
    cat("No of strata: ", x$ns, "\n")
    if (x$ns > 1){
        for (i in 1:x$ns){
            cat("Stratum No. ", i, "\n")
            cat("--------------\n")
            cat("Number of risksets: ", x$antrs[i], "\n")
            cat("Number of events per risk set: ", summary(x$n.events), "\n")
        }
    }
}

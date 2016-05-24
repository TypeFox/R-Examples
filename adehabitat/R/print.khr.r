"print.khr" <- function(x, ...)
{
    if (!inherits(x, "khr"))
        stop("x should be an object of class khr")
    cat("********** Utilization distribution of Animals ************\n\n")
    if (inherits(x, "khrud"))
        cat("Type: probability density\n")
    if (inherits(x, "kbbhrud"))
        cat("Type: probability density estimated with the Brownian bridge approach\n")
    if (inherits(x, "khrvol"))
        cat("Type: volume under UD (only used to compute home ranges)\n")
    cat("\nUD have been estimated using the kernel method for the following animals:\n")

    print(names(x), quote=FALSE)
    th<-x[[1]]$hmeth
    if (th=="LSCV")
        cat("\nThe smoothing parameter was estimated by cross validation\n")
    if (th=="href")
        cat("\nThe smoothing parameter was estimated by the reference method (ad hoc)\n")
    if (is.numeric(th))
        cat("\nThe smoothing parameter was set to", th, "\n")

    cat("\nEach animal is a component of the list, and for each animal,\n")
    cat("the following elements are available:\n")
    cat("$UD       The utilization distribution (object of class \"asc\")\n")
    cat("$locs     The relocations of the animal\n")
    if (th=="LSCV") {
        cat("$h        A list with the following components:\n")
        cat("          $CV   The results of cross-validation\n")
        cat("          $h    The value of the smoothing parameter\n")
    }
    if (th=="href") {
        cat("$h        The value of the smoothing parameter\n")
    }
    if (th=="bb") {
        cat("$h        The values of the smoothing parameters\n")
    }

    if (th=="LSCV") {
        m<-TRUE
        for (i in 1:length(x))
            m[i]<-x[[i]]$h$convergence
        names(m)<-names(x)
        if (!all(m)) {
            cat("\nWARNING!! No convergence in cross-validation for the following animals:\n")
            print(names(m)[!m], quote=FALSE)
            cat("Consider a new fit of UD using the ad hoc method for h.\n")
        }
    }
}


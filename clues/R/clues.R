# v0.3.1 created on Feb. 15, 2009 by Weiliang Qiu
#  (1) if number of clusters = number of data points,
#      then we set number of clusters = 1, 
#      and  index (or silhouette index) = NULL
#
# v0.2.4 created on Feb. 4, 2009 by Weiliang Qiu
#  (1) added code to output final results
#  
# y - data matrix. Rows are observations; Columns are variables
# n0 - initial guess of the number of clusters
clues <- function(y, 
    n0 = 5, 
    alpha = 0.05, 
    eps = 1.0e-4, 
    itmax = 20, 
    K2.vec = NULL, 
    strengthMethod = "sil", 
    strengthIni = -3, 
    disMethod = "Euclidean",
    quiet = TRUE)
{
    if(is.null(K2.vec))
    {
      K2.vec=n0
    }

    strengthMethod <- match.arg(strengthMethod, c("sil", "CH"))
    disMethod <- match.arg(disMethod, c("Euclidean", "1-corr"))
 
    # y-- data matrix; rows are observations; columns are variables
    if(!is.matrix(y))
    { y <- matrix(y, ncol = 1) }
 
    nr <- nrow(y)
 
    if(strengthMethod == "CH")
    {
        res <- clues_CH(y=y, n0=n0, alpha=alpha, eps=eps, 
          itmax=itmax, K2.vec=K2.vec,  CH2=strengthIni,
            disMethod=disMethod, quiet=quiet)
    } else if (strengthMethod == "sil") {
        res <- clues_sil(y=y, n0=n0, alpha=alpha, eps=eps, itmax=itmax, 
             K2.vec=K2.vec, s2=strengthIni, 
             disMethod=disMethod, quiet=quiet)
    } else {
        stop("'strengthMethod' does not match 'CH' or 'sil'!\n")
    }
 
    if(nr == res$g)
    { res$g <- 1
      res$mem <- rep(1, nr)
      res$size <- nr
      if(strengthMethod == "CH")
      { 
          res$CH <- NULL
      } else {
          res$avg.s <- NULL
          res$s <- rep(NULL, nr)
      }
    }
 
    if(!quiet)
    {
        cat("********* output results ********\n")
        if(strengthMethod == "CH")
        {
            cat("CH = ", res$CH, "\n")
        } else {
            cat("avg Silhouette = ", res$avg.s, "\n")
        }
        cat("Estimated number of clusters = ", res$g, "\n")
        cat("*********************************\n")
    }
 
    res$y <- y
    res$strengthMethod <- strengthMethod
    res$disMethod <- disMethod

    class(res) <- c("clues", "partition")
    invisible(res)
}

###########
# written based on the functions in the file 'pam.q' in the package 'cluster'
###########

## non-exported:
.print.clues <- function(x, ...) {
    cat("Number of data points:\n")
    print(nrow(x$y.old1))
    cat("\n")
    cat("Number of variables:\n")
    print(ncol(x$y.old1))
    cat("\n")

    cat("Number of clusters:\n")
    print(x$g)
    cat("\n")
    cat("Cluster sizes:\n")
    print(x$size, ...)
    cat("\n")
    cat("Strength method:\n")
    if(x$strengthMethod == "CH")
    {
        print("CH")
        cat("\n")
        cat("CH:\n")
        print(x$CH)
        cat("\n")
    } else {
        print("sil")
        cat("\n")
        cat("avg Silhouette:\n")
        print(x$avg.s)
        cat("\n")
    }

    cat("dissimilarity measurement:\n")
    print(x$disMethod, ...)
    cat("\n")

}

print.clues <- function(x, ...)
{
    .print.clues(x, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.clues <- function(object, ...)
{
    class(object) <- "summary.clues"
    object
}

print.summary.clues <- function(x, ...)
{
    .print.clues(x, ...)

    cat("\nAvailable components:\n")
    print(names(x), ...)
    cat("\n")

    invisible(x)
}


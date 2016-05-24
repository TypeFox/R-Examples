`summary.mvloc` <-
function(object,..., digits=4)
    {   
    cat("The", object$est.name, "of", object$dname, "is:\n") 
    print(format(round(object$location,digits)),quote=F)
    cat("\n")
    cat("And has the covariance matrix:")
    cat("\n")
    print(round(object$vcov,digits))
    invisible(object)
    }

#############################################################
#
#	print.localdepth function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: December, 17, 2008
#	Version: 0.1-2
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

print.localdepth <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Local Depth Statistics: \n")
    print(summary(x$localdepth), digits=digits, ...)
    cat("Depth Statistics: \n")
    print(summary(x$depth), digits=digits, ...)
    cat("\n")
    cat("Information \n")
    cat("Method: ", x$method, "\n")
    cat("Threshold parameter 'tau' is set to: ", x$tau, "\n")
    cat("The dimension is measured by: ", x$use, "\n")
    cat("Membership evaluation: ", x$type, "\n")
    if (x$nsamp=='all')
      cat("All objects were explored \n")
    else {
      if (x$method!='mahalanobis') {
        cat("A Monte Carlo approach was used on: ", x$num[1], " objects in order to have ", x$num[2])
        if (x$method=='simplicial')
          cat(" simplicies ")
        else if (x$method=='ellipsoid')
          cat(" ellipsoids ")
        cat("with size smaller than 'tau'\n")
      } else {
        cat("A Monte Carlo approach was used on: ", x$num[1], " objects\n")
      }
    }
    cat("\n")
    invisible(x)
}

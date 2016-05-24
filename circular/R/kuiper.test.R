
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   kuiper.test function                                    #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 27, 2006                                     #
#   Version: 0.3                                            #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

kuiper.test <- function(x, alpha=0) {

    # Handling missing values  
    x <- na.omit(x)
    if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "circularp") <- attr(x, "class") <- NULL
    result <- list()
    result$call <- match.call()
    result$statistic <- KuiperTestRad(x, alpha)
    result$alpha <- alpha
    class(result) <- "kuiper.test"
    return(result)
}

KuiperTestRad <- function(x, alpha) {
   if (!any(c(0, 0.01, 0.025, 0.05, 0.1, 0.15)==alpha)) stop("'alpha' must be one of the following values: 0, 0.01, 0.025, 0.05, 0.1, 0.15")
   x <- sort(x %% (2 * pi))/(2 * pi)
   n <- length(x)
   i <- 1:n
   D.P <- max(i/n - x)
   D.M <- max(x - (i - 1)/n)
   V <- (D.P + D.M) * (sqrt(n) + 0.155 + 0.24/sqrt(n))
   return(V)
}

#############################################################
#                                                           #
#   print.kuiper.test function                              #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 19, 2003                                #
#   Version: 0.1-1                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.kuiper.test <- function(x, digits=4, ...) {
    V <- x$statistic
    alpha <- x$alpha
    kuiper.crits <- cbind(c(0.15, 0.1, 0.05, 0.025, 0.01), c(1.537, 1.62, 1.747, 1.862, 2.001))
    cat("\n", "      Kuiper's Test of Uniformity", "\n", "\n")
    cat("Test Statistic: ", round(V, digits=digits), "\n")
    if (alpha == 0) {
        if (V < 1.537)
        cat("P-value > 0.15", "\n", "\n")
    else if (V < 1.62)
         cat("0.10 < P-value < 0.15", "\n", "\n")
    else if (V < 1.747)
         cat("0.05 < P-value < 0.10", "\n", "\n")
    else if (V < 1.862)
         cat("0.025 < P-value < 0.05", "\n", "\n")
    else if (V < 2.001)
         cat("0.01 < P-value < 0.025", "\n", "\n")
    else cat("P-value < 0.01", "\n", "\n")
    } else {
    Critical <- kuiper.crits[(1:5)[alpha == c(kuiper.crits[, 1])],2]
    cat("Level", alpha, "Critical Value:", round(Critical, 4), "\n")
    if (V > Critical)
        cat("Reject Null Hypothesis", "\n", "\n")
    else cat("Do Not Reject Null Hypothesis", "\n", "\n")
    }
    invisible(x)
}


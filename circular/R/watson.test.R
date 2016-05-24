
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   watson.test function                                    #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 31, 2006                                     #
#   Version: 0.3-1                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

watson.test <- function(x, alpha = 0, dist = c("uniform", "vonmises")) {
    # Handling missing values
    x <- na.omit(x)
    if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }
    dist <- match.arg(dist)
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "circularp") <- attr(x, "class") <- NULL
    if (!any(c(0, 0.01, 0.025, 0.05, 0.1)==alpha))
       stop("'alpha' must be one of the following values: 0, 0.01, 0.025, 0.05, 0.10")
    result <- WatsonTestRad(x, dist)
    result$call <- match.call()
    result$n <- length(x)
    result$alpha <- alpha
    result$dist <- dist
    class(result) <-"watson.test"
    return(result)
}

WatsonTestRad <- function(x, dist) {
    n <- length(x)
    if (dist == "uniform") {
    u <- sort(x)/(2 * pi)
    u.bar <- mean.default(u)
    i <- 1:n
    sum.terms <- (u - u.bar - (2 * i - 1)/(2 * n) + 0.5)^2
    u2 <- sum(sum.terms) + 1/(12 * n)
    u2 <- (u2 - 0.1/n + 0.1/(n^2)) * (1 + 0.8/n)
        result <- list(statistic=u2, row=NA)
    } else {
        res <- MlevonmisesRad(x, bias=FALSE)
    mu.hat <- res[1]
    kappa.hat <- res[4]
    x <- (x - mu.hat) %% (2 * pi)
    x <- matrix(x, ncol = 1)
    z <- apply(x, 1, PvonmisesRad, mu=0, kappa=kappa.hat, tol=1e-020)
    z <- sort(z)
    z.bar <- mean.default(z)
    i <- 1:n
    sum.terms <- (z - (2 * i - 1)/(2 * n))^2
    Value <- sum(sum.terms) - n * (z.bar - 0.5)^2 + 1/(12 * n)                
    if (kappa.hat < 0.25)
        row <- 1
        else if (kappa.hat < 0.75)
            row <- 2
        else if (kappa.hat < 1.25)
            row <- 3
        else if (kappa.hat < 1.75)
            row <- 4
        else if (kappa.hat < 3)
            row <- 5
        else if (kappa.hat < 5)
            row <- 6
        else row <- 7   
        result <- list(statistic=Value, row=row)
    }
    return(result)
}

#############################################################
#                                                           #
#   print.watson.test function                              #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 19, 2003                                #
#   Version: 0.1-1                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.watson.test <- function(x, digits=4, ...) {
    dist <- x$dist
    n <- x$n
    alpha <- x$alpha

    if (dist == "uniform") {
        u2 <- x$statistic
    cat("\n", "      Watson's Test for Circular Uniformity", "\n", "\n")
    crits <- c(99, 0.267, 0.221, 0.187, 0.152)
    if (n < 8) {
        warning("Total Sample Size < 8:  Results may not be valid", "\n", "\n")
    }
        cat("Test Statistic:", round(u2, digits=digits), "\n")
    if (alpha == 0) {
        if (u2 > 0.267)
        cat("P-value < 0.01", "\n", "\n")
        else if (u2 > 0.221)
             cat("0.01 < P-value < 0.025", "\n", "\n")
        else if (u2 > 0.187)
             cat("0.025 < P-value < 0.05", "\n", "\n")
        else if (u2 > 0.152)
             cat("0.05 < P-value < 0.10", "\n", "\n")
        else cat("P-value > 0.10", "\n", "\n")
    } else {
        index <- (1:5)[alpha == c(0, 0.01, 0.025, 0.05, 0.1)]
        Critical <- crits[index]
        if (u2 > Critical)
        Reject <- "Reject Null Hypothesis"
        else Reject <- "Do Not Reject Null Hypothesis"
        cat("Level", alpha, "Critical Value:", round(Critical, digits=digits), "\n")
        cat(Reject, "\n\n")
        }
    
    } else if (dist=="vonmises") {
               Value <- x$statistic
               row <- x$row
           cat("\n", "      Watson's Test for the von Mises Distribution \n\n")
               u2.crits <- cbind(c(0, 0.5, 1, 1.5, 2, 4, 100), c(0.052, 0.056, 0.066, 0.077, 0.084, 0.093, 0.096), c(0.061, 0.066, 0.079, 0.092, 0.101, 0.113, 0.117), c(0.081, 0.09, 0.11, 0.128, 0.142, 0.158, 0.164))
 
           if (alpha != 0) {
           if (alpha == 0.1)
               col <- 2
           else if (alpha == 0.05)
               col <- 3
           else if (alpha == 0.01)
            col <- 4
           Critical <- u2.crits[row, col]
           if (Value > Critical)
               Reject <- "Reject Null Hypothesis"
           else Reject <- "Do Not Reject Null Hypothesis"
           cat("Test Statistic:", round(Value, digits=digits), "\n")
           cat("Level", alpha, "Critical Value:", round(Critical, digits=digits), "\n")
           cat(Reject, "\n\n")
           } else {
           cat("Test Statistic:", round(Value, digits=digits), "\n")
           if (Value < u2.crits[row, 2])
               cat("P-value > 0.10", "\n", "\n")
           else if ((Value >= u2.crits[row, 2]) && (Value < u2.crits[row, 3]))
                cat("0.05 < P-value > 0.10", "\n", "\n")
           else if ((Value >= u2.crits[row, 3]) && (Value < u2.crits[row, 4]))
                cat("0.01 < P-value > 0.05", "\n", "\n")
           else cat("P-value < 0.01", "\n", "\n")
           }

           }
    invisible(x)
}


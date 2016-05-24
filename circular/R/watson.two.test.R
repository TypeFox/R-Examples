
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   watson.two.test function                                #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: December, 16, 2009                                #
#   Version: 0.3-2                                          #
#                                                           #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#############################################################

watson.two.test <- function(x, y, alpha=0) {

    # Handling missing values
    x <- na.omit(x)
    if (length(x)==0) {
        warning("'x': No observations (at least after removing missing values)")
        return(NULL)
    }      
    y <- na.omit(y)
    if (length(y)==0) {
        warning("'y': No observations (at least after removing missing values)")
        return(NULL)
    }      

    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "circularp") <- attr(x, "class") <- NULL
    y <- conversion.circular(y, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(y, "circularp") <- attr(y, "class") <- NULL

    result <- WatsonTwoTestRad(x, y)
      
    result$call <- match.call()
    result$alpha <- alpha
    class(result) <- "watson.two.test"
    return(result)
}

WatsonTwoTestRad <- function(x, y) {
    n1 <- length(x)
    n2 <- length(y)
    n <- n1 + n2
    x <- cbind(sort(x %% (2 * pi)), rep(1, n1))
    y <- cbind(sort(y %% (2 * pi)), rep(2, n2))
    xx <- rbind(x, y)
    rank <- order(xx[, 1])
    xx <- cbind(xx[rank,  ], 1:n)
    a <- 1:n
    b <- 1:n
    for (i in 1:n) {
       a[i] <- sum(xx[1:i, 2] == 1)
       b[i] <- sum(xx[1:i, 2] == 2)
    }
    d <- b/n2 - a/n1
    dbar <- mean.default(d)
    u2 <- (n1 * n2)/n^2 * sum((d - dbar)^2)
    result <- list(statistic=u2, nx=n1, ny=n2)
    return(result)
}

#############################################################
#                                                           #
#   print.watson.two.test function                          #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: November, 19, 2003                                #
#   Version: 0.1-1                                          #
#                                                           #
#   Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.watson.two.test <- function(x, digits=4, ...) {
    u2 <- x$statistic
    n1 <- x$nx
    n2 <- x$ny
    alpha <- x$alpha
    n <- n1 + n2
    cat("\n      Watson's Two-Sample Test of Homogeneity \n\n")
    if (n < 18)
    warning("Total Sample Size < 18:  Consult tabulated critical values \n\n")
    crits <- c(99, 0.385, 0.268, 0.187, 0.152)
    cat("Test Statistic:", round(u2, digits=digits), "\n")
    
    if (alpha == 0) {
    if (u2 > 0.385)
        cat("P-value < 0.001", "\n", "\n")
    else if (u2 > 0.268)
         cat("0.001 < P-value < 0.01", "\n", "\n")
    else if (u2 > 0.187)
         cat("0.01 < P-value < 0.05", "\n", "\n")
    else if (u2 > 0.152)
         cat("0.05 < P-value < 0.10", "\n", "\n")
    else cat("P-value > 0.10", "\n", "\n")
    } else {
    index <- (1:5)[alpha == c(0, 0.001, 0.01, 0.05, 0.1)]
    Critical <- crits[index]
    if (u2 > Critical)
        Reject <- "Reject Null Hypothesis"
    else Reject <- "Do Not Reject Null Hypothesis"
    cat("Level", alpha, "Critical Value:", round(Critical, digits=digits), "\n")
    cat(Reject, "\n\n")
    }
invisible(x)
}

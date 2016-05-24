
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   rao.spacing.test function                               #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 31, 2006                                     #
#   Version: 0.3-1                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

rao.spacing.test <- function(x, alpha = 0) {

    # Handling missing values
    x <- na.omit(x)
    if ((n <- length(x))==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }  
    if (!any(c(0, 0.01, 0.025, 0.05, 0.1, 0.15)==alpha))
       stop("'alpha' must be one of the following values: 0, 0.01, 0.025, 0.05, 0.1, 0.15")

    x <- conversion.circular(x, units="degrees", zero=0, rotation="counter", modulo="2pi")
    attr(x, "circularp") <- attr(x, "class") <- NULL

    statistic <- RaoSpacingTestDeg(x)
    result <- list()
    result$call <- match.call()
    result$statistic <- statistic
    result$alpha <- alpha
    result$n <- n
    class(result) <- "rao.spacing.test"
    return(result)
}

RaoSpacingTestDeg <- function(x) {
    x <- sort(x %% 360)
    n <- length(x)
    if (n < 4) {
       warning("Sample size too small")
       U <- NA
    } else {
       spacings <- c(diff(x), x[1] - x[n] + 360)
       U <- 1/2 * sum(abs(spacings - 360/n))
    }
    return(U)
}

#############################################################
#                                                           #
#   print.rao.spacing.test function                         #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: August, 15, 2007                                  #
#   Version: 0.2-2                                          #
#                                                           #
#   Copyright (C) 2007 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.rao.spacing.test <- function(x, digits=4, ...) {
    U <- x$statistic
    alpha <- x$alpha
    n <- x$n
    data(rao.table, package='circular', envir=sys.frame(which=sys.nframe()))
    if (n <= 30)
    table.row <- n - 3
    else if (n <= 32)
         table.row <- 27
    else if (n <= 37)
         table.row <- 28
    else if (n <= 42)
         table.row <- 29
    else if (n <= 47)
         table.row <- 30
    else if (n <= 62)
         table.row <- 31
    else if (n <= 87)
         table.row <- 32
    else if (n <= 125)
         table.row <- 33
    else if (n <= 175)
         table.row <- 34
    else if (n <= 250)
         table.row <- 35
    else if (n <= 350)
         table.row <- 36
    else if (n <= 450)
         table.row <- 37
    else if (n <= 550)
         table.row <- 38
    else if (n <= 650)
         table.row <- 39
    else if (n <= 750)
         table.row <- 40
    else if (n <= 850)
         table.row <- 41
    else if (n <= 950)
         table.row <- 42
    else table.row <- 43
        
    cat("\n")
    cat("       Rao's Spacing Test of Uniformity", "\n", "\n")
    cat("Test Statistic =", round(U, digits=digits), "\n")
    
    if (alpha == 0) {
        if (U > rao.table[table.row, 1])
        cat("P-value < 0.001", "\n", "\n")
        else if (U > rao.table[table.row, 2])
             cat("0.001 < P-value < 0.01", "\n", "\n")
        else if (U > rao.table[table.row, 3])
             cat("0.01 < P-value < 0.05", "\n", "\n")
        else if (U > rao.table[table.row, 4])
             cat("0.05 < P-value < 0.10", "\n", "\n")
        else cat("P-value > 0.10", "\n", "\n")
        x$accepted <- NA
    } else {
        table.col <- (1:4)[alpha == c(0.001, 0.01, 0.05, 0.1)]
        critical <- rao.table[table.row, table.col]
        cat("Level", alpha, "critical value =", critical, "\n")
        if (U > critical) {
           cat("Reject null hypothesis of uniformity \n\n")
           x$accepted <- FALSE
        } else {
           cat("Do not reject null hypothesis of uniformity \n\n")
           x$accepted <- TRUE
        }
  }
invisible(x)
}

#############################################################
#                                                           #
#	wle.var.test function                               #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April 3, 2003                                 #
#	Version: 0.2                                        #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#	Based on var.test function in                       #
#       ctest package version 1.2.0                         #
#                                                           #
#############################################################

wle.var.test <- function(x, y, ratio = 1, alternative = c("two.sided", "less", "greater"), conf.level = 0.95, x.root=1, y.root=1) {

    if (!((length(ratio) == 1) && is.finite(ratio) && (ratio > 0)))
        stop("ratio must be a single positive number")

    alternative <- match.arg(alternative)

    if (!((length(conf.level) == 1) && is.finite(conf.level) &&
          (conf.level > 0) && (conf.level < 1)))
        stop("conf.level must be a single number between 0 and 1")

    DNAME <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))

    if (inherits(x, "wle.lm") && inherits(y, "wle.lm")) {

        x.tot.sol <- x$tot.sol
        if (x.tot.sol<x.root) {
            stop(paste("'x' Root ",x.root," not found"))
        }
        if (x.tot.sol!=1) {
            x.res <- x$residuals[x.root,]
            x.c <- x$coefficients[x.root,]
            x.w <- x$weights[x.root,]
        } else {
            x.res <- x$residuals
            x.c <- x$coefficients
            x.w <- x$weights
        }

            x.n <- length(x.res)
            x.p <- length(x.c)
        
        DF.x <- x$tot.weights[x.root]*x.n - x.p 

        y.tot.sol <- y$tot.sol
        if (y.tot.sol<y.root) {
            stop(paste("'y' Root ",y.root," not found"))
        }

        if (y.tot.sol!=1) {
            y.res <- y$residuals[y.root,]
            y.c <- y$coefficients[y.root,]
            y.w <- y$weights[y.root,]
        } else {
            y.res <- y$residuals
            y.c <- y$coefficients
            y.w <- y$weights
        }

            y.n <- length(y.res)
            y.p <- length(y.c)

        DF.y <- y$tot.weights[y.root]*y.n - y.p 

    } else {

    if (inherits(x, "wle.normal") && inherits(y, "wle.normal")) {

       x.tot.sol <- x$tot.sol
        if (x.tot.sol<x.root) {
            stop(paste("'x' Root ",x.root," not found"))
        }
        if (x.tot.sol!=1) {
            x.res <- x$residuals[x.root,]
            x.w <- x$weights[x.root,]
        } else {
            x.res <- x$residuals
            x.w <- x$weights
        }
        x.n <- length(x.w)
       
        DF.x <- x$tot.weights[x.root]*x.n - 1 

        y.tot.sol <- y$tot.sol
        if (y.tot.sol<y.root) {
            stop(paste("'y' Root ",y.root," not found"))
        }

        if (y.tot.sol!=1) {
            y.res <- c(y$residuals[y.root,])
            y.w <- y$weights[y.root,]
        } else {
            y.res <- y$residuals
            y.w <- y$weights
        }
        y.n <- length(y.w)
       
        DF.y <- y$tot.weights[y.root]*y.n - 1
    
      } else {
          stop("'x' and 'y' must be of class wle.lm or wle.normal")
      }
      }
      
    V.x <- c(x.w%*%(x.res^2) / DF.x)
    V.y <- c(y.w%*%(y.res^2) / DF.y)

    ESTIMATE <- V.x / V.y
    STATISTIC <- ESTIMATE / ratio
    PARAMETER <- c(DF.x, DF.y)

    PVAL <- pf(STATISTIC, DF.x, DF.y)
    if (alternative == "two.sided") {
        PVAL <- 2 * min(PVAL, 1 - PVAL)
        BETA <- (1 - conf.level) / 2
        CINT <- c(ESTIMATE / qf(1 - BETA, DF.x, DF.y),
                  ESTIMATE / qf(BETA, DF.x, DF.y))
    }
    else if (alternative == "greater") {
        PVAL <- 1 - PVAL
        CINT <- c(ESTIMATE / qf(conf.level, DF.x, DF.y), Inf)
    }
    else
        CINT <- c(0, ESTIMATE / qf(1 - conf.level, DF.x, DF.y))
    names(STATISTIC) <- "WF"
    names(PARAMETER) <- c("num df", "denom df")
    names(ESTIMATE) <- names(ratio) <- "ratio of variances"
    attr(CINT, "conf.level") <- conf.level
    RVAL <- list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 p.value = PVAL,
                 conf.int = CINT,
                 estimate = ESTIMATE,
                 null.value = ratio,
                 alternative = alternative,
                 method = "WF test to compare two variances",
                 data.name = DNAME)
    attr(RVAL, "class") <- "htest"
    return(RVAL)
}



## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

##' Goodness-of-fit test for extreme value copulas
##' See Bernoulli 2011
##'
##' @title Goodness-of-fit test for extreme value copulas
##' @param copula a copula of the desired family
##' @param x the data
##' @param N the number of parameteric bootstrap iterations
##' @param method parameter estimation method
##' @param estimator nonparametric estimator of the Pickands dependence function
##' @param m grid size
##' @param verbose display progress bar if TRUE
##' @param print.every is deprecated
##' @param optim.method for fitCopula
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
gofEVCopula <- function(copula, x, N = 1000, method = "mpl",
                        estimator = "CFG", m = 1000,
                        verbose = TRUE, print.every = NULL,
                        optim.method = "Nelder-Mead")
{
    n <- nrow(x)
    p <- ncol(x)

    if (n < 2) stop("There should be at least 2 observations")

    if (copula@dimension != 2 || p != 2)
      stop("The copula and the data should be of dimension two")

    if (!is.null(print.every)) {
        warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
        verbose <- print.every > 0
    }

    ## make pseudo-observations
    u <- pobs(x)

    ## fit the copula
    fcop <- fitCopula(copula, u, method, estimate.variance=FALSE,
                      optim.method=optim.method)@copula

    ## where to compute A
    g <- seq(0,1-1/m,by=1/m)

    ## compute the test statistic
    s <- .C(cramer_vonMises_Afun,
            as.integer(n),
            as.integer(m),
            as.double(-log(u[,1])),
            as.double(-log(u[,2])),
            as.double(A(fcop,g)),
            stat = double(2),
            as.integer(estimator == "CFG"))$stat

    ## simulation of the null distribution
    s0 <- matrix(NA, N, 2)
    if (verbose) {
	pb <- txtProgressBar(max = N, style = if(isatty(stdout())) 3 else 1) # setup progress bar
	on.exit(close(pb)) # and close it on exit
    }
    for (i in 1:N)
    {
        u0 <- pobs(rCopula(n, fcop))

        ## fit the copula
        fcop0 <-  fitCopula(copula, u0, method, estimate.variance=FALSE,
                            optim.method=optim.method)@copula

        s0[i,] <- .C(cramer_vonMises_Afun,
                     as.integer(n),
                     as.integer(m),
                     as.double(-log(u0[,1])),
                     as.double(-log(u0[,2])),
                     as.double(A(fcop0,g)),
                     stat = double(2),
                     as.integer(estimator == "CFG"))$stat

        if (verbose) setTxtProgressBar(pb, i) # update progress bar
    }

    ## corrected version only
    structure(class = "htest",
              list(method = paste("Parametric bootstrap based GOF test for EV copulas with argument 'method' set to '",
                   method, "' and argument 'estimator' set to '", estimator, "'", sep = ""),
                   parameter = c(parameter = fcop@parameters),
                   statistic = c(statistic = s[1]),
                   p.value=(sum(s0[,1] >= s[1])+0.5)/(N+1),
                   data.name = deparse(substitute(x))))

}


### Version for simulations; was named gofEVCopula before; not exported
gofA <- function(copula, x, N = 1000, method = "mpl", # estimator = "CFG",
                    m = 1000, verbose = TRUE, optim.method = "Nelder-Mead")
{
    n <- nrow(x)
    p <- ncol(x)
    ## make pseudo-observations
    u <- pobs(x)

    ## fit the copula
    fcop <- fitCopula(copula, u, method, estimate.variance=FALSE,
                      optim.method=optim.method)@copula

    ## statistic based on Cn
    sCn <- .C(cramer_vonMises,
              as.integer(n),
              as.integer(p),
              as.double(u),
              as.double(pCopula(u,fcop)),
              stat = double(1))$stat

    ## where to compute A
    g <- seq(0,1-1/m,by=1/m)

    ## compute the CFG test statistic
    sCFG <- .C(cramer_vonMises_Afun,
               as.integer(n),
               as.integer(m),
               as.double(-log(u[,1])),
               as.double(-log(u[,2])),
               as.double(A(fcop,g)),
               stat = double(2),
               as.integer(1)            # estimator == "CFG"
               )$stat
    ## compute the Pickands test statistic
    sPck <- .C(cramer_vonMises_Afun,
               as.integer(n),
               as.integer(m),
               as.double(-log(u[,1])),
               as.double(-log(u[,2])),
               as.double(A(fcop,g)),
               stat = double(2),
               as.integer(0)# estimator == Pickands
               )$stat

    s <- c(sCn, sCFG, sPck)

    ## simulation of the null distribution
    s0 <- matrix(NA_real_, N, 5)
    if (verbose) {
	pb <- txtProgressBar(max = N, style = if(isatty(stdout())) 3 else 1) # setup progress bar
	on.exit(close(pb)) # and close it on exit
    }

    ## set starting values for fitCopula
    copula@parameters <- fcop@parameters
    for (i in 1:N)
    {
        u0 <- pobs(rCopula(n, fcop))

        ## fit the copula
        fcop0 <-  fitCopula(copula, u0, method, estimate.variance=FALSE,
                            optim.method=optim.method)@copula

        sCn0 <- .C(cramer_vonMises,
                   as.integer(n),
                   as.integer(p),
                   as.double(u0),
                   as.double(pCopula(u0, fcop0)),
                   stat = double(1))$stat

        sCFG0 <-  .C(cramer_vonMises_Afun,
                     as.integer(n),
                     as.integer(m),
                     as.double(-log(u0[,1])),
                     as.double(-log(u0[,2])),
                     as.double(A(fcop0,g)),
                     stat = double(2),
                     as.integer(1)      # estimator == "CFG"
                     )$stat
        sPck0 <-  .C(cramer_vonMises_Afun,
                     as.integer(n),
                     as.integer(m),
                     as.double(-log(u0[,1])),
                     as.double(-log(u0[,2])),
                     as.double(A(fcop0,g)),
                     stat = double(2),
                     as.integer(0)      # estimator == Pickands
                     )$stat
        s0[i,] <- c(sCn0, sCFG0, sPck0)

        if (verbose) setTxtProgressBar(pb, i) # update progress bar
    }

    list(statistic = s,
         pvalue = sapply(1:5, function(i) (sum(s0[,i] >= s[i])+0.5)/(N+1)),
         parameters = fcop@parameters)
}

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

##' Computes the sum of binomial coefficients
binom.sum <- function(n,k) sum(choose(n, 0:k))

##' Multivariate independence test based on the empirical
##' copula process as proposed by Christian Genest and Bruno
##' Rémillard (2004), Test 13:2, pages 335-369.
##'
##' @title Simulate the distribution of TA for A of cardinality 2 to m under independence
##' @param n sample size
##' @param p dimension
##' @param m up to subsets of cardinality m
##' @param N number of simulations
##' @param verbose display progress bar if TRUE
##' @param print.every is deprecated
##' @return an object of class 'indepTestDist'
##' @author Ivan Kojadinovic
indepTestSim <- function(n, p, m=p, N=1000, verbose = TRUE, print.every = NULL)
{
    if (!is.numeric(n) || (n <- as.integer(n)) < 2)
        stop("n should be an integer greater than 2")
    if (!is.numeric(p) || (p <- as.integer(p)) < 2)
        stop("p should be an integer greater than 2")
    if (!is.numeric(m) || (m <- as.integer(m)) < 2 || m > p)
        stop(paste("m should be an integer greater than 2 and smaller than",p))
    if (!is.numeric(N) || (N != as.integer(N)) || N <= 5)
	stop("N should be an integer greater than 5")
    if(N < 100 && getOption("copula:warn.idTS", TRUE))
	warning("N should be at least 100")
    if (!is.null(print.every)) {
        warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
        verbose <- print.every > 0
    }

    sb <- binom.sum(p,m)

    r. <- .C(simulate_empirical_copula,
             n = as.integer(n),
             N = as.integer(N),
             p = as.integer(p),
             m = as.integer(m),
             TA0 = double(N * (sb - p - 1)),
             G0 = double(N),
             subsets = integer(sb),
             subsets.char = character(sb),
             fisher0 = double(N),
             tippett0 = double(N),
             as.integer(verbose))

    structure(class = "indepTestDist",
              list(sample.size = n,
                   data.dimension = p,
                   max.card.subsets = m,
                   number.repetitions = N,
                   subsets = r.$subsets.char[(p+2):sb],
                   subsets.binary = r.$subsets[(p+2):sb],
                   dist.statistics.independence = matrix(r.$TA0,N,sb-p-1),
                   dist.global.statistic.independence = r.$G0,
                   dist.fisher.independence = r.$fisher0,
                   dist.tippett.independence = r.$tippett0))
}

##' Multivariate independence test based on the empirical
##' copula process as proposed by Christian Genest and Bruno
##' Rémillard (2004), Test 13:2, pages 335-369.
##'
##' @title Test in itself
##' @param x the data
##' @param d an object of class 'indepTestDist'
##' @param alpha
##' @return an object of class 'indepTest'
##' @author Ivan Kojadinovic
indepTest <- function(x, d, alpha=0.05)
{
    if (!is.numeric(x <- as.matrix(x)))
        stop("data should be numerical")
    if (any(is.na(x)))
        stop("data cannot contain missing values")
    if (nrow(x) < 2)
        stop("data should contain more than 2 rows")
    if (!inherits(d, "indepTestDist"))
	stop("d should be obtained by indepTestSim()")
    if (ncol(x) != d$data.dimension)
      stop("d was not obtained from simulations based on data whose dimension is equal to that of x")
    if (nrow(x) != d$sample.size)
      warning("d was not obtained from simulations based on data whose size is equal to that of x")
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
        stop("the significance level alpha is not properly set")

    m <- d$max.card.subsets
    p <- ncol(x)
    N <- as.integer(d$number.repetitions)

    ## transform data to ranks
    x <- apply(x,2,rank)

    sb <- binom.sum(p,m)

    ## perform test
    r. <- .C(empirical_copula_test,
             as.double(x),
             nrow(x),
             as.integer(p),
             as.integer(m),
             as.double(d$dist.statistics.independence),
             as.double(d$dist.global.statistic.independence),
             as.integer(N),
             as.integer(d$subsets.binary),
             TA = double(sb - p - 1),
             G = double(1),
             pval = double(sb - p - 1),
             fisher = double(1),
             tippett = double(1),
             globpval = double(1),
             as.double(d$dist.fisher.independence),
             as.double(d$dist.tippett.independence))

    ## compute critical values at the alpha level
    beta <- (1 - alpha)^(1 / (sb - p -1))
    critical <- numeric(0)
    from <- 1
    for (j in 2:m)
    {
        to <- from + choose(p,j) - 1
        crit <- sort(as.double(d$dist.statistics.independence[,from:to]))[round(beta * N * (to - from + 1))]
        critical <- c(critical,rep(crit,choose(p,j)))
        from <- to + 1
    }

    structure(class = "indepTest",
              list(subsets=d$subsets,statistics=r.$TA, critical.values=critical,
                   pvalues = r.$pval, fisher.pvalue=r.$fisher, tippett.pvalue=r.$tippett, alpha=alpha,
                   beta=beta, global.statistic=r.$G, global.statistic.pvalue=r.$globpval))
}


##' Serial independence test based on the empirical
##' copula process as proposed by Christian Genest and Bruno
##' Rémillard (2004), Test 13:2, pages 335-369.
##'
##' @title Simulate the distribution of the TAs for A of cardinality 2 to m
##' containing 1 under serial independence
##' @param n sample size
##' @param lag.max maximum lag
##' @param m subsets of cardinality up to m
##' @param N number of simulations
##' @param verbose display progress bar if TRUE
##' @param print.every is deprecated
##' @return an object of class 'serialIndepTestDist'
##' @author Ivan Kojadinovic
serialIndepTestSim <- function(n, lag.max, m=lag.max+1, N=1000, verbose = TRUE, print.every = NULL)
{
    if (!is.numeric(n) || (n <- as.integer(n)) < 2)
        stop("n should be an integer greater than 2")
    if (!is.numeric(lag.max) || (lag.max <- as.integer(lag.max)) < 1)
        stop("lag.max should be an integer greater than 1")

    p <- lag.max + 1 ## dimension of the virtual data

    if (n-p+1 < 2)
      stop("wrong number of lags with respect to the sample size")
    if (!is.numeric(m) || (m <- as.integer(m)) < 2 || m > p)
        stop(gettextf("m should be an integer greater than 2 and smaller than %d", p))
    if (!is.numeric(N) || (N != as.integer(N)) || N <= 10)
	stop("N should be an integer greater than 10")
    if(N < 100) warning("N should be at least 100")
    if (!is.numeric(lag.max) || (p <- as.integer(lag.max) + 1) <= 1 || n-p+1 < 2)
      stop("wrong number of lags")
    if (!is.null(print.every)) {
        warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
        verbose <- print.every > 0
    }


    sb <- binom.sum(p-1,m-1)

    R <- .C(simulate_empirical_copula_serial,
             n = as.integer(n-p+1),
             N = as.integer(N),
             p = as.integer(p),
             m = as.integer(m),
             TA0 = double(N * (sb - 1)),
             G0 = double(N),
             subsets = integer(sb),
             subsets.char = character(sb),
             fisher0 = double(N),
             tippett0 = double(N),
             as.integer(verbose))

    structure(class = "serialIndepTestDist",
	      list(sample.size = n,
		   lag.max = lag.max,
		   max.card.subsets = m,
		   number.repetitions = N,
		   subsets = R $subsets.char[2:sb],
		   subsets.binary = R $subsets[2:sb],
		   dist.statistics.independence = matrix(R $TA0,N,sb-1),
		   dist.global.statistic.independence = R $G0,
		   dist.fisher.independence = R $fisher0,
		   dist.tippett.independence = R $tippett0))
} ## end{ serialIndepTestSim }


##' Serial independence test based on the empirical
##' copula process as proposed by Christian Genest and Bruno
##' Rémillard (2004), Test 13:2, pages 335-369.
##'
##' @title Test in itself
##' @param x the data
##' @param d an object of class 'serialIndepTestDist'
##' @param alpha the asymptotic level of the test
##' @return an object of class 'indepTest'
##' @author Ivan Kojadinovic
serialIndepTest <- function(x, d, alpha=0.05)
{
    if (!is.numeric(x))
      stop("data should be numerical")
    x <- as.numeric(x)
    if (any(is.na(x)))
      stop("data cannot contain missing values")
    if (length(x) < 2)
      stop("data should contain more than 2 observations")
    if (!inherits(d, "serialIndepTestDist"))
      stop("d should result from  serialIndepTestSim()")
    n <- length(x)
    if (n != as.integer(d$sample.size))
      warning("d was not obtained from simulations based on data whose size is equal to that of x")
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
      stop("the significance level alpha is not properly set")

    ## transform data to pseudo-obs
    x <- rank(x)/n

    m <- d$max.card.subsets
    p <- as.integer(d$lag.max) + 1
    n <- n - p + 1
    N <- as.integer(d$number.repetitions)

    sb <- binom.sum(p-1,m-1)

    ## perform test
    r. <- .C(empirical_copula_test_serial,
            as.double(x),
            as.integer(n),
            as.integer(p),
            as.integer(m),
            as.double(d$dist.statistics.independence),
            as.double(d$dist.global.statistic.independence),
            as.integer(N),
            as.integer(d$subsets.binary),
            TA = double(sb - 1),
            G = double(1),
            pval = double(sb - 1),
            fisher = double(1),
            tippett = double(1),
            globpval = double(1),
            as.double(d$dist.fisher.independence),
            as.double(d$dist.tippett.independence))

    ## compute critical values at the alpha level
    beta <- (1 - alpha)^(1 / (sb -1))
    critical <- numeric(0)
    from <- 1
    for (j in 1:(m-1))
    {
        to <- from + choose(p-1,j) - 1
        crit <- sort(as.double(d$dist.statistics.independence[,from:to]))[round(beta * N * (to - from + 1))]
        critical <- c(critical,rep(crit,choose(p-1,j)))
        from <- to + 1
    }

    structure(class = "indepTest",
	      list(subsets=d$subsets,statistics=r.$TA, critical.values=critical,
		   pvalues = r.$pval, fisher.pvalue=r.$fisher, tippett.pvalue=r.$tippett, alpha=alpha,
		   beta=beta, global.statistic=r.$G, global.statistic.pvalue=r.$globpval))
}

##' Independence test among random vectors based on the empirical
##' copula process
##'
##' @title
##' @param x data
##' @param d dimension of the data
##' @param m consider subsets up to cardinality m
##' @param N number of bootstrap replicates
##' @param alpha asymptotic nominal level
##' @param verbose display progress bar if TRUE
##' @param print.every is deprecated
##' @return an object of class 'indepTest'
##' @author Ivan Kojadinovic
multIndepTest <- function(x, d, m=length(d), N=1000, alpha=0.05,
                          verbose = TRUE, print.every = NULL)
{
    if (!is.numeric(x <- as.matrix(x)))
        stop("data should be numerical")
    if (any(is.na(x)))
        stop("data cannot contain missing values")
    if ((n <- nrow(x)) < 2)
        stop("data should contain more than 2 rows")
    if (!is.numeric(d) || any((d <- as.integer(d)) < 1) || sum(d) != (nc <- ncol(x)))
      stop(paste("wrong vector of dimensions"))

    p <- length(d)

    if (!is.numeric(m) || (m <- as.integer(m)) < 2 || m > p)
        stop(paste("m should be an integer greater than 2 and smaller than",p))
    if (!is.numeric(N) || (N != as.integer(N)) || N <= 10)
	stop("N should be an integer greater than 10")
    if(N < 100) warning("N should be at least 100")
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
        stop("the significance level alpha is not properly set")
    if (!is.null(print.every)) {
        warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
        verbose <- print.every > 0
    }

    ## transform data to pseudo-observations
    x <- apply(x,2,rank)/n

    sb <- binom.sum(p,m)

    ## block indices
    b <- c(0,cumsum(d))

    ## bootstrap
    bootstrap <- .C(bootstrap_MA_I,
                    n = as.integer(n),
                    N = as.integer(N),
                    p = as.integer(p),
                    b = as.integer(b),
                    U = as.double(x),
                    m = as.integer(m),
                    MA0 = double(N * (sb - p - 1)),
                    I0 = double(N),
                    subsets = integer(sb),
                    subsets.char = character(sb),
                    as.integer(verbose))

    subsets <- bootstrap$subsets.char[(p+2):sb]
    subsets.binary <- bootstrap$subsets[(p+2):sb]
    dist.statistics.independence <- matrix(bootstrap$MA0,N,sb-p-1)

    ## perform test
    r. <- .C(empirical_copula_test_rv,
             as.double(x),
             as.integer(n),
             as.integer(p),
             as.integer(b),
             as.integer(m),
             as.double(bootstrap$MA0),
             as.double(bootstrap$I0),
             as.integer(N),
             as.integer(subsets.binary),
             MA = double(sb - p - 1),
             I = double(1),
             pval = double(sb - p - 1),
             fisher = double(1),
             tippett = double(1),
             Ipval = double(1)
	     )

    ## compute critical values at the alpha level
    beta <- (1 - alpha)^(1 / (sb - p - 1))
    critical <- numeric(sb - p - 1)

    for (k in 1:(sb - p - 1))
      critical[k] <- sort(dist.statistics.independence[,k])[round(beta * N)]

     structure(class = "indepTest",
	       list(subsets=subsets,statistics=r.$MA, critical.values=critical,
		    pvalues = r.$pval, fisher.pvalue=r.$fisher, tippett.pvalue=r.$tippett,
		    alpha=alpha, beta=beta,
		    global.statistic=r.$I, global.statistic.pvalue=r.$Ipval))
}

##' Multivariate serial independence test based on the empirical
##' copula process -- see AISM paper for more details
##'
##' @title Multivariate serial independence test based
##' @param x the data
##' @param lag.max the maximum lag
##' @param m consider subsets of cardinality up to m
##' @param N number of permutations
##' @param alpha asymptotic level of the test
##' @param verbose display progress bar if TRUE
##' @param print.every is deprecated
##' @return an object of class 'indepTest'
##' @author Ivan Kojadinovic
multSerialIndepTest <- function(x, lag.max, m=lag.max+1, N=1000, alpha=0.05,
                                verbose = TRUE, print.every = NULL)
{
    if (!is.numeric(x <- as.matrix(x)))
        stop("data should be numerical")
    if (any(is.na(x)))
        stop("data cannot contain missing values")
    if ((n <- nrow(x)) < 2)
        stop("data should contain more than 2 rows")
    if (!is.numeric(lag.max) || (p <- as.integer(lag.max) + 1) <= 1 || n-p+1 < 2)
      stop("wrong number of lags")
    if (!is.numeric(m) || ((m <- as.integer(m)) < 2) || m > p)
        stop(paste("m should be an integer greater than 2 and smaller than",p))
    if (!is.numeric(N) || (N != as.integer(N)) || N <= 10)
	stop("N should be an integer greater than 10")
    if(N < 100) warning("N should be at least 100")
    if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
        stop("the significance level alpha is not properly set")
    if (!is.null(print.every)) {
        warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
        verbose <- print.every > 0
    }

    q <- ncol(x)

    ## transform data to ranks
    x <- apply(x,2,rank)/n

    n <- n-p+1 ## number of rows of p-dimensional virtual data

    sb <- binom.sum(p-1,m-1)

    ## bootstrap
    bootstrap <- .C(bootstrap_serial,
                    n = as.integer(n),
                    N = as.integer(N),
                    p = as.integer(p),
                    q = as.integer(q),
                    U = as.double(x),
                    m = as.integer(m),
                    MA0 = double(N * (sb - 1)),
                    I0 = double(N),
                    subsets = integer(sb),
                    subsets.char = character(sb),
                    as.integer(verbose))

    subsets <- bootstrap$subsets.char[2:sb]
    subsets.binary <- bootstrap$subsets[2:sb]
    dist.statistics.independence <- matrix(bootstrap$MA0,N,sb-1)

    ## perform test
    r. <- .C(empirical_copula_test_rv_serial,
             as.double(x),
             as.integer(n),
             as.integer(p),
             as.integer(q),
             as.integer(m),
             as.double(bootstrap$MA0),
             as.double(bootstrap$I0),
             as.integer(N),
             as.integer(subsets.binary),
             MA = double(sb - 1),
             I = double(1),
             pval = double(sb - 1),
             fisher = double(1),
             tippett = double(1),
             Ipval = double(1))

    ## compute critical values at the alpha level
    beta <- (1 - alpha)^(1 / (sb - 1))
    critical <- numeric(0)
    from <- 1
    for (j in 1:(m-1))
    {
        to <- from + choose(p-1,j) - 1
        crit <- sort(as.double(dist.statistics.independence[,from:to]))[round(beta * N * (to - from + 1))]
        critical <- c(critical,rep(crit,choose(p-1,j)))
        from <- to + 1
    }

    structure(class = "indepTest",
	      list(subsets=subsets,statistics=r.$MA, critical.values=critical,
		   pvalues = r.$pval, fisher.pvalue=r.$fisher, tippett.pvalue=r.$tippett, alpha=alpha,
		   beta=beta, global.statistic=r.$I, global.statistic.pvalue=r.$Ipval))
}

##' Displays a 'dependogram' from an object of class 'indepTest'
##'
##' @title Displays a 'dependogram'
##' @param test an object of class 'indepTest'
##' @param pvalues logical; if TRUE display p-values, otherwise statistics
##' @param print additional printed information
##' @return nothing
##' @author Ivan Kojadinovic
dependogram <- function(test, pvalues=FALSE, print=FALSE)
{
  if (!inherits(test, "indepTest"))
    stop("'test' should be obtained by means of the functions indepTest, multIndepTest, serialIndepTest, or multSerialIndepTest")
  stopifnot(is.logical(pvalues), is.logical(print))

  op <- par(las=3); on.exit(par(op))

  l <- length(test$statistics)
  if (pvalues)
    {
      plot(c(1,1:l),c(0,test$pvalues), type="h", ylab="p-value per subset", xlab="",
           axes=FALSE, main="Dependogram")
      axis(1, at=1:l, labels=test$subsets)
      axis(2)
      points(1:l,rep(1 - test$beta, l),pch=20)
    }
  else ## not pvalues
    {
      plot(c(1,1:l),c(0,test$statistics), type="h", ylab="statistic per subset", xlab="",
           axes=FALSE, main="Dependogram")
      axis(1, at=1:l, labels=test$subsets)
      axis(2)
      points(1:l,test$critical.values,pch=20)
    }
  if (print)
    {
      cat("The subset statistics, p-values and critical values are:\n")
      print(data.frame(subset=test$subsets, statistic=test$statistics,
                       pvalue=test$pvalues, critvalue=test$critical.values))
      cat("The critical values are such that the simultaneous acceptance region \n",
          "has probability 1 -", test$alpha, "under the null.\n")
      cat("The individual rejection probability for any statistic obtained from the Mobius \n",
          "decomposition is 1 -", test$beta, "under the null.\n\n")
    }
}


## summary and show methods for class 'indepTest'
print.indepTest <- function(x, ...)
{
  cat("\nGlobal Cramer-von Mises statistic:", x$global.statistic,
      "with p-value", x$global.statistic.pvalue, "\n")
  cat("Combined p-values from the Mobius decomposition:\n")
  cat("  ", x$fisher.pvalue, " from Fisher's rule,\n")
  cat("  ", x$tippett.pvalue, " from Tippett's rule.\n")
  invisible(x)
}

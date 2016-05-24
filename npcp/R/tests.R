#################################################################################
##
##   R package npcp by Ivan Kojadinovic Copyright (C) 2014
##
##   This file is part of the R package npcp.
##
##   The R package npcp is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package npcp is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package npcp. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

#################################################################################
## Some change-point tests based on the empirical cdfs
#################################################################################

cpTestFn <- function(x, statistic = c("cvmmax", "cvmmean", "ksmax", "ksmean"),
                     method = c("nonseq", "seq"), b = 1,
                     weights = c("parzen", "bartlett"),
                     m = 5, L.method=c("max","median","mean","min"),
                     N = 1000, init.seq = NULL)
{
    statistic <- match.arg(statistic)
    method <- match.arg(method)
    weights <- match.arg(weights)
    L.method <- match.arg(L.method)

    stopifnot(is.matrix(x))
    d <- ncol(x)
    n <- nrow(x)
    npb <- n - 1 # number of possible breakpoints

    if (is.null(b))
        b <- bOptEmpProc(x, m=m, weights=weights,
                         L.method=L.method)
    stopifnot(b >= 1)

    ## initial standard normal sequence for generating dependent multipliers
    if (is.null(init.seq))
        init.seq <- rnorm(N * (n + 2 * (b - 1)))
    else
        stopifnot(length(init.seq) == N * (n + 2 * (b - 1)))

    ## test
    out <- .C("cptestF",
              as.double(x),
              as.integer(n),
              as.integer(d),
              cvm = double(npb),
              ks = double(npb),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(method == "seq"),
              cvm0 = double(N * npb),
              ks0 = double(N * npb),
              as.double(init.seq),
              PACKAGE = "npcp")

    ## for CVM statistics
    cvm <- out$cvm
    cvm0 <- matrix(out$cvm0,N,npb)

    ## for KS statistics
    ks <- out$ks
    ks0 <- matrix(out$ks0,N,npb)

    pval <- function(phi,s0,s)
        ( sum( apply(s0,1,phi) >= phi(s) ) + 0.5 ) / (N + 1)

    statistics <- c(cvmmax = max(cvm), cvmmean = sum(cvm)/n,
                    ksmax = max(ks), ksmean = sum(ks)/n)
    p.values <- c(cvmmax = pval(max,cvm0,cvm),
                  cvmmean = pval(sum,cvm0,cvm),
                  ksmax = pval(max,ks0,ks),
                  ksmean = pval(sum,ks0,ks))

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection based on the multivariate empirical c.d.f. with 'method'=\"%s\"", method),
                   statistic = statistics[statistic],
                   p.value =  p.values[statistic],
                   cvm = c(cvm=cvm), ks = c(ks=ks),
                   all.statistics = statistics,
                   all.p.values = p.values, b = c(b=b),
                   data.name = deparse(substitute(x))))
}


#################################################################################
## Change-point tests based on the empirical copula
#################################################################################

cpTestCn <- function(x, method = c("seq", "nonseq"), b = 1,
                     weights = c("parzen", "bartlett"), m = 5,
                     L.method=c("max","median","mean","min"),
                     N = 1000, init.seq = NULL)
{
    method <- match.arg(method)
    weights <- match.arg(weights)
    L.method <- match.arg(L.method)

    stopifnot(is.matrix(x))
    d <- ncol(x)
    stopifnot(d > 1)
    n <- nrow(x)
    npb <- n - 1 # number of possible breakpoints

    if (is.null(b))
        b <- bOptEmpProc(x, m=m, weights=weights,
                         L.method=L.method)
    stopifnot(b >= 1)

    ## initial standard normal sequence for generating dependent multipliers
    if (is.null(init.seq))
        init.seq <- rnorm(N * (n + 2 * (b - 1)))
    else
        stopifnot(length(init.seq) == N * (n + 2 * (b - 1)))

    ## test
    out <- .C("cpTestC",
              as.double(x),
              as.integer(n),
              as.integer(d),
              cvm = double(npb),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(method == "seq"),
              cvm0 = double(N * npb),
              as.double(init.seq),
              PACKAGE = "npcp")

    cvm <- out$cvm
    cvm0 <- matrix(out$cvm0,N,npb)
    statistic <- c(cvmmax=max(cvm))

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection based on the empirical copula with 'method'=\"%s\"", method),
                   statistic = statistic,
                   p.value =  ( sum( apply(cvm0,1,max) >= statistic ) + 0.5 ) / (N + 1),
                   cvm = c(cvm=cvm), b = c(b=b),
                   data.name = deparse(substitute(x))))
}

#################################################################################
## Related to the cdf of the KS statistic
## From the source of ks.test -- credit to Rcore
#################################################################################

pkolmogorov1x <- function(x, n)
{
    if (x <= 0)
        return(0)
    if (x >= 1)
        return(1)
    j <- seq.int(from = 0, to = floor(n * (1 - x)))
    1 - x * sum(exp(lchoose(n, j) + (n - j) * log(1 - x - j/n) + (j - 1) *
                        log(x + j/n)))
}

#################################################################################
## Change-point tests based on multivariate extensions of Spearman's rho
#################################################################################

cpTestRho <- function(x, method = c("mult", "asym.var"),
                      ## method = c("seq", "nonseq", "asym.var"),
                      statistic = c("pairwise", "global"),
                      ## L.method = c("pseudo","max","median","mean","min"),
                      b = 1, weights = c("parzen", "bartlett"),
                      N = 1000, init.seq = NULL)
{


    method <- match.arg(method)
    statistic <- match.arg(statistic)
    weights <- match.arg(weights)
    ## L.method <- match.arg(L.method)
    L.method <- "pseudo"

    stopifnot(is.matrix(x))
    d <- ncol(x)
    stopifnot(d > 1)
    n <- nrow(x)
    npb <- n - 1 # number of possible breakpoints

    if (is.null(b))
    {
        res <- bOptRho(x, statistic=statistic, weights=weights, L.method=L.method)
        b <- res$b
        influ <- res$influnonseq
        fbin <- res$fbin
        influest <- TRUE
    }
    else
    {
        ## f in natural order
        f <- switch(statistic,
                    global = c(rep(0,2^d - 1),1),
                    pairwise = c(rep(0, d + 1), rep(1, choose(d,2)),
                    rep(0, 2^d - choose(d,2) - d - 1)))

        ## convert f into binary order
        powerset <-  .C("k_power_set",
                        as.integer(d),
                        as.integer(d),
                        powerset = integer(2^d),
                        PACKAGE="npcp")$powerset

        fbin <- .C("natural2binary",
                   as.integer(d),
                   as.double(f),
                   as.integer(powerset),
                   fbin = double(2^d),
                   PACKAGE="npcp")$fbin

        influ <- double(n)
        influest <- FALSE
    }
    stopifnot(b >= 1)

    m <- switch(method,
                "mult" = 1, ##"seq" = 1,
                ## "nonseq" = 2,
                "asym.var" = 3)

    ## if (method %in% c("seq", "nonseq"))
    if (method == "mult")
    {
        ## initial standard normal sequence for generating dependent multipliers
        if (is.null(init.seq))
            init.seq <- rnorm(N * (n + 2 * (b - 1)))
        else
            stopifnot(length(init.seq) == N * (n + 2 * (b - 1)))
    }

    ## test
    out <- .C("cpTestRho",
              as.double(x),
              as.integer(n),
              as.integer(d),
              rho = double(npb),
              as.double(fbin),
              influ = as.double(influ),
              as.integer(influest),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(m),
              rho0 = double(N * npb),
              avar = double(1),
              as.double(init.seq),
              PACKAGE = "npcp")

    rho <- out$rho
    stat <- max(rho)

    ## if (method %in% c("seq", "nonseq"))
    if (method == "mult")
    {
        rho0 <- matrix(out$rho0,N,npb)
        p.value <- ( sum( apply(rho0,1,max) >= stat ) + 0.5 ) / (N + 1)
    }
    else
    {
        if (!(out$avar > .Machine$double.eps)) ## OK?
        {
            cat("b =",b,"\n")
            cat("influ",out$influ,"\n")
            stop("The asymptotic variance is numerically equal to zero.")
        }
        else
        {
            #if (n <= 100)
            p.value <- 2 * (1 - pkolmogorov1x(stat / sqrt(out$avar), n))
            #else
            #p.value <- 2 * exp(-2 * n * stat^2 / out$avar)
            p.value <- min(1, max(0, p.value)) ## guard as in ks.test
        }
    }

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection based on multivariate extensions of Spearman's rho with 'statistic'=\"%s\" and 'method'=\"%s\"", statistic, method),
                   statistic = c(rhomax=stat),
                   p.value = p.value,
                   rho = c(rho=rho), b = c(b=b),
                   data.name = deparse(substitute(x))))
}

#################################################################################
## Change-point tests based on U-statistics
#################################################################################

cpTestU <- function(x, statistic = c("kendall", "variance", "gini"),
                    method = c("seq", "nonseq", "asym.var"),
                    b = 1, weights = c("parzen", "bartlett"),
                    N = 1000, init.seq = NULL)
{
    method <- match.arg(method)
    statistic <- match.arg(statistic)
    weights <- match.arg(weights)

    stopifnot(is.matrix(x))
    n <- nrow(x)
    d <- ncol(x)
    if (statistic %in% c("variance", "gini"))
        stopifnot(d==1)
    if (statistic == "kendall")
        stopifnot(d > 1)

    npb <- n - 3 # number of possible breakpoints

    ## kernel
    h.func <- switch(statistic,
                     variance = function(x, y) (x - y)^2/2,
                     gini = function(x, y) abs(x - y),
                     kendall = function(x, y) prod(x < y) + prod(y < x))

    h <- matrix(0,n,n)
    for (i in seq_len(n))
        for (j in  seq_len(i))
            if (i != j)
            {
                h[i,j] <- h.func(x[i,],x[j,])
                h[j,i] <- h[i,j]
            }

    influ <- colSums(h) / (n-1) ## h1.n without centering term

    if (is.null(b))
        b <- bOptU(influ, weights=weights)

    m <- switch(method,
                "seq" = 1,
                "nonseq" = 2,
                "asym.var" = 3)

    if (method %in% c("seq", "nonseq"))
    {
        ## initial standard normal sequence for generating dependent multipliers
        if (is.null(init.seq))
            init.seq <- rnorm(N * (n + 2 * (b - 1)))
        else
            stopifnot(length(init.seq) == N * (n + 2 * (b - 1)))
    }

    ## test
    out <- .C("cpTestU",
              as.double(h),
              as.integer(n),
              as.double(influ),
              u = double(npb),
              as.integer(N),
              as.integer(weights == "bartlett"),
              as.integer(b),
              as.integer(m),
              u0 = double(N * npb),
              avar = double(1),
              as.double(init.seq),
              PACKAGE = "npcp")

    u <- out$u
    names(u) <- paste0("u",2:(n-2))
    stat <- max(u)

    if (method %in% c("seq", "nonseq"))
    {
        u0 <- matrix(out$u0,N,npb)
        p.value <- ( sum( apply(u0,1,max) >= stat ) + 0.5 ) / (N + 1)
    }
    else
    {
        if (!(out$avar > .Machine$double.eps)) ## OK?
        {
            cat("b =",b,"\n")
            cat("h1.n",influ,"\n")
            stop("The asymptotic variance is numerically equal to zero.")
        }
        else
        {
            #if (n <= 100)
            p.value <- 2 * (1 - pkolmogorov1x(stat / (2 * sqrt(out$avar)), n))
            #else
            #p.value <- 2 * exp(-2 * n * stat^2 / out$avar)
            p.value <- min(1, max(0, p.value)) ## guard as in ks.test
        }
    }

    structure(class = "htest",
              list(method = sprintf("Test for change-point detection based on U-statistics with 'statistic'=\"%s\" and 'method'=\"%s\"", statistic, method),
                   statistic = c(umax=stat),
                   p.value = p.value,
                   u = u, b = c(b=b),
                   data.name = deparse(substitute(x))))

}


#################################################################################
## NOT EXPORTED: for testing
#################################################################################

fitGEV <- function(x, method = c("pwm", "gpwm"), gamma=-0.35, delta=0,
                   landwehr = TRUE, noties = TRUE)
{
    method <- match.arg(method)
    stopifnot(is.vector(x, "numeric"))
    stopifnot(is.vector(gamma, "numeric"))
    stopifnot(is.vector(delta, "numeric"))
    n <- length(x)
    meth <- switch(method, "pwm" = 1, "gpwm" = 2)
    p <- 3 # number of statistics

    out <- .C("fitGEV",
              as.double(x),
              as.integer(n),
              as.double(gamma),
              as.double(delta),
              as.integer(meth),
              as.integer(landwehr),
              as.integer(noties),
              param = double(p),
              avar = double(p),
              PACKAGE = "npcp")

    sderr <- sqrt(out$avar/n)

    list(parameters=c(loc = out$param[1], scale = out$param[2], shape = out$param[3]),
         sderrs = c(loc = sderr[1], scale = sderr[2], shape = sderr[3]))
}

#################################################################################
## Change-point tests based on probability weighted moments
## for block maxima with the GEV in mind
#################################################################################

cpTestBM <- function(x, method = c("pwm", "gpwm"), r=10)
{
    stopifnot(is.vector(x, "numeric"))
    n <- length(x)
    r <- as.integer(r)
    stopifnot(r > 0 && n - (2 * r - 1) >= 1)
    method <- match.arg(method)

    ## internal parameters (see paper)
    r <- 10 # omit breakpoints at the beginning and the end
    gamma.var <- -0.35 # for cdf from entire sample
    delta.var <- 0 # for cdf from entire sample
    gamma.stat <- 0 # for cdfs from subsamples
    delta.stat <- 0 # for cdfs from subsamples
    landwehr <- TRUE # use Landwehr's PWM estimates if method=="pwm"
    noties <- TRUE # continous distribution assumed
    center <- TRUE # center wrt to location parameter estimate

    npb <- n - (2 * r - 1) # number of possible breakpoints


    meth <- switch(method, "pwm" = 1, "gpwm" = 2)
    p <- 3 # number of statistics

    out <- .C("cptestBM",
              as.double(x),
              as.integer(n),
              as.integer(r), # omit breakpoints at the beginning and the end
              stat = double(npb * p),
              as.double(gamma.var),
              as.double(delta.var),
              as.double(gamma.stat),
              as.double(delta.stat),
              as.integer(meth),
              as.integer(landwehr),
              as.integer(noties),
              as.integer(center),
              param = double(p),
              avar = double(p),
              PACKAGE = "npcp")

    ## statistics
    stat <- matrix(out$stat, npb, p)
    stat.loc <- stat[,1]
    names(stat.loc) <- paste0("loc",r:(n-r))
    stat.scale <- stat[,2]
    names(stat.scale) <- paste0("scale",r:(n-r))
    stat.shape <- stat[,3]
    names(stat.shape) <- paste0("shape",r:(n-r))
    maxstat.loc <- max(stat.loc)
    maxstat.scale <- max(stat.scale)
    maxstat.shape <- max(stat.shape)

    ## corresponding p-values
    if (any(out$avar <= .Machine$double.eps)) ## OK?
        stop("Some of the asymptotic variances are numerically equal to zero.")
    else
    {
        var.offset <- if (method=="pwm") c(0,10,20) else rep(0,3) ## variance offset
        out$avar <- out$avar * (var.offset + n)/n
        p.value.loc <- 2 * (1 - pkolmogorov1x(maxstat.loc / sqrt(out$avar[1]), n))
        p.value.loc <- min(1, max(0, p.value.loc)) ## guard as in ks.test
        p.value.scale <- 2 * (1 - pkolmogorov1x(maxstat.scale / sqrt(out$avar[2]), n))
        p.value.scale <- min(1, max(0, p.value.scale)) ## guard as in ks.test
        p.value.shape <- 2 * (1 - pkolmogorov1x(maxstat.shape / sqrt(out$avar[3]), n))
        p.value.shape <- min(1, max(0, p.value.shape)) ## guard as in ks.test
    }


    structure(class = "htest",
              list(method = sprintf("Test for change-point detection in the distribution of block maxima with 'method'=\"%s\"", method),
                   statistic = c(loc.stat = maxstat.loc, scale.stat = maxstat.scale, shape.stat = maxstat.shape),
                   #parameter = c(loc.param = out$param[1], scale.param = out$param[2], shape.param = out$param[3]),
                   parameter = c(p.value.loc = p.value.loc, p.value.scale = p.value.scale, p.value.shape = p.value.shape), ## improve
                   pvalues = c(p.value.loc = p.value.loc, p.value.scale = p.value.scale, p.value.shape = p.value.shape),
                   stats.loc = stat.loc,
                   stats.scale = stat.scale,
                   stats.shape = stat.shape,
                   data.name = deparse(substitute(x))))
}

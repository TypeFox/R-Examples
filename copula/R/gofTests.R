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


### Various goodness-of-fit tests ##############################################


### Test statistics ############################################################

##' Test statistics for various goodness-of-fit tests of (supposedly) U[0,1]^d
##' distributed vectors of random variates
##'
##' @title Test statistics for U[0,1]^d
##' @param u (n, d)-matrix of supposedly U[0,1]^d observations
##' @param method various test statistics. Available are:
##'        "Sn"     : the test statistic S_n (Cramer-von Mises) in Genest, Remillard, Beaudoin (2009)
##'        "SnB"    : the test statistic S_n^{(B)} in Genest, Remillard, Beaudoin (2009)
##'        "SnC"    : the test statistic S_n^{(C)} in Genest, Remillard, Beaudoin (2009)
##'        "AnChisq": Anderson-Darling test statistic after map to a chi-square
##'                   distribution using the standard normal quantile function
##'        "AnGamma": Anderson-Darling test statistic after map to an Erlang/Gamma
##'                   distribution using the logarithm
##' @param useR logical indicating whether R or C implementations are used
##' @param ... additional arguments for the different methods
##' @return n-vector of values of the chosen test statistic
##' @author Marius Hofert and Martin Maechler
gofTstat <- function(u, method=c("Sn", "SnB", "SnC", "AnChisq", "AnGamma"),
		     useR=FALSE, ...)
{
    if(!is.matrix(u)) u <- rbind(u, deparse.level=0L)
    d <- ncol(u)
    n <- nrow(u)
    method <- match.arg(method)
    switch(method,
           "Sn" =
       { ## S_n
           if(!hasArg(copula)) stop("object 'copula' required for pCopula() call")
           copula <- list(...)$copula
           ## compute \sum_{i=1}^n (\hat{C}_n(\bm{u}_i) - C_{\theta_n}(\bm{u}_i))^2
           ## typically \bm{u}_i ~> pobs \hat{\bm{U}}_i
           C. <- pCopula(u, copula=copula) # C_{\theta_n}(\bm{u}_i), i=1,..,n
           if(useR) {
               C.n <- C.n(u, U=u, ...) # \hat{C}_n(\bm{u}_i), i=1,..,n
               sum((C.n - C.)^2)
           } else {
               .C(cramer_vonMises,
                  as.integer(n),
                  as.integer(d),
                  as.double(u),
                  as.double(C.),
                  stat=double(1))$stat
           }
       },
	   "SnB" =
       { ## S_n(B)
           lu2 <- log1p(-u^2) # n x d matrix of log(1-u_{ij}^2)
           ## Note (modulo rowSums/colSums):
           ## Idea: sum1 = sum(prod(1-u^2)) = sum(exp(sum(lu2)))
           ## = exp(log( sum(exp(rowSums(lu2))) )) = exp(lsum(rowSums(lu2)))
           slu2 <- rowSums(lu2) # vector of length n
	   sum1 <- exp(lsum(cbind(slu2, deparse.level=0L))) # lsum() needs a matrix; result: 1 value
           ## 2) The notation here is similar to Genest, Remillard,
           ## Beaudoin (2009) but we interchange k and j (since j always runs
           ## in 1:d). That being said...
	   lu <- t(log1p(-u)) # t(n x d matrix of log(1-u_{ij})) --> accessing columns
           ln <- log(n)
           ## Idea:
           ##   1/n sum_i sum_k prod_j (1-max(u_{ij},u_{kj}))
           ## = 1/n sum_i sum_k exp( sum_j log(1-max{.,.}) )
           ## = 1/n sum_i sum_k exp( sum_j log(min{1-u_{ij},1-u_{kj}}) )
           ## = 1/n sum_i sum_k exp( sum_j min{ log(1-u_{ij}), log(1-u_{kj}) })
           ## = 1/n sum_i sum_k exp( sum(pmin{ lu[i,], lu[k,]}) )
           ## = 1/n sum_i exp( log(sum_k exp( sum(pmin{ lu[i,], lu[k,]}) )) )
           ## = 1/n sum_i exp( lsum( sum(pmin{ lu[i,], lu[k,]}) ) )
           ## = sum_i exp(-log(n) + lsum( sum(pmin{ lu[i,], lu[k,]}) ))
           ## = sum_i exp(-log(n) + lsum_{over k in 1:n}( sum(pmin{ lu[i,], lu[k,]}) ))
           ## => for each fixed i, (l)apply lsum()
	   sum2mands <- unlist(lapply(1:n, function(i){
	       lu.i <- lu[,i] ## now i is fixed
	       sum.k <- vapply(1:n, function(k)# sum over k (n-dim. vector)
			       sum(pmin(lu.i, lu[,k])), NA_real_)
	       ls.i <- lsum(cbind(sum.k, deparse.level=0L)) # lsum( sum(pmin(...)) ) for fixed i; 1 value
	       exp(-ln + ls.i)
	   }))
	   n/3^d - sum1/2^(d-1) + sum(sum2mands)
       },
	   "SnC" =
       { ## S_n(C)
	   Dn <- apply(u, 1, function(u.){ # Dn is a vector of length n
	       ## u. is one row. We want to know the number of rows of u
	       ## that are (all) componentwise <= u.
	       mean(apply(t(u) <= u., 2, all)) # TRUE <=> the whole row in u is <= u.
	   })
           Cperp <- apply(u, 1, prod)
	   sum((Dn-Cperp)^2)
       },
           "AnChisq" = ad.test( pchisq(rowSums(qnorm(u)^2), d) )$statistic,
	   "AnGamma" = ad.test( pgamma(rowSums(-log(u)), shape=d) )$statistic,
	   stop("unsupported method ", method))
}


### Computing different goodness-of-fit tests ##################################

##' Goodness-of-fit tests based on the parametric bootstrap
##'
##' @title Goodness-of-fit tests based on the parametric bootstrap
##' @param copula object of type 'copula' representing the H_0 copula
##'        (if necessary, parameters will be used as starting values for fitCopula())
##' @param x (n, d)-matrix containing the data
##' @param N number of bootstrap replications
##' @param method goodness-of-fit test statistic to be used; see ?gofTstat
##' @param estim.method estimation method for the unknown parameter vector; see ?fitCopula
##' @param trafo.method transformation to U[0,1]^d
##' @param verbose logical indicating whether a progress bar is shown
##' @param trafoArgs a list of optional arguments passed to the transformation method
##' @param ... additional arguments passed to \code{fitCopula}
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic, Marius Hofert
##' Note: - different '...' should be passed to different *.method
gofPB <- function(copula, x, N, method = eval(formals(gofTstat)$method),
                  estim.method = eval(formals(fitCopula)$method),
                  trafo.method = c("none", "rtrafo", "htrafo"),
		  trafoArgs = list(), verbose=TRUE, ...)
{
    ## checks
    stopifnot(is(copula, "copula"), N >= 1)
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    trafo.method <- match.arg(trafo.method)
    if(trafo.method=="htrafo") {
	if(!is(copula, "outer_nacopula"))
	    stop("trafo.method='htrafo' only implemented for copula objects of type 'outer_nacopula'")
	if(length(copula@childCops))
	    stop("currently, only Archimedean copulas are supported")
    }
    if(method=="Sn" && is(copula, "tCopula"))
        stop("The parametric boostrap with method=\"Sn\" is not available for t copulas as pCopula() cannot be computed for non-integer degrees of freedom yet.")

    ## 1) Compute the pseudo-observations
    uhat <- pobs(x)

    ## 2) Fit the copula
    C.th.n <- fitCopula(copula, uhat, method=estim.method,
			estimate.variance=FALSE, ...)@copula
    ## 3) Compute the realized test statistic
    doTrafo <- (method == "Sn" && trafo.method != "none")
    u <- if(doTrafo) {
	stopifnot(is.list(trafoArgs))
	if(length(names(trafoArgs)) != length(trafoArgs))
	    stop("'trafoArgs' must be a fully named list")
	switch(trafo.method,
	       "rtrafo"= do.call(rtrafo, c(list(uhat, cop=C.th.n), trafoArgs)),
	       "htrafo"= do.call(htrafo, c(list(uhat, cop=C.th.n), trafoArgs)),
	       stop("wrong transformation method"))
    } else uhat
    T <- if(method=="Sn") gofTstat(u, method=method, copula=C.th.n)
         else gofTstat(u, method=method)

    ## 4) Simulate the test statistic under H_0
    if(verbose) {
	pb <- txtProgressBar(max=N, style=if(isatty(stdout())) 3 else 1) # setup progress bar
	on.exit(close(pb)) # on exit, close progress bar
    }
    T0 <- vapply(1:N, function(k) {
        ## 4.1) Sample the fitted copula
        Uhat <- pobs( rCopula(n, C.th.n) )

        ## 4.2) Fit the copula
        C.th.n. <- fitCopula(copula, Uhat, method=estim.method,
                             estimate.variance=FALSE, ...)@copula
        ## 4.3) Compute the test statistic
        u. <- if(doTrafo) { # (no checks needed; all done above)
	    switch(trafo.method,
		   "rtrafo"= do.call(rtrafo, c(list(Uhat, cop=C.th.n.), trafoArgs)),
		   "htrafo"= do.call(htrafo, c(list(Uhat, cop=C.th.n.), trafoArgs)))
        } else Uhat
        T0. <- if(method=="Sn") gofTstat(u., method=method, copula=C.th.n.)
               else gofTstat(u., method=method)

        if(verbose) setTxtProgressBar(pb, k) # update progress bar
        T0. # return
    }, NA_real_)

    ## 5) return result object
    structure(class = "htest",
	      list(method = sprintf(
                   "Parametric bootstrap goodness-of-fit test with 'method'=\"%s\", 'estim.method'=\"%s\"",
		   method, estim.method),
                   parameter = c(parameter = C.th.n@parameters),
                   statistic = c(statistic = T),
                   p.value = (sum(T0 >= T) + 0.5) / (N + 1), # typical correction => p-values in (0, 1)
                   data.name = deparse(substitute(x))))
}

##' J-score \hat{J}_{\theta_n} for the multiplier bootstrap
##'
##' @title J-score \hat{J}_{\theta_n} for the multiplier bootstrap
##' @param copula object of type 'copula'
##' @param u (n, d)-matrix of (pseudo-)observations
##' @param method "mpl" or one of "itau", "irho"
##' @return n-vector containing \hat{J}_{\theta_n}
##' @author Marius Hofert (based on ideas of Ivan Kojadinovic)
##' Note: References:
##'       * For estim.method="mpl":
##'         I. Kojadinovic and J. Yan (2011), A goodness-of-fit test for multivariate
##'         multiparameter copulas based on multiplier central limit theorems, Statistics
##'         and Computing 21:1, pages 17-30
##'       * For estim.method="itau" or ="irho" (d=2):
##'         I. Kojadinovic, J. Yan and M. Holmes (2011), Fast large-sample goodness-of-fit
##'         tests for copulas, Statistica Sinica 21:2, pages 841-871
##'       * for the function in general: van der Vaart "Asymptotic Statistics" (2000, p. 179 top)
Jscore <- function(copula, u, method)
{
    ## checks
    stopifnot(is(copula, "copula"))
    if(!is.matrix(u)) u <- rbind(u, deparse.level=0L)
    stopifnot((d <- ncol(u)) > 1, (n <- nrow(u)) > 0, dim(copula) == d)

    ## deal with different methods
    switch(method,
           "mpl"=
       { ## see page 7 in Kojadinovic and Yan (2011)
           ## integrals computed from n realizations by Monte Carlo
           dcop <- dcopwrap(copula, u)
           influ0 <- derPdfWrtParams(copula, u) / dcop
           derArg <- derPdfWrtArgs  (copula, u) / dcop
           influ <- lapply(1:copula@dimension, function(i) influ0*derArg[,i]) # FIXME improve
           n <- nrow(u)
           d <- ncol(u)
           p <- length(copula@parameters)
           S <- matrix(0, n, p)
           for(j in 1:d) { # FIXME remove for() and apply()
               ij <- order(u[,j], decreasing=TRUE)
               ijb <- vapply(1:n, function(i) sum(t(u[,j]) <= u[i,j]), NA_real_)
               S <- S + rbind(rep.int(0, p),
                              apply(influ[[j]][ij,,drop=FALSE], 2, cumsum))[n+1-ijb,,drop=FALSE] / n -
                                  matrix(colMeans(influ[[j]]*u[,j]), n, p, byrow=TRUE)
           }
           solve(crossprod(influ0)/n, t(derPdfWrtParams(copula, u)/dcopwrap(copula, u)-S))
       },
           "ml"=
       {
           stop("method='ml' not available")
       },
           "itau"=
       { ## see page 849 in Kojadinovic, Yan, and Holmes (2011)
           stopifnot(dim(copula)==2)
           (4/dTau(copula)) * ( 2*pCopula(u, copula) - rowSums(u) + (1-tau(copula))/2 )
       },
           "irho"=
       { ## see Equation (3.5) on page 847 in Kojadinovic, Yan, and Holmes (2011)
           stopifnot(dim(copula)==2)
           i1 <- order(u[,1], decreasing=TRUE)
           i2 <- order(u[,2], decreasing=TRUE)
           n <- nrow(u)
           i1b <- vapply(1:n, function(k) sum(t(u[,1]) <= u[k,1]), NA_real_)
           i2b <- vapply(1:n, function(k) sum(t(u[,2]) <= u[k,2]), NA_real_)
           rmu <- rowMeans(u)
           term <- c(0, cumsum(u[,2][i1]))[n+1-i1b] / n - rmu +
               c(0, cumsum(u[,1][i2]))[n+1-i2b] / n - rmu
           (1/dRho(copula)) * ( 12*(apply(u, 1, prod) + term) - 3 - rho(copula) )
       },
           stop("wrong method"))
}

##' Goodness-of-fit tests based on the multiplier bootstrap
##'
##' @title Goodness-of-fit tests based on the multiplier bootstrap
##' @param copula object of type 'copula' representing the H_0 copula
##' @param x (n, d)-matrix containing the data
##' @param N number of bootstrap replications
##' @param method test statistics. Available are:
##'        "Sn"     : see gofTstat()
##'        "Rn"     : test statistic R_n in Genest, Huang, Dufour (2013)
##' @param estim.method estimation method for the unknown parameter vector; see
##'        ?fitCopula
##' @param verbose indicating whether verbose output is shown
##' @param useR logical indicating whether an R or the C implementation is used
##' @param ... additional arguments passed to the different methods
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic, Marius Hofert
##' Note: - Ugly C code switch
##'       - Using gofTstat() for the multiplier bootstrap seems difficult, since
##'         the realized test statistic and the bootstrapped one are computationally different
gofMB <- function(copula, x, N, method=c("Sn", "Rn"),
                  estim.method=eval(formals(fitCopula)$method),
                  verbose=TRUE, useR=FALSE, m = 1/2, zeta.m = 0, b = 0.05, ...)
{
    ## checks
    stopifnot(is(copula, "copula"), N>=1)
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)
    if(estim.method == "ml") stop("estim.method='ml' not available")
    if(estim.method %in% c("irho", "itau") && d > 2)
        stop("only bivariate case possible for estim.method='irho' or ='itau'")

    ## 1) Compute the pseudo-observations
    u. <- pobs(x)

    ## 2) Fit the copula (copula is the H_0 copula; C.th.n = C_{\theta_n}(.))
    C.th.n <- fitCopula(copula, u., method=estim.method,
			estimate.variance=FALSE, ...)@copula
    ## big, ugly switch due to C code version (although most likely not required)
    switch(method,
           "Sn"=
       {
           ## 3) Compute the realized test statistic
           T <- if(useR) { # R version
               ## compute \sum_{i=1}^n (\hat{C}_n(\hat{\bm{U}}_i) - C_{\theta_n}(\hat{\bm{U}}_i))^2
               gofTstat(u., method="Sn", useR=TRUE, copula=C.th.n)
           } else { # C version
               ## former code (MH: I don't see any reason why we should not use
               ## cramer_vonMises (via gofTstat()) here)
               ## .C(cramer_vonMises_grid,
               ##    as.integer(d),
               ##    as.double(u.),
               ##    as.integer(n),
               ##    as.double(u.),
               ##    as.integer(n),
               ##    as.double(pCopula(u., C.th.n)), # C_{\theta_n}(\hat{\bm{U}}_i), i=1,..,n
               ##    stat=double(1))$stat
               gofTstat(u., method="Sn", copula=C.th.n)
           }

           ## 4) Simulate the test statistic under H_0
           T0 <- .C(multiplier,## -> ../src/gof.c
                    p = as.integer(d),
                    U = as.double(u.),
                    n = as.integer(n),
                    G = as.double(u.),
                    g = as.integer(n),
                    influ = as.double(dCdtheta(C.th.n, u.) %*%
                                      Jscore(C.th.n, u=u., method=estim.method)),
                    N = as.integer(N),
                    s0 = double(N))$s0
       },
           "Rn"=
       {
           ## checks
	   if(estim.method!="itau")
	       stop("Currently only estim.method='itau' available for method='Rn'")
           ## 3) Compute the realized test statistic
           C.th.n. <- pCopula(u., C.th.n) # n-vector
           denom <- (C.th.n.*(1-C.th.n.) + zeta.m)^m # n-vector
           Cn. <- C.n(u., U=u.) # n-vector
           T <- sum( ((Cn. - C.th.n.)/denom)^2 ) # test statistic R_n

           ## 4) Simulate the test statistic under H_0

           ## 4.1) Draw Z's with mean 0 and variance 1
           Z <- matrix(rnorm(N*n), nrow=N, ncol=n) # (N, n)-matrix
           Zbar <- rowMeans(Z) # N-vector

           ## 4.2) Compute \hat{C}_n^{(h)}, h=1,..,N
           factor <- (Z-Zbar)/sqrt(n) # (N, n)-matrix
           ## now use a trick similar to C.n()
           ## note: - u. is an (n, d)-matrix
           ##       - t(u.)<=u.[i,] is an (d, n)-matrix
           ##       - colSums(t(u.)<=u.[i,])==d is an n-vector of logicals with kth entry
           ##         I_{\{\hat{\bm{U}}_k <= \bm{u}_i\}}
           ind <- vapply(1:n, function(i) colSums(t(u.) <= u.[i,]) == d, logical(n)) # (n, n)-matrix
           hCnh. <- factor %*% ind # (N, n)-matrix * (n, n)-matrix = (N, n)-matrix
           for(j in 1:d) { # bivariate only here
               ## build a matrix with 1s only, except for the jth column, which contains u.[,j]
               u.j <- matrix(1, nrow=n, ncol=d)
               u.j[,j] <- u.[,j]
               ## evaluate the empirical copula there
               ind.u.j <- vapply(1:n, function(i) colSums(t(u.) <= u.j[i,]) == d, logical(n)) # (n, n)-matrix
               Cnh <- factor %*% ind.u.j # (N, n)-matrix * (n, n)-matrix = (N, n)-matrix
               ## compute the "derivative" of the empirical copula at u.
               Cjn. <- dCn(u., U=u., j.ind=j, b=b) # n vector
               ## compute the sum
               Cjn.. <- matrix(rep(Cjn., each=N), nrow=N, ncol=n) # auxiliary (N, n)-matrix
               hCnh. <- hCnh. - Cjn.. * Cnh # (N, n)-matrix
           }

           ## 4.3) compute score function J and Theta
           J.th.n <- Jscore(C.th.n, u=u., method=estim.method) # n-vector
           Theta <- rowSums(Z * rep(J.th.n, each=N) / sqrt(n)) # rowSums((N, n)-matrix) = N-vector

           ## 4.4) Compute the multiplier replicates of the test statistic
           d.dth.C.th.n <- as.vector(dCdtheta(C.th.n, u=u.)) # n-vector; Note: we need to have pseudo-observations here, since for a Gumbel (H_0) copula, for example, d.dth.C=NaN for those u which have at least one component equal to 1!!!
           num <- t(hCnh. - outer(Theta, d.dth.C.th.n)) # (n, N)-matrix
           T0 <- colSums( (num/denom)^2 )/n # (n, N)-matrix => N vector; test statistic R_n^{(h)}
       },
           stop("wrong method"))

    ## return result object
    structure(class = "htest",
	      list(method = sprintf(
		   "Multiplier bootstrap goodness-of-fit test with 'method'=\"%s\", 'estim.method'=\"%s\"",
		   method, estim.method),
                   parameter = c(parameter = C.th.n@parameters),
                   statistic = c(statistic = T),
                   p.value = (sum(T0 >= T) + 0.5) / (N + 1), # typical correction =>  p-values in (0, 1)
                   data.name = deparse(substitute(x))))
}


### Wrapper ####################################################################

##' Goodness-of-fit test wrapper function
##'
##' @title Goodness-of-fit test wrapper function
##' @param copula object of type 'copula' representing the H_0 copula
##' @param x (n, d)-matrix containing the data
##' @param N the number of bootstrap (parametric or multiplier) replications
##' @param method goodness-of-fit test statistic to be used
##' @param estim.method estimation method for the unknown parameter vector
##' @param simulation parametric bootstrap ('pb') or multiplier method ('mult')
##' @param verbose logical indicating whether a progress bar is shown
##' @param print.every deprecated
##' @param ... additional arguments passed to the internal auxiliary functions
##'        gofPB() and gofMB()
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic, Marius Hofert
gofCopula <- function(copula, x, N=1000, method = "Sn",
                      estim.method = eval(formals(fitCopula)$method),
                      simulation = c("pb", "mult"),
		      verbose=TRUE, print.every=NULL, ...)
{
    ## checks
    stopifnot(is(copula, "copula"), N>=1)
    if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
    stopifnot((d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
    ## 'method' is checked inside gofPB() / gofMB()
    estim.method <- match.arg(estim.method)
    simulation <- match.arg(simulation)
    ## deprecation
    if (!is.null(print.every)) {
        warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
        verbose <- print.every > 0
    }
    ## back-compatibility
    if(missing(estim.method) && !missing(method)) {
        eMeth <- eval(formals()$estim.method)
	if(!is.na(i <- pmatch(method, eMeth))) {
	    warning("old (pre 0.999-*) argument 'method' is now called 'estim.method'")
	    estim.method <- eMeth[i]
	    method <- "Sn"
	}
    }

    ## distinguish the methods
    switch(simulation,
           "pb" = { ## parametric bootstrap
               gofPB(copula, x, N=N, method=method, estim.method=estim.method,
                     verbose=verbose, ...)
           },
           "mult" = { ## multiplier bootstrap
               gofMB(copula, x=x, N=N, method=method, estim.method=estim.method, ...)
           },
           stop("Invalid simulation method ", simulation))
}

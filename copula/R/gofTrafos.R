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


### Goodness-of-fit test transformations #######################################


### multivariate transformations ###############################################

##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to U[0,1]^d vectors of random variates
##' via Rosenblatt's transformation
##'
##' @title Rosenblatt transformation for a (nested) Archimedean copula
##' @param u data matrix in [0,1]^(n, d) ((pseudo-/copula-)observations if
##'        inverse==TRUE and U[0,1] observations if inverse==FALSE)
##' @param cop object of class Copula
##' @param j.ind index j for which C(u_j | u_1,..,u_{j-1}) is computed
##'        (in which case only this result is returned) or NULL
##'        in which case C(u_j | u_1,..,u_{j-1}) is computed for all j in
##'        {2,...,d} (in which case the non-transformed first column is returned
##'        as well)
##' @param n.MC parameter n.MC for evaluating the derivatives via Monte Carlo;
##'        if 0, the available (theoretical) formula is used.
##' @param inverse logical indicating whether the inverse of rtrafo is computed
##'        (this is known as 'conditional distribution method' for sampling)
##' @param log logical indicating whether the log-transform is computed
##' @return (n,d) matrix U supposedly U[0,1]^d (if inverse==FALSE) or 'copula'
##'         (if inverse=TRUE) realizations (if j.ind==NULL) or
##'         C(u_j | u_1,..,u_{j-1}) (if j.ind in {2,..,d}) [or the log of
##'         the result if log=TRUE]
##' @author Marius Hofert and Martin Maechler
rtrafo <- function(u, cop, j.ind=NULL, n.MC=0, inverse=FALSE, log=FALSE)
{
    ## checks
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(0 <= u, u <= 1, is(cop, "Copula"), n.MC >= 0)
    n <- nrow(u)
    d <- ncol(u)
    stopifnot(d >= 2)
    is.null.j.ind <- is.null(j.ind)
    stopifnot(is.null.j.ind || (length(j.ind)==1 && 2 <= j.ind && j.ind <= d))
    jj     <- if(is.null.j.ind) 2:d else j.ind
    jj.res <- if(is.null.j.ind) 1:d else j.ind

    ## distinguish the families (as they are quite different)

    ## (nested) Archimedean case ###############################################

    if((nac <- is(cop, "outer_nacopula")) || is(cop, "archmCopula")) {
	if(nac) {
	    if(length(cop@childCops))
		stop("currently, only Archimedean copulas are supported")
	    cop <- cop@copula
	    th <- cop@theta
	} else { # "archmCopula"
	    th <- cop@parameters
	    cop <- getAcop(cop)
	}
	stopifnot(cop@paraConstr(th))

        ## compute inverse
        if(inverse) {
            U <- u # consider u as U[0,1]^d
            max.col <- if(is.null(j.ind)) d else j.ind
            ## Clayton case (explicit)
            if(cop@name=="Clayton") {
                sum. <- U[,1]^(-th)
                for(j in 2:max.col) {
                    U[,j] <- log1p((1-j+1+sum.)*(u[,j]^(-1/(j-1+1/th))-1))/(-th)
                    eUj <- exp(U[,j])
                    sum. <- sum. + eUj^(-th)
                    if(!log) U[,j] <- eUj
                }
            } else { # non-Clayton
                ## general case
                if(n.MC) stop("The inverse of Rosenblatt's transformation with uniroot() requires n.MC=0")
                U <- u # treat u as U[0,1]^d
                for(i in 1:n)
                    for(j in j.ind)
                        U[i,j] <- uniroot(function(x)
                                          rtrafo(c(U[i,1:(j-1)], x), cop=cop, j.ind=j) - u[i,j],
                                          interval=c(0,1))$root
                if(log) U <- log(U)
            }
            return(U[,jj.res])
        }

        ## compute conditional probabilities
	psiI <- cop@iPsi(u, theta=th)	   # (n,d) matrix of psi^{-1}(u)
	psiI. <- t(apply(psiI, 1, cumsum)) # corresponding (n,d) matrix of row sums
	if(n.MC==0){
	    ## Note: C(u_j | u_1,...,u_{j-1}) = \psi^{(j-1)}(\sum_{k=1}^j \psi^{-1}(u_k)) /
            ##                                  \psi^{(j-1)}(\sum_{k=1}^{j-1} \psi^{-1}(u_k))
	    C.j <- function(j) {
                ## computes C(u_j | u_1,...,u_{j-1}) with the same idea as cacopula()
                ## but faster
		logD <- cop@absdPsi(as.vector(psiI.[,c(j,j-1)]), theta=th,
				    degree=j-1, n.MC=0, log=TRUE)
		res <- logD[1:n]-logD[(n+1):(2*n)]
		if(log) res else exp(res)
	    }
	} else { ## n.MC > 0
	    ## draw random variates
	    V <- cop@V0(n.MC, th)
	    C.j <- function(j) {
                ## computes C(u_j | u_1,...,u_{j-1}) with the same idea as
                ## default method of absdPsiMC (only difference: draw V's only once)
		arg <- c(psiI.[,j], psiI.[,j-1])
		iInf <- is.infinite(arg)
		logD <- numeric(2*n)
		logD[iInf] <- -Inf
		if(any(!iInf)) logD[!iInf] <- lsum(-V %*% t(arg[!iInf]) +
						   (j-1) * log(V) - log(n.MC))
		res <- logD[1:n]-logD[(n+1):(2*n)]
		if(log) res else exp(res)
	    }
	}
        trafo <- vapply(jj, C.j, numeric(n))
        if(is.null.j.ind) cbind(trafo, if(log) log(u[,1]) else u[,1]) else trafo

    } else if(is(cop, "normalCopula")) {

    ## Gauss copula (see, e.g., Cambou, Hofert, Lemieux) #######################

        P <- getSigma(cop)
        stopifnot(dim(P) == c(d,d)) # defensive programming

        ## compute inverse
        if(inverse) {
            U <- u # consider u as U[0,1]^d
            max.col <- if(is.null.j.ind) d else j.ind
            x <- qnorm(u[,1:max.col, drop=FALSE]) # will be updated with previously transformed U's
            for(j in 2:max.col) {
                P. <- P[j,1:(j-1), drop=FALSE] %*% solve(P[1:(j-1),1:(j-1), drop=FALSE]) # (1,j-1) %*% (j-1,j-1) = (1,j-1)
                mu.cond <- as.numeric(P. %*% t(x[,1:(j-1), drop=FALSE])) # (1,j-1) %*% (j-1,n) = (1,n) = n
                P.cond <- P[j,j] - P. %*% P[1:(j-1),j, drop=FALSE] # (1,1) - (1,j-1) %*% (j-1,1) = (1,1)
                U[,j] <- pnorm(qnorm(u[,j], mean=mu.cond, sd=sqrt(P.cond)), log.p=log)
                x[,j] <- qnorm(if(log) exp(U[,j]) else U[,j]) # based on previously transformed U's
            }
            if(log) U[,1] <- log(U[,1]) # adjust the first column for 'log'
            return(U[,jj.res])
        }

        ## compute conditional probabilities (more efficient due to vapply())
        max.col.ran <- if(is.null.j.ind) jj.res else 1:j.ind
        x <- qnorm(u[, max.col.ran, drop=FALSE])
        C.j <- function(j) {
            P. <- P[j,1:(j-1), drop=FALSE] %*% solve(P[1:(j-1),1:(j-1), drop=FALSE]) # (1,j-1) %*% (j-1,j-1) = (1,j-1)
            mu.cond <- as.numeric(P. %*% t(x[,1:(j-1), drop=FALSE])) # (1,j-1) %*% (j-1,n) = (1,n) = n
            P.cond <- P[j,j] - P. %*% P[1:(j-1),j, drop=FALSE] # (1,1) - (1,j-1) %*% (j-1,1) = (1,1)
            pnorm(x[,j], mean=mu.cond, sd=sqrt(P.cond), log.p=log)
        }
        trafo <- vapply(jj, C.j, numeric(n))
        if(is.null.j.ind) cbind(trafo, if(log) log(u[,1]) else u[,1]) else trafo

    } else if(is(cop, "tCopula")) {

    ## t Copula (see, e.g., Cambou, Hofert, Lemieux) ###########################

	P <- getSigma(cop)
        stopifnot(dim(P) == c(d,d)) # defensive programming
	nu <- getdf(cop)
        n <- nrow(u)

        ## compute inverse
        if(inverse) {
            U <- u # consider u as U[0,1]^d
            max.col <- if(is.null(j.ind)) d else j.ind
            x <- qt(u[,1:max.col, drop=FALSE], df=nu) # will be updated with previously transformed U's
            for(j in 2:max.col) {
                P1.inv <- solve(P[1:(j-1),1:(j-1), drop=FALSE])
                x1 <- x[,1:(j-1), drop=FALSE]
                g  <- vapply(1:n, function(i) x1[i, ,drop=FALSE] %*% P1.inv %*%
                             t(x1[i, ,drop=FALSE]), numeric(1))
                P.inv <- solve(P[1:j, 1:j, drop=FALSE])
                s1 <- sqrt((nu+j-1)/(nu+g))
                s2 <- (x1 %*% P.inv[1:(j-1),j, drop=FALSE]) / sqrt(P.inv[j,j])
                U[,j] <- pt((qt(u[,j], df=nu+j-1)/s1-s2) / sqrt(P.inv[j,j]),
                            df=nu, log.p=log)
                x[,j] <- qt(if(log) exp(U[,j]) else U[,j], df=nu) # based on previously transformed U's
            }
            if(log) U[,1] <- log(U[,1]) # adjust the first column for 'log'
            return(U[,jj.res])
        }

        ## compute conditional probabilities (more efficient due to vapply())
        max.col.ran <- if(is.null.j.ind) jj.res else 1:j.ind
        x <- qt(u[, max.col.ran, drop=FALSE], df=nu)
        C.j <- function(j) {
            P1.inv <- solve(P[1:(j-1),1:(j-1), drop=FALSE])
            x1 <- x[,1:(j-1), drop=FALSE]
            g  <- vapply(1:n, function(i) x1[i, ,drop=FALSE] %*% P1.inv %*%
                         t(x1[i, ,drop=FALSE]), numeric(1))
            P.inv <- solve(P[1:j, 1:j, drop=FALSE])
            s1 <- sqrt((nu+j-1)/(nu+g))
            s2 <- (x1 %*% P.inv[1:(j-1),j, drop=FALSE]) / sqrt(P.inv[j,j])
            lres <- pt(s1 * ( sqrt(P.inv[j,j]) * x[,j, drop=FALSE] + s2),
                       df=nu+j-1, log.p=TRUE)
            if(log) lres else exp(lres)
        }
        trafo <- vapply(jj, C.j, numeric(n))
        if(is.null.j.ind) cbind(trafo, if(log) log(u[,1]) else u[,1]) else trafo

    } else {
	stop("not yet implemented for copula class ", class(cop))
    }
}

##' Transforms vectors of random variates following the given (nested) Archimedean
##' copula (with specified parameters) to [0,1]^d (or [0,1]^(d-1)) vectors of
##' random variates via the transformation of Hering and Hofert (2014) or its
##' inverse
##'
##' @title Transformation of Hering and Hofert (2014) or its inverse
##' @param u data matrix in [0,1]^(n, d)
##' @param cop an outer_nacopula
##' @param include.K boolean indicating whether the last component, K, is also
##'        used (include.K = TRUE); ignored when inverse=TRUE since K is crucial
##'        there
##' @param n.MC parameter n.MC for K
##' @param inverse logical indicating whether the inverse of htrafo is computed
##'        (this transformation can also be found in Wu, Valdez, Sherris (2006)).
##' @param method method to compute qK() (see there)
##' @param u.grid argument of qK() (for method "discrete")
##' @param ... additional arguments passed to qK() (see there)
##' @return matrix of transformed realizations
##' @author Marius Hofert and Martin Maechler
htrafo <- function(u, cop, include.K=TRUE, n.MC=0, inverse=FALSE,
                   method=eval(formals(qK)$method), u.grid, ...)
{
    ## checks
    stopifnot(is(cop, "outer_nacopula"))
    if(length(cop@childCops))
        stop("currently, only Archimedean copulas are provided")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot((d <- ncol(u)) >= 2, 0 <= u, u <= 1)
    ## trafos
    th <- cop@copula@theta
    if(inverse){ # "simulation trafo" of Wu, Valdez, Sherris (2006)
        ## ingredient 1: log(psi^{-1}(K^{-1}(u_d)))
        KI <- qK(u[,d], cop=cop@copula, d=d, n.MC=n.MC, method=method,
                 u.grid=u.grid, ...) # n-vector K^{-1}(u_d)
        lpsiIKI <- cop@copula@iPsi(KI, th, log=TRUE) # n-vector log(psi^{-1}(K^{-1}(u_d)))
        n <- nrow(u)
        ## ingredient 2: sum_{k=j}^{d-1} log(u_k)/k) for j=1,..,d
        lu. <- log(u[,-d, drop=FALSE]) * (ik <- 1/rep(1:(d-1), each=n))
        cslu. <- apply(lu.[,(d-1):1, drop=FALSE], 1, cumsum) ## note: we apply cumsum to reversed columns,
        ## because we need the "upper partial sums"
        cslu. <- if(d==2) as.matrix(cslu.) else t(cslu.) # get a result of the right dimensions
        cslu <- cbind(cslu.[,(d-1):1, drop=FALSE], 0) # revert the column order, and bind last column to it
        ## => n x d matrix
        ## ingredient 3:
        l1p <- cbind(0, log1p(-u[,1:(d-1), drop=FALSE]^ ik)) # n x (d-1) matrix + dummy 0's in the first col
        ## finally, compute the transformation
        expo <- rep(lpsiIKI, d) + cslu + l1p
        cop@copula@psi(exp(expo), th)
    } else { # "goodness-of-fit trafo" of Hering and Hofert (2014)
        lpsiI <- cop@copula@iPsi(u, th, log=TRUE) # matrix log(psi^{-1}(u))
        lcumsum <- matrix(unlist(lapply(1:d, function(j)
                                        lsum(t(lpsiI[,1:j, drop=FALSE])))),
                          ncol=d)
        u. <- matrix(unlist(lapply(1:(d-1),
                                   function(k) exp(k*(lcumsum[,k]-
                                                      lcumsum[,k+1])) )),
                     ncol=d-1) # transformed components (uniform under H_0)
	if(include.K) u. <- cbind(u., pK(cop@copula@psi(exp(lcumsum[,d]), th),
                                         cop=cop@copula, d=d, n.MC=n.MC))
        u.
    }
}


### Gof wrapper (working but deprecated) #######################################

##' Conducts a goodness-of-fit test for the given H0 copula cop based on the
##' (copula) data u
##'
##' @title Goodness-of-fit testing for (nested) Archimedean copulas
##' @param u (copula-)data matrix
##' @param cop outer_nacopula with specified H0 parameters
##' @param n.bootstrap number of bootstrap replications
##' @param estim.method estimation method, see enacopula
##' @param include.K boolean whether K is included in the transformation
##' @param n.MC parameter n.MC for K
##' @param trafo multivariate goodness-of-fit transformation; available are:
##'        "Hering.Hofert" the transformation of Hering, Hofert (2014)
##'        "Rosenblatt" the transformation of Rosenblatt (1952)
##' @param method test statistic for the test of U[0,1]^d; see gofTstat()
##' @param verbose if TRUE, the progress of the bootstrap is displayed
##' @param ... additional arguments to enacopula
##' @return htest object
##' @author Marius Hofert and Martin Maechler
gnacopula <- function(u, cop, n.bootstrap,
		      estim.method=eval(formals(enacopula)$method),
		      include.K=TRUE, n.MC=0,
		      trafo=c("Hering.Hofert", "Rosenblatt"),
		      method=eval(formals(gofTstat)$method),
		      verbose=TRUE, ...)
{
    .Deprecated("gofCopula")
    gofCopula(cop, x=u, N=n.bootstrap, method=method,
              estim.method=estim.method,
              verbose=verbose, trafo.method=if(trafo == "Rosenblatt")
              "rtrafo" else "htrafo", n.MC=n.MC, if(trafo == "htrafo")
              include.K=include.K)

## working (but deprecated)

##     ## setup
##     if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
##     stopifnot(0 <= u, u <= 1, is(cop, "outer_nacopula"), (d <- ncol(u)) >= 2,
##               max(cop@comp) == d, n.bootstrap >= 0, n.MC >= 0)
##     if(length(cop@childCops))
## 	stop("currently, only Archimedean copulas are provided")
##     if(n.bootstrap == 0)
## 	stop("Choose a reasonable number of bootstrap replications or
## apply the transformations yourself,  see ?gnacopula.")
##     u.name <- deparse(substitute(u))

##     ## additional warnings for now
##     estim.method <- match.arg(estim.method)
##     if(estim.method != "mle"){
## 	if(estim.method == "smle") warning("'estim.method = \"smle\"' may be time-consuming!") else
## 	warning("Consistency for the chosen estim.method is not clear. Additionally, numerical problems might appear.")
##     }

##     ## build multivariate transformation
##     trafo <- match.arg(trafo)
##     method <- match.arg(method)
##     gtrafomulti <-
##         switch(trafo,
##                "Hering.Hofert" = {
##                    function(u, cop) htrafo(u, cop=cop, include.K=include.K, n.MC=n.MC)
##                },
##                "Rosenblatt" = {
##                    function(u, cop) rtrafo(u, cop=cop, n.MC=n.MC)
##                },
##                stop("invalid 'trafo' argument"))

##     ## build test statistic function and 'meth' string describing the method
##     meth <- paste0("Bootstrapped (B =", n.bootstrap,") test of ")
##     meth2 <- paste0(method,", est.method = ", estim.method)
##     meth <-
## 	switch(method,
##                "Sn" = {
##                    paste0(meth, meth2)
##                },
##                "SnB" =, "SnC" = {
## 		   paste0(meth, meth2," (with trafo = ", trafo, ")")
## 	       },
## 	       "AnChisq" =, "AnGamma" = {
## 		   paste0(meth, "Anderson and Darling (with trafo = ",
##                           trafo, " and method = ", meth2, ")")
## 	       },
## 	       stop("wrong 'method' argument"))

##     ## main part --- Bootstrapping ------------------

##     ## (1) estimate the parameter by the provided estimation method and
##     ##	   define the estimated copula
##     theta.hat <- enacopula(u, cop, method=estim.method, ...)
##     cop.hat <- onacopulaL(cop@copula@name, list(theta.hat, 1:d)) # copula with theta.hat

##     ## (2) transform the data with the copula with estimated parameter
##     u.prime <- gtrafomulti(u, cop=cop.hat) # transformed data in the unit hypercube

##     ## (3) conduct the Anderson-Darling test or compute the test statistic (depends on method)
##     T <- if(method=="Sn") gofTstat(u.prime, method=method, copula=cop.hat)
##          else gofTstat(u.prime, method=method)

##     ## (4) conduct the parametric bootstrap
##     theta.hat. <- numeric(n.bootstrap) # vector of estimators
##     T. <- numeric(n.bootstrap)# vector of gofTstat() results
##     if(verbose) {	     # setup progress bar and close it on exit
## 	pb <- txtProgressBar(max = n.bootstrap, style = if(isatty(stdout())) 3 else 1)
## 	on.exit(close(pb))
##     }
##     for(k in 1:n.bootstrap) {

## 	## (4.1) sample from the copula with estimated parameter and build
## 	##	     the corresponding pseudo-observations
## 	u. <- pobs(rnacopula(nrow(u), cop.hat))

## 	## (4.2) estimate the parameter by the provided method and define
## 	##	     the estimated copula
## 	theta.hat.[k] <- enacopula(u., cop, method=estim.method, ...)
## 	cop.hat. <- onacopulaL(cop@copula@name, list(theta.hat.[k], 1:d))

## 	## (4.3) transform the data with the copula with estimated parameter
## 	u.prime. <- gtrafomulti(u., cop=cop.hat.)

## 	## (4.4) compute the test statistic
##         T.[k] <- if(method=="Sn") gofTstat(u.prime., method=method, copula=cop.hat.)
##         else gofTstat(u.prime., method=method)
##         if(verbose) setTxtProgressBar(pb, k) # update progress bar
##     }

##     ## (5) build and return results
##     structure(class = "htest",
## 	      list(p.value= (sum(T. > T) + 0.5)/(n.bootstrap+1),
##                    statistic = T, data.name = u.name,
## 		   method=meth, estimator=theta.hat,
## 		   bootStats = list(estimator=theta.hat., statistic=T.)))
}

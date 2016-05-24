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


### Implementation of function evaluations and random number generation for
### nested Archimedean copulas

##' Returns the copula density at u
##'
##' @title Density of nested Archimedean copulas
##' @param x nacopula
##' @param u argument of the copula x (can be a matrix)
##' @param log if TRUE the log-density is evaluated
##' @param ... further arguments passed to the copula specific 'dacopula' function;
##'   typically  'n.MC' (Monte Carlo size), but potentially more
##' @author Marius Hofert and Martin Maechler
.dnacopula <- function(u, copula, log=FALSE, ...) {
    if(length(copula@childCops))
	stop("currently, only Archimedean copulas are supported")
    stopifnot(ncol(u) >= dim(copula)) # will be larger for children
    C <- copula@copula
    C@dacopula(u, C@theta, log=log, ...)
}
dnacopula <- function(x, u, log=FALSE, ...) {
    stopifnot(is(x, "outer_nacopula"))
    .Deprecated("dCopula")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    .dnacopula(u, x, log=log, ...)
}

##' Returns the copula density at u. Generic algorithm
##'
##' @title Density of nested Archimedean copulas (generic form)
##' @param x nacopula
##' @param u argument of the copula x
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'        to n.MC; otherwise the exact formula is used
##' @param log if TRUE the log-density is evaluated
##' @author Marius Hofert
dnacopulaG <- function(x, u, n.MC=0, log = FALSE) {
    stopifnot(is(x, "outer_nacopula"))
    if(length(x@childCops))
	stop("currently, only Archimedean copulas are supported")
    dacopulaG(x@copula, u=u, n.MC=n.MC, log=log)
}

##' Returns the copula density at u. Generic algorithm
##'
##' @title Density of Archimedean copulas (Generic form)
##' @param acop acopula
##' @param u argument of the copula x
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'        to n.MC; otherwise the exact formula is used
##' @param log if TRUE the log-density is evaluated
##' @author Martin Maechler, based on Marius' code
dacopulaG <- function(acop, u, n.MC=0, log = FALSE) {
    stopifnot(is(acop, "acopula"), 0 <= u, u <= 1)
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    if((d <- ncol(u)) < 2) stop("u should be at least bivariate")
    th <- acop@theta
    res <- rep(NaN,n <- nrow(u)) # density not defined on the margins
    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
    if(any(n01)) {
        u. <- u[n01,, drop=FALSE]
	psiI <- acop@iPsi(u.,th)
	psiID <- acop@absdiPsi(u.,th)
        res[n01] <- acop@absdPsi(rowSums(psiI),theta = th,degree = d, n.MC = n.MC, log = log)
        res[n01] <- if(log) res[n01] + rowSums(log(psiID)) else res[n01] * apply(psiID,1,prod)
    }
    res
}

##' Returns the copula value at u
##'
##' @title Evaluation of nested Archimedean copula
##' @param u argument of the copula x (can be a matrix)
##' @param copula nacopula
##' @return f_x(u)
##' @author Marius Hofert, Martin Maechler
.pnacopula <- function(u, copula, ...) {
    stopifnot(ncol(u) >= dim(copula)) # will be larger for children
    C <- copula@copula
    th <- C@theta
    C@psi(rowSums(## use u[,j, drop=FALSE] for the direct components 'comp':
		  cbind(C@iPsi(u[,copula@comp, drop=FALSE], theta=th),
			## and recurse down for the children:
			C@iPsi(unlist(lapply(copula@childCops, .pnacopula, u=u)),
			       theta=th))),
	  theta=th)
}
pnacopula <- function(x,u) {
    .Deprecated("pCopula")
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    .pnacopula(u,x)
}


##' Returns the copula value at u
##'
##' @title CDF / Evaluation of Archimedean copula
##' @param u argument of the copula x (can be a matrix)
##' @param C acopula
##' @param theta parameter
##' @return C_\theta(u)
##' @author Martin Maechler
.pacopula <- function(u, C, theta = C@theta)
    C@psi(rowSums(C@iPsi(u, theta=theta)), theta=theta)
pacopula <- function(u, C, theta = C@theta) {
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    .pacopula(u,C,theta)
}

##' Compute the probability P[l < U <= u]  where U ~ copula x
##'
##' @title Compute the probability P[l < U <= u]  where U ~ copula x
##' @param x copula object
##' @param l d-vector of lower "integration" limits
##' @param u d-vector of upper "integration" limits
##' @return the probability that a random vector following the given copula
##'         falls in the hypercube with lower and upper corner l and u, respectively.
##' @author Marius Hofert, Martin Maechler
setGeneric("prob", function(x, l, u) standardGeneric("prob"))

setMethod("prob", signature(x ="Copula"),
          function(x, l,u) {
              d <- dim(x)
              stopifnot(is.numeric(l), is.numeric(u),
                        length(u) == d, d == length(l),
                        0 <= l, l <= u, u <= 1)
              if(d > 30)
		  stop("prob() for copula dimensions > 30 are not supported (yet)")
              D <- 2^d
              m <- 0:(D - 1)
              ## digitsBase() from package 'sfsmisc' {slightly simplified} :
              ## Purpose: Use binary representation of 0:N
              ## Author: Martin Maechler, Date:  Wed Dec  4 14:10:27 1991
              II <- matrix(0, nrow = D, ncol = d)
              for (i in d:1L) {
                  II[,i] <- m %% 2L + 1L
                  if (i > 1) m <- m %/% 2L
              }
              ## Sign: the ("u","u",...,"u") case has +1; = c(2,2,...,2)
              Sign <- c(1,-1)[1L + (- rowSums(II)) %% 2]
              U <- array(cbind(l,u)[cbind(c(col(II)), c(II))], dim = dim(II))
              sum(Sign * pCopula(U, x))
          })

##' Returns (n x d)-matrix of random variates
##'
##' @title Random number generation for nested Archimedean copulas
##' @param n number of vectors of random variates to generate
##' @param copula outer_nacopula
##' @param x  (for back compatibility)
##' @param ...
##' @return matrix of random variates
##' @author Marius Hofert, Martin Maechler
rnacopula <- function(n, copula, x, ...)
{
    if(!missing(x)) {
	if(missing(copula)) {
	    warnings("the argument 'x' has been renamed to 'copula' and is deprecated")
	    copula <- x
	}
	stop("cannot specify both 'copula' and 'x'")
    }
    Cp <- copula@copula			# outer copula
    theta <- Cp@theta			# theta for outer copula
    V0 <- Cp@V0(n,theta)		# generate V0's
    childL <- lapply(copula@childCops, rnchild, # <-- start recursion
		     theta0=theta,V0=V0,...)
    dns <- length(copula@comp)	 # dimension of the non-sectorial part
    r <- matrix(rexp(n*dns), n, dns) # generate the non-sectorial part
    ## put pieces together
    mat <- Cp@psi(r/V0, theta=theta)	# transform
    mat <- cbind(mat, do.call(cbind,lapply(childL, `[[`, "U")))
    ## get correct sorting order:
    j <- c(copula@comp, unlist(lapply(childL, `[[`, "indCol")))
    ## extra check
    stopifnot(length(j) == ncol(mat))
    mat[, order(j), drop=FALSE] # permute data and return
}

##' Returns (n x d)-matrix of random variates
##'
##' @title Random number generation for Archimedean copulas
##' @param n number of vectors of random variates to generate
##' @param Cp acopula
##' @param theta copula parameter
##' @param d dimension
##' @return matrix of random variates
##' @author Martin Maechler
racopula <- function(n, Cp, theta = Cp@theta, d)
{
    stopifnot(n == as.integer(n), d == as.integer(d),
	      n >= 0, d >= 2, is.finite(theta))
    V0 <- Cp@V0(n, theta)	# generate V0's
    E <- matrix(rexp(n*d), n, d)# .. the exponentials
    Cp@psi(E/V0, theta=theta)
}


##' Returns a list with an (n x d)-matrix of random variates and a vector of
##' indices.
##'
##' @title Random number generation for children of nested Archimedean copulas
##' @param x nacopula
##' @param n number of vectors of random variates to generate
##' @param theta0 parameter theta0
##' @param V0 vector of V0's
##' @return list(U = matrix(*,n,d), indCol = vector of length d)
##' @author Marius Hofert, Martin Maechler
rnchild <- function(x, theta0, V0,...)
{
    n <- length(V0)
    Cp <- x@copula # inner copula
    ## Consistency checks -- for now {comment later} :
    stopifnot(is(Cp, "acopula"), is.numeric(n), n == as.integer(n),
              is.numeric(V0), length(V0) == n, is.numeric(theta0))
    theta1 <- Cp@theta # theta_1 for inner copula
    ## generate V01's (only for one sector since the
    ## recursion in rnacopula() takes care of all sectors):
    V01 <- Cp@V01(V0, theta0,theta1,...)
    childL <- lapply(x@childCops, rnchild, # <-- recursion
                     theta0=theta1, V0=V01,...)
    dns <- length(x@comp)	# dimension of the non-sectorial part
    r <- matrix(rexp(n*dns), n, dns) # generate the non-sectorial part
    ## put pieces together: first own comp.s, then the children's :
    mat <- Cp@psi(r/V01, theta1) # transform
    if(length(childL) && length(U <- lapply(childL, `[[`, "U")))
	mat <- cbind(mat, do.call(cbind, U))
    ## get correct sorting order:
    j <- c(x@comp, unlist(lapply(childL, `[[`, "indCol")))
    list(U = mat, indCol = j) # get list and return
}

if(FALSE) { # evaluate the following into your R session if you need debugging:
    trace(rnacopula, browser, exit=browser, signature=signature(x ="outer_nacopula"))

    debug(rnchild)
}

##' @title Constructor for outer_nacopula
##' @param family either character string (short or longer form of
##'	 copula family name) or an "acopula" family object
##' @param nacStructure a "formula" of the form C(th, comp, list(C(..), C(..)))
##' @return a valid outer_nacopula object
##' @author Martin Maechler
onacopula <- function(family, nacStructure) {
    ## , envir = ... , enclos=parent.frame() or
    ## , envir=environment()
    ##
### FIXME: base this on  onacopulaL() -- replacing nacStructure by nacList
### ----- : (1) replacing  C() with list() *AND* by
###         (2) wrapping the 3rd argument with list(.) if it's not already
### Use a *recursive* function like this one (but add the "list(.)" wrapping:
###    repC <- function(e) { sC <- as.symbol("C"); if(identical(e[[1]], sC)) e[[1]] <- as.symbol("list"); if(length(e) == 4) e[[4]] <- repC(e[[4]]); e}
    nacl <- substitute(nacStructure)
    C. <- as.symbol("C")
    if(!is.call(nacl) || !length(nacl) || !identical(nacl[[1]], C.)) {
        ## assume nacStructure is rather a call or expression *object*
        nacl <- if(is.expression(nacStructure)) nacStructure[[1]]
                else if(is.call(nacStructure)) nacStructure # else NULL
        if(!length(nacl) || !identical(nacl[[1]], C.))
            stop("invalid 'nacStructure'.  Maybe use 'onacopulaL()' which is more robust albeit more verbose")
    }
    COP <- getAcop(family)
    nacl[[1]] <- as.symbol("oC")
### does not work ..>>>>>>>>>>> we should use onacopulaL() inside functions! <<<<
    ## needed, e.g., when onacopula() is called from inside a user function:
    ## for(j in 2:3) if(is.language(nacl[[j]]))
    ##     nacl[[j]] <- eval(nacl[[j]], parent.frame())
    ## pframe <- parent.frame()
    mkC <- function(cClass, a,b,c) {
	if(missing(b) || length(b) == 0) b <- integer()
	if(missing(c) || length(c) == 0) c <- list()
	else if(length(c) == 1 && !is.list(c)) c <- list(c)
### does not work ...
	## else if(!is.list(c)) {
        ##     c <- eval(c, pframe) ; stopifnot(is.list(c))
        ## }
	## if(!is.numeric(a)) {
        ##     a <- eval(a, pframe) ; stopifnot(is.numeric(a), length(a) == 1)
        ## }
	## if(!is.numeric(b)) {
        ##     b <- eval(b, pframe) ; stopifnot(is.numeric(b))
        ## }
	else stopifnot(is.list(c))
	stopifnot(is.numeric(a) || is.na(a), length(a) == 1, is.numeric(b))
	if(any(sapply(c, class) != "nacopula"))
	    stop("third entry of 'nacStructure' must be NULL or list of 'C(..)' terms")
	new(cClass, copula = setTheta(COP, a),
	    comp = as.integer(b), childCops = c)
    }
    C <- function(a,b,c) mkC("nacopula", a,b,c)
    oC <- function(a,b,c) mkC("outer_nacopula", a,b,c)
    eval(nacl)
}

##' @title Constructor for outer_nacopula - with list() input and *recursively*
##' @param family
##' @param nacList
##' @return
##' @author Martin Maechler
onacopulaL <- function(family, nacList) {
    COP <- getAcop(family)
    mkC <- function(cClass, abc) {
        stopifnot(is.list(abc), 2 <= (nL <- length(abc)))
        if(nL == 2) abc[[3]] <- list() else
        if(nL > 3) stop("nacLists must be of length 3 (or 2)")
	new(cClass, copula = setTheta(COP, abc[[1]]),
	    comp = as.integer(abc[[2]]),
            childCops = lapply(abc[[3]], function(.) mkC("nacopula", .)))
    }
    mkC("outer_nacopula", nacList)
}

##' @title The *inverse* of onacopulaL
##' @param x an (outer) nacopula
##' @return a list like the 'nacList' argument of onacopulaL()
##' @author Martin Maechler
nac2list <- function(x) {
    stopifnot(is(x, "nacopula"))
    if(length(kids <- x@childCops))# recurse
        list(x@copula@theta, x@comp, lapply(kids, nac2list))
    else list(x@copula@theta, x@comp)
}

##'*randomly* construct a nested Archimedean copula model
##' @title Random nacopula Model
##' @param family the Archimedean family
##' @param d integer >=2; the dimension
##' @param pr.comp probability of a direct component on each level
##' @param rtau0 a \code{\link{function}} to generate a (random) tau,
##'   corresponding to theta0, the outermost theta.
##' @return an 'outer_nacopula'
##' @author Martin Maechler,  10 Feb 2012
rnacModel <- function(family, d, pr.comp, rtau0 = function() rbeta(1, 2,4),
                      order=c("random", "each", "seq"), digits.theta = 2)
{
### TODO:  1) depth.max [default "Inf" = d]
### ----   2) distribution instead of uniform {1,..., n%/% 2}  for nkids

    COP <- getAcop(family)
    stopifnot(length(d) == 1, d == round(d), d >= 2,
	      length(pr.comp) == 1, 0 < pr.comp, pr.comp < 1,
	      is.function(rtau0),
	      length(digits.theta) == 1, is.numeric(digits.theta))
    ## Do: construct a random nacList  and call onacopulaL()
    ## Care: theta values must obey family specific nesting constraints
    ##	  we use paraInterval +	 theta1 >= theta0  [ok for the 5 families]

    RR <- if(do.round <- digits.theta < 17)
	function(t) round(t, digits.theta) else identity
    stopifnot(is.numeric(t0 <- rtau0()), length(t0) == 1, abs(t0) <= 1)
    theta0 <- RR(COP@iTau(t0))
    bound.th <- as.numeric(COP@paraInterval)[2]
    rtheta <- { ## function to generate child copula thetas
	if(bound.th == Inf) function(min.th) # in [min.th, Inf)
	    RR(1/runif(1, 0, 1/min.th))
       else function(min.th) RR(runif(1, min.th, bound.th))
    }

    rComp <- switch(match.arg(order),
                    "random" = function(n, size)      sample.int(n, size=size),
                    "each"   = function(n, size) sort(sample.int(n, size=size)),
                    "seq"    = function(n, size) seq_len(size))
    ## mkL(): The function to be called recursively
    ## cs (originally = 1:d) := the set of component IDs to choose from
    mkL <- function(cs, theta) {
	nc <- length(cs)
        if(nc < 2) {
            stop("nc = ",nc," < 2 -- should not happen")
        } else if(nc == 2) {
            ncomp <- 2
            comp <- cs
            cs <- integer(0)
            nc <- length(cs)
        } else { ## nc >= 3
            ncomp <- rbinom(1, size=nc, prob = pr.comp)
            if(nc == 3) ## two cases: either 3 comp, or 1 comp + 1 childcop
                if(ncomp == 0) ncomp <- 1
            if(ncomp) {
                if(ncomp+1 == nc) ## cannot leave only one for children
                    ncomp <- ncomp+1
                ii <- rComp(nc, size=ncomp)
                comp <- cs[ ii]
                cs   <- cs[-ii]
                nc <- length(cs)
            } else comp <- integer(0)   # and keep cs
        }
	if(nc) {
	    if(nc < 2) stop("internal error - nc=", nc, " < 2")
	    ## choose number of child copulas (>= 1)
	    nkids <- sample.int(nc %/% 2, size=1)
	    if(!ncomp && nkids == 1)## if no 'comp's, at least *two* child copulas:
		nkids <- 2
	    ## now, at least two components (from cs) must go to each "kid"
	    if(nkids * 2 > nc)
		stop("internal error: 'nkids * 2 > nc': nkids=",nkids)
	    nc. <- nc - 2*nkids
	    n.each <- 2 + as.vector(rmultinom(1, size= nc.,
					      prob=rep.int(1/nkids, nkids)))
	    stopifnot(sum(n.each) == nc)
            ## i.each: The coordinate IDs for each child copula
	    i.each <- as.list(n.each)
	    ind <- cs # the indices from which to choose
	    for(i in seq_along(i.each)) {
                ii <- rComp(length(ind), size=n.each[i])
		i.each[[i]] <- ind[ii]
		ind <- ind[-ii]
	    }
	    ## now recurse for each child copula
	    LL <- lapply(i.each, function(ie)
			 mkL(ie, theta= rtheta(min.th = theta)))
	    list(theta, comp, LL)
	}
	else list(theta, comp)
    }## mkL()

    nacl <- mkL(cs = 1:d, theta0)
    onacopulaL(COP, nacl)
}

##' @title Pairwise thetas of Nested Archimedean Copula
##' @param x an (outer) nacopula (with thetas sets)
##' @return matrix [d x d] of thetas, T,  T[j,k] = theta of C(U_j,U_k)
##' @author Martin Maechler
##' @note This version is simple (to program) but "not so pretty":
##'  1) .allComp() recurses every time
##'  2) we overwrite some matrix element several times ..
##' --> instead using the *longer* nacPairthetas() below
nacPairthetas0 <- function(x) {
    stopifnot(is(x, "nacopula"))
    d <- dim(x)
    T <- matrix(NA_real_, d,d)

    ## setT() : Set thetas for "sub copula" x
    ## x : nacopula
    ## iB: subset {1:d}, the components *below*
    setT <- function(x, iB) {
	T[iB, iB] <<- x@copula@theta
	for(kk in x@childCops)# recurse
	    setT(kk, .allComp(kk))
    }
    setT(x, 1:d)
    diag(T) <- rep.int(NA_real_, d)
    T
}## nacPairthetas0

##' @title Pairwise thetas of Nested Archimedean Copula
##' @param x an (outer) nacopula (with thetas sets)
##' @return matrix [d x d] of thetas, T,  T[j,k] = theta of C(U_j,U_k)
##' @author Martin Maechler
nacPairthetas <- function(x) {
    stopifnot(is(x, "nacopula"))
    d <- dim(x)
    T <- matrix(NA_real_, d,d)
    ## setT() : Set thetas for "sub copula" x; return its 'comp'
    setT <- function(x) {
        ii <- x@comp
        th <- x@copula@theta
        T[ii,ii] <<- th
	for(kk in x@childCops) {# recurse
	    ik <- setT(kk)
            T[ii,ik] <<- T[ik,ii] <<- th
            ii <- c(ii,ik)
        }
        ii
    }
    setT(x)
    diag(T) <- rep.int(NA_real_, d)
    T
}## nacPairthetas

###-- methods - glue  former "copula" <--> former "nacopula" ---------

setMethod("pCopula", signature("matrix", "nacopula"), .pnacopula)
setMethod("pCopula", signature("numeric", "nacopula"),
	  function(u, copula, ...)
	  .pnacopula(rbind(u, deparse.level = 0L), copula))

setMethod("dCopula", signature("matrix", "nacopula"), .dnacopula)
setMethod("dCopula", signature("numeric", "nacopula"),
	  function(u, copula, log=FALSE, ...)
	  .dnacopula(rbind(u, deparse.level = 0L), copula, log=log))

setMethod("rCopula", signature("numeric", "nacopula"), rnacopula)

setMethod("tailIndex", "acopula",
	  function(copula, ...) {
	      th <- copula@theta
	      if(any(is.na(th)))
		  warning("'theta' is NA -- maybe rather apply to setTheta(.)")
	      c(copula@lambdaL(th), copula@lambdaU(th))
	  })
setMethod("tailIndex", "nacopula", function(copula, ...) tailIndex(copula@copula, ...))

setMethod("tau", "acopula", function(copula) copula@tau(copula@theta))
setMethod("tau", "nacopula", function(copula) tau(copula@copula))

setMethod("rho", "nacopula", function(copula) rho(copula@copula))
setMethod("rho", "acopula", function(copula)
    stop(gettextf("%s() method for class \"%s\" not implemented;",
                                    "rho", class(copula)),
         "\nconsider contacting  maintainer(\"copula\")")
)

setMethod("iTau", "acopula", function(copula, tau, ...) copula@iTau(tau))


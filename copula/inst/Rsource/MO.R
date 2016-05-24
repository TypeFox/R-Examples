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


### Sampling of MO and estimation of MO and MSMVE distributions ################


### Sampling ###################################################################

##' @title Sampling algorithm for the Marshall Olkin distribution / copula
##' @param n sample size
##' @param lambda vector of length 2^d-1
##' @param method indicating whether copula data or data from the multivariate
##'        distribution is generated
##' @return (n, d)-matrix with samples from a Marshall-Olkin distribution
##' @author Marius Hofert
##' @note Algorithm based on Arnold (1975) "A characterization of the exponential
##'       distribution by multivariate geometric compounding";
##'       see Lemma 3.4 and Algorithm 3.3 in Mai, Scherer (2012, p. 114)
##'       Idea:
##'       1) draw Y_1, Y_2, ... ~ iid in {1,..,2^d-1} (indicates a line
##'          in J and thus a subset) according to probabilities p (s. below)
##'       2) Draw E_1, E_2, ... ~ iid Exp(sum(lambda))
##'       3) Define X_j = sum_{k=1}^{min{i: j in J[Y_i,]}} E_k
##'       => (X_1,..,X_d) ~ MO distribution(lambda)
rMO <- function(n, lambda, method=c("copula", "mvd"))
{
    ## checks
    stopifnot(n >= 1, is.numeric(lambda), lambda >= 0)
    method <- match.arg(method)
    D <- length(lambda)
    d <- log(D+1, base=2) # D = 2^d-1
    if((d != as.integer(d)) || (d < 2)) stop("length(lambda) must be 2^d-1 for d >= 2")

    ## compute set J of indices of all non-empty subsets of {1,..,d}
    J <- expand.grid(rep(list(c(FALSE, TRUE)), d))[-1,] # remove row with all FALSE
    rownames(J) <- NULL; names(J) <- paste("Idx", 1:d, sep="")

    ## Step 1
    rate <- apply(J, 2, function(Jk) sum(lambda[Jk]))

    ## Step 2
    lambdaSum <- sum(lambda)
    p <- lambda/lambdaSum # ... corresponding to J

    ## Step 3 + 4
    U <- matrix(, nrow=n, ncol=d)
    for(i in seq_len(n)) { # unfortunately inefficient
        E <- 0
        tofill <- rep(TRUE, d)
        while(any(tofill)) {
            ## (randomly) pick out a line in J (according to probs in p)
            ii <- J[sample(seq_len(D), size=1, replace=TRUE, prob=p),] # J[Y,]
            idx <- ii & tofill # cols of U[i,] we work on / fill
            E <- E + rexp(1, lambdaSum)
            U[i,idx] <- rep(E, sum(idx))
            if(method=="copula") U[i,idx] <- exp(-rate[idx]*U[i,idx]) # apply marginal survival functions
            tofill <- tofill & !idx # adjust tofill
        }
    }
    U
}


### Estimation #################################################################

##' @title Parameter estimation for Marshall-Olkin distributions
##' @param x (n, d)-matrix with n >= 2, d >= 2
##' @return (2^d-1, d+1)-matrix with the first d columns indicating
##'         the subset J of the estimated lambda_J and the last column
##'         containing lambda_J
##' @author Marius Hofert
fitMO <- function(x)
{
    ## checks
    stopifnot(length(dim(x)) == 2L,
              (d <- ncol(x)) >= 2L, (n <- nrow(x)) >= 2L)

    ## compute estimators Lambda_J = hat(l)_{theta,n}(J)
    ## compute set J of indices of all non-empty subsets of {1,..,d}
    J <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), d))[-1,]) # remove row with all FALSE
    rownames(J) <- NULL; colnames(J) <- paste("Idx", 1:d, sep="")
    ## compute Lambda_J corresponding to a non-empty subset J of {1,..,d}
    Lam <- 1 / apply(J, 1, function(ii) mean(apply(x[,ii, drop=FALSE], 1, min))) # 2^d-1 vector
    ## => comparably expensive (!)

    ## Compute the inverse Ainv of the (2^d-1, 2^d-1)-matrix A with (i,j) element 1
    ## iff the set corresponding to J[i,] has a non-empty intersection with the set
    ## corresponding to J[j,]
    ## Note: - for checking purposes: Ainv <- solve(pmin(J %*% t(J), 1)),
    ##         but we have an explicit formula for Ainv here
    ##       - this part does *not* depend on 'n' anymore
    Ainv <- (-1)^(J %*% t(J) + 1) # J %*% t(J) = J_i cap J_j
    cupJ <- function(i, j) sum(J[i,] | J[j,]) < d # = 1 <=> |J_i cup J_j| < d
    D <- 2^d-1
    Ainv[outer(1:D, 1:D, FUN=Vectorize(cupJ))] <- 0 # => Ainv = A^{-1}
    cbind(J, lambda = as.vector(Ainv %*% Lam)) # solve A lambda = Lambda w.r.t. lambda
}

##' @title Parameter estimation for stable tail dependence functions in MSMVEs and EVCs
##' @param x (n, d)-matrix with n >= 2, d >= 2
##' @param stdf function(t, theta) which returns the value l_theta(t) for a matrix t
##' @param weights d-vector of weights (corresponding to non-empty subsets of same size)
##' @param interval interval where to optimize stdf in the 1d-case or used as 'lower'
##'        and 'upper' in the 'general d case'
##' @param start initial value for the optimization of stdf with optim() for general d
##' @param verbose logical indicating whether verbose messages are printed
##' @param ... additional arguments passed to optimize() (one-dimensional case)
##'        or optim() (general case)
##' @return estimate (return value of optim()) [typically 'init']
##' @author Marius Hofert
fitStdf <- function(x, stdf, weights, interval, start, verbose=TRUE, ...)
{
    ## checks
    stopifnot(length(dim(x)) == 2L,
              (d <- ncol(x)) >= 2L, (n <- nrow(x)) >= 2L,
              is.function(stdf),
              length(weights) == d)
    if(verbose && d > 10) warning("Note that d > 10 can become quite time consuming (order O(2^d))")

    ## compute estimators Lambda_J = hat(l)_{theta,n}(J)
    ## compute set J of indices of all non-empty subsets of {1,..,d}
    J <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), d))[-1,]) # remove row with all FALSE
    rownames(J) <- NULL; colnames(J) <- paste("Idx", 1:d, sep="")
    ## compute Lambda_J corresponding to a non-empty subset J of {1,..,d}
    Lam <- 1 / apply(J, 1, function(ii) mean(apply(x[,ii, drop=FALSE], 1, min))) # 2^d-1 vector

    ## main (choose equal weights for sets J of equal size)
    f <- function(th) sum(weights[rowSums(J)] * (stdf(J, th) - Lam)^2)
    has.interval <- !missing(interval)
    has.start <- !missing(start)
    if(has.start & has.interval) # general case with boundaries
        optim(par=start, fn=f, lower=interval[1], upper=interval[2], ...)
    else if(has.start) # general case
        optim(par=start, fn=f, ...)
    else if(has.interval) # one-dimensional case
        optimize(f=f, interval=interval, ...)
    else stop("Either 'interval' or 'start' have to be provided")
}


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


### Densities of two-level nested Archimedean copulas ##########################

## Reference {in *.Rd style}:
##     Marius Hofert and David Pham (2013).
##     Densities of nested Archimedean copulas.
##     \emph{Journal of Multivariate Analysis}, \bold{118}, 37--52.
## doi = http://dx.doi.org/10.1016/j.jmva.2013.03.006
##
## --> http://arxiv.org/abs/1204.2410 is the preprint of the above, see for more details


### functions for the likelihood ###############################################

##' Compute the coefficients a_{s, d_s k}(t) for proper child copulas
##'
##' @title Compute the coefficients a_{s, d_s k}(t) for all t and k in 1:d_s
##' @param t vector t_s(u_s)
##' @param cops sector s copula list (of type nacopula with no children)
##' @param th0 parameter theta of the root copula
##' @return n x d_s matrix containing a_{s, d_s k}(t) = for all t and all k in 1:d_s
##' @author Marius Hofert
##' Note: s_{nk}(x) = \sum_{l=k}^n s(n,l)S(l,k)x^l = (-1)^{n-k} copula:::coeffG(n, x)
##'       but we need s_{nk}(x) vectorized in k and (possibly; for AMH) x
a.coeff <- function(t, cops, th0)
{
    stopifnot(is(cops, "nacopula"), length(cops@childCops)==0)
    ds <- dim(cops) # dimension d_s of the sector copula cops
    k. <- 1:ds
    ths <- cops@copula@theta # parameter theta_s of the sector copula cops
    switch(cops@copula@name, # sector copula family
           "AMH"={
               stop("a.coeff: AMH's family is currently not supported")
               ## Argh... the following code fails because copula:::coeffG does
               ## not support alpha > 1
               th0s <- (ths-th0)/(1-ths)
               n <- length(t)
               coeff <- sapply(t, function(t.) copula:::coeffG(ds,
                                                               1/(1-th0s*exp(-t.))))
               ## => coeff is an d_s x n matrix (rows = k., cols = t)
               t(coeff * (-1)^(ds-k.)) # n x d_s matrix
           },
           "Clayton"={
               as <- th0/ths
               n <- length(t)
               t. <- outer(1+t, as*k.-ds, FUN="^") # n x d_s matrix
               t. * rep((-1)^(ds-k.)*copula:::coeffG(ds, as), each=n)
           },
           "Frank"={
               stop("a.coeff: Frank's family is currently not supported")
           },
           "Gumbel"={
               as <- th0/ths
               n <- length(t)
               t. <- outer(t, as*k.-ds, FUN="^") # n x d_s matrix
               t. * rep((-1)^(ds-k.)*copula:::coeffG(ds, as), each=n)
           },
           "Joe"={
               stop("a.coeff: Joe's family is currently not supported")
           },
           stop("Family not supported"))
}

##' Coefficients b_k(t), k=d_0,..,d for all t's (sample size n) where
##'
##'     b_k(t) = \sum_{j in Q_{d.,k}^{d_0} \prod_{s=1}^{d_0} a_{s, d_s j_s}(t_s(u_s))
##'
##' where t_s is the sum over iPsi of sector s,
##'       a_{s, d_s j_s}(t) is a function specific to each family (and = 1 for
##'       non-sectorial parts),
##'   and Q_{d.,k}^{d_0} is the set of all vectors j=(j_1,..,j_d_0) such that the sum
##'       of the components j_1,..,j_d_0 equals k and j_s <= d_s for all s=1,..,d_0.
##' @title Coefficients b_k(t), k=d_0,..,d
##' @param t n x S matrix where the sth column contains the sum of psi_s^{-1}'s
##'        of the corresponding (proper) sector s
##' @param cop nacopula object of nesting depth 2 with S sectors and possibly
##'        non-sectorial parts (d_0 = dim(C0))
##' @return n x (d-d_0+1) matrix where the kth column contains b_k(t), k=d_0,..,d
##' @author Marius Hofert
b.coeff <- function(t, cop)
{
    stopifnot(is.function(blockparts <- partitions::blockparts),
              ## requireNamespace("partitions"),
              t >= 0, is(cop, "outer_nacopula"), nesdepth(cop) == 2,
              (S <- ncol(t)) == (lchild <- length(cop@childCops)))
    if(!is.matrix(t)) t <- rbind(t) # convert t to a matrix (if necessary)

    ## get copula dimensions
    C <- cop
    d <- dim(C) # dimension of C
    d0 <- length(C@comp) + lchild # dimension of C0
    d. <- rep(1, d0) # vector of dimension d0 giving the dimensions of each sector (including the degenerate ones in the end with dimension 1 each)
    d.[1:S] <- unlist(lapply(C@childCops, dim)) # build in the dimensions of each proper child cop (in the beginning)
    stopifnot(sum(d.)==d)

    ## build the matrices Q for each k=d0,..,d (according to an email
    ## conversation with Robin Hankin = maintainer("partitions") from 2012-01-21)
    Q <- lapply(d0:d, function(k) blockparts(d.-rep(1L, d0), k-d0)+1L)
    ## => the kth element of this list Q now contains an d_0 x ? matrix where ? is
    ##    the number of elements of Q_{d.,k}^{d_0}. Each column thus contains one
    ##    vector j in Q_{d.,k}^{d_0}

    ## build the list A of length d0 containing n x d_s matrices
    A <- lapply(1:d0, function(s){
        if(d.[s]==1) { # non-sectorial components
            matrix(1, n,1L)
        } else { # proper child copulas
            a.coeff(t[,s], cops=C@childCops[[s]], th0=C@copula@theta)
        }
    })
    ## => the sth (in 1:d0) element of this list contains a n x d_s matrix with
    ##    the values a_{s, d_s k}, k=1,..,d_s.

    ## build b_k, k=d0,..,d
    n <- nrow(t)
    b <- matrix(NA, nrow=n, ncol=d-d0+1)
    for(k in d0:d){ # compute b_k=b[,k-d0+1] (first arg = t)
        ind <- k-d0+1
        Q.k <- Q[[ind]] # pick out Q for this fixed k
        ## walk over the rows of Q.k (they contain all j[s] for s fixed)
        factor.t.s.js <- array(, dim=c(n, d0, ncol(Q.k))) # contains the factors
        for(s in 1:d0){
            j.s <- Q.k[s, ] # sth row of Q.k
            factor.t.s.js[, s, ] <- A[[s]][, j.s] # n x ncol(Q.k) matrix
        }
        prods <- apply(factor.t.s.js, MARGIN=c(1,3), FUN=prod) # n x ncol(Q.k) matrix
        b[,ind] <- rowSums(prods)
    }
    b # return
}


### likelihood #################################################################

##' Log-likelihood of a two-level nested Archimedean copula
##'
##' @title Log-likelihood of a two-level nested Archimedean copula
##' @param cop two-level nested Archimedean copula ("outer_nacopula")
##' @param u matrix of realizations/observations u
##' @return -log-likelihood
##' @author Marius Hofert
nacLL <- function(cop, u)
{
    if(!is.matrix(u)) u <- rbind(u)
    stopifnot(is(cop, "outer_nacopula"), nesdepth(cop) <= 2,
              (d <- ncol(u)) == dim(cop))
    C <- cop
    S <- length(C@childCops) # number of proper child copulas (<= d0)
    if(nesdepth(C) == 1 || S == 0) stop("nacLL: Have non-nested Archimedean copula;
 use the more efficient algorithms implemented for this special case.")
    ## FIXME?: rather return +Inf than an error ?
    if(C@copula@theta > min(unlist(lapply(C@childCops,
                                          function(cc) cc@copula@theta))))
	stop("Nested archimedean copula parameters do not fulfill nesting condition")
    ## setup
    th0 <- C@copula@theta
    lcomp <- length(C@comp) # if > 0 there is a non-sectorial part
    d0 <- S+lcomp # dim(C_0)
    n <- nrow(u)
    liPsiD1. <- matrix(NA, nrow=n, ncol=if(lcomp>0) S+1 else S) # n x S(+1) matrix of log(absdiPsi()); the "+1" comes from the non-sectorial part (if there is one)
    eta0. <- matrix(NA, nrow=n, ncol=d0) # n x d_0 matrix containing C_s(u_s)'s and u_0's (if there is a non-sectorial part)

    ## walk over the non-sectorial part (if available)
    if(lcomp > 0) {
        u0 <- u[,C@comp, drop=FALSE]
        liPsiD1.[,S+1] <- rowSums(C@copula@absdiPsi(u0, theta=th0, log=TRUE))
        eta0.[,(S+1):d0] <- u0 # just set to the u's of the non-sectorial part
    }

    ## now walk over all sectors
    t. <- matrix(NA, nrow=n, ncol=S) # n x S matrix of t's
    for(s in 1:S) {
        Cs <- C@childCops[[s]] # sector s copula list
        us <- u[, Cs@comp] # col-indices of u that belong to sector s
        c.s <- Cs@copula
        t.[,s] <- rowSums(c.s@iPsi(us, theta=c.s@theta))
        liPsiD1.[,s] <- rowSums(c.s@absdiPsi(us, theta=c.s@theta, log=TRUE))
        eta0.[,s] <- c.s@psi(t.[,s], theta=c.s@theta) # C_s(u_s)
    }

    ## finish computations
    liPsiD1sum <- rowSums(liPsiD1.) # sum(log(absdiPsi)); vector of length n
    eta0 <- rowSums(C@copula@iPsi(eta0., theta=th0)) # eta0; vector of length n; = C@copula@iPsi(pnacopula(C, u), theta=th0)

    ## compute b_k's, k=d0,..,d
    ## note: it suffices to give b.coeff only the sectorial t's
    b.mat <- b.coeff(t., cop=C) # n x (d-d_0+1) matrix

    ## compute x_k, k=d0,..,d, a n x (d-d_0+1) matrix
    x <- sapply(d0:d, function(k){
        log((-1)^(d-k) * b.mat[,k-d0+1]) +
            C@copula@absdPsi(eta0, theta=th0, degree=k, log=TRUE)
    })

    ## return
    sum(copula:::lsum(t(x)) + liPsiD1sum) # sum over all t's
}



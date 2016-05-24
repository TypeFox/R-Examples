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


### Empirical estimator of psi by Genest, Neslehova, Ziegel (2011) #############

##' @title Compute w_1,..,w_m and p_1,..,p_m as in Algorithm 1 in Genest, Neslehova, Ziegel (2011)
##' @param x data matrix (not necessarily in the unit hypercube)
##' @return (m, 2)-matrix containing the w's (1st col) and p's (2nd col)
##' @author Marius Hofert
w.p <- function(x)
{
    stopifnot((n <- nrow(x)) >= 1, (d <- ncol(x)) >= 1)
    W <- vapply(seq_len(n), function(i) sum( colSums(t(x)<x[i,])==d ) / (n+1), NA_real_)
    ptab <- table(W, dnn=NULL, deparse.level=0)/n # dimnames contain w_1,..,w_m in increasing order; entries p_m,..,p_1
    w <- as.numeric(dimnames(ptab)[[1]]) # w_1,..,w_m in increasing order; w[1] > 0, w[m] < 1
    p. <- as.numeric(ptab) # p_m,..,p_1
    cbind(w=w, p=rev(p.)) # dim = c(m,2)
}

##' @title Compute r_1,..,r_m, see page 233 in Genest, Neslehova, Ziegel (2011)
##' @param wp (m, 2)-matrix as returned by w.p()
##' @param d dimension
##' @param interval uniroot() initial interval
##' @param r.max (arbitrary) scaling for r_m
##' @param ... additional arguments passed to uniroot()
##' @return m-vector containing r_1,..,r_m=1
##' @author Marius Hofert
##' @note Apart from r_m=1, the r's can be quite small
r.solve <- function(wp, d, interval=c(0,1), r.max=1, ...)
{
    stopifnot(ncol(wp)==2, (m <- nrow(wp)) >= 1, diff(wp[1,]) > 0, d > 1)
    w <- wp[,1]
    p <- wp[,2]

    ## m >= 1
    if(m==1) return(1)

    ## m >= 2
    s <- numeric(m-1)
    s[m-1] <- 1-(w[2]/p[m])^(1/(d-1)) # well-defined (m >= 2)
    if(m==2) return(c(s[m-1], 1))

    ## m >= 3
    for(i in 3:m) { # entered at least once (at least one 'i')
        ## goal: compute s[m-i+1]
        js <- (m-i+2):m # js[1] ranges in m-1 to *2*
        ## a = (a_k) for k = 1,..,i-1 (k = j - (m-i+1))
        ## a = \prod_{l=m-i+2}^{j-1} s_l for j in {m-i+2,..,m}
        a <- c(1, # j = m-i+2 => empty product; deal with this separately
               vapply(js[-1], # (m-i+3) : m
                      function(j) prod(s[(m-i+2):(j-1)]), NA_real_)) # has length at least *2*
        f <- function(s) sum((1-s*a)^(d-1) * p[js]) - w[i] # maybe use log-scale here
        s[m-i+1] <- uniroot(f, interval=interval, ...)$root
    }

    ## compute r_1,.., r_m
    rev(r.max*c(1, cumprod(rev(s)))) # r_m = r.max, r_i = s_i * r_{i+1}, i in {m-1,..,1}; TODO
}

##' @title Compute Pseudo-observations R_1, .., R_n from r_1,..,r_m and p_1,..,p_m
##' @param r r_1,..,r_m as in Genest, Neslehova, Ziegel (2011)
##' @param p p_1,..,p_m as in Genest, Neslehova, Ziegel (2011)
##' @param n sample size
##' @param sample logical indicating whether the result is shuffled
##' @return n-vector of pseudo-observations from the radial parts
##' @author Marius Hofert
R.pobs <- function(r, p, n, sample=FALSE) {
    stopifnot((m <- length(r)) == length(p), n >= 1)
    res <- rep(r, n*p)
    if(sample) sample(res) else res
}

##' @title Non-parametric Generator Estimator
##' @param x evaluation points
##' @param r r_1,..,r_m as in Genest, Neslehova, Ziegel (2011)
##' @param p p_1,..,p_m as in Genest, Neslehova, Ziegel (2011)
##' @param d dimension
##' @return psi_{n,d}(x)
##' @author Marius Hofert
psi.n <- function(x, r, p, d) {
    stopifnot((m <- length(r)) == length(p), d > 1)
    vapply(x, function(x.) if(any(ii <- r > x.)) sum((1-x./r[ii])^(d-1) * p[ii]) else 0, NA_real_)
}

##' @title Inverse of the Non-parametric Generator Estimator
##' @param u evaluation points (n-vector in [0,1])
##' @param r r_1,..,r_m as in Genest, Neslehova, Ziegel (2011)
##' @param p p_1,..,p_m as in Genest, Neslehova, Ziegel (2011)
##' @param d dimension
##' @param ... additional arguments passed to uniroot()
##' @return psi_{n,d}^{-1}(u)
##' @author Marius Hofert
##' Note: 1) psi.n(r_{m-i}) = w_{i+1}, i = 1, .., m-1
##'          x-axis: r_1, r_2, .., r_m [psi.n(r_m)=0]
##'          y-axis: w_1=0, ..., w_{m-2}, w_{m-1}, w_m < 1; psi.n(0)=1, w_m
##'       2) \hat{U}_i are also j/(n+1) as w_k's => if iPsi.n(\hat{U}_i) = w_k for some k
##'          one could return the corresponding r
iPsi.n <- function(u, r, p, d, interval=c(0, max.r), ...) {
    stopifnot((m <- length(r)) == length(p), d > 1, 0 <= u, u <= 1, (n <- length(u)) >= 1)
    max.r <- max(r)
    res <- rep(max.r, n) # => psi_{n,d}^{-1}(u) = max.r for all u <= psi_{n,d}(max.r)
    f <- function(x, u) psi.n(x, r=r, p=p, d=d) - u
    ii <- u > psi.n(max.r, r=r, p=p, d=d)
    if(any(ii)) {
        res[ii] <- vapply(u[ii], function(u.)
                          uniroot(f, interval=interval, u=u., ...)$root, NA_real_)
    }
    res
}

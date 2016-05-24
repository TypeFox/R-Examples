###############################################################################
## Utilities for use with TLNise
## Copyright (C) 2005, Roger D. Peng <rpeng@jhsph.edu>
##     
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
###############################################################################

## B0 is a p x p x N array
transformB0 <- function(B0, V) {
        V0 <- apply(V, c(1, 2), mean)
        eig.V0 <- eigen(V0)
        rtV0 <- eig.V0$vectors %*% sqrt(diag(eig.V0$values, nrow = length(eig.V0$values))) %*% t(eig.V0$vectors)
        lapply(1:(dim(B0)[3]), function(i) {
                rtV0 %*% (solve(B0[, , i]) - diag(1, NROW(B0[, , i]))) %*% rtV0
        })
}

## Y is an n x p matrix of city-specific coefficients
## V is a p x p x n array of city-specific covariance matrices
## W is a n x k matrix of k second stage regression variables (without
## the intercept)

## If you have single coefficients for each city then p == 1; Y and V
## are both vectors of length n.  If there are no second stage
## variables, then an intercept is assumed.

sampleGamma <- function(tlnise.out, V, Y, W = NULL) {
        if(is.null(dim(V)))
                V <- array(V, c(1, 1, NROW(V)))
        if(is.null(dim(Y)))
                Y <- matrix(Y, nrow = NROW(Y), ncol = 1)
        A <- transformB0(tlnise.out$B0, V)
        stopifnot(is.list(A))
        n <- dim(V)[3]  ## Number of cities
        wt <- exp(tlnise.out$lr - max(tlnise.out$lr))

        ## These formulas are taken directly from the Everson & Morris
        ## JRSS-B paper (2000)
        gammastar <- lapply(A, function(Aj) {
                Z <- lapply(1:n, function(k) {
                        if(!is.null(W))
                                w <- kronecker(diag(NCOL(Aj)), W[k, ])
                        else
                                w <- diag(NCOL(Aj))
                        Z1 <- w %*% solve(V[, , k] + Aj) %*% t(w)
                        Z2 <- drop(w %*% solve(V[,,k] + Aj) %*% Y[k, ])
                        list(Z1 = Z1, Z2 = Z2)
                })
                Z1 <- lapply(Z, "[[", "Z1")
                r <- array(unlist(Z1), c(dim(Z1[[1]]), n))
                Dstar <- solve(apply(r, c(1,2), sum))
                Z2 <- lapply(Z, "[[", "Z2")
                r <- matrix(unlist(Z2), nrow = ncol(Y) * NCOL(W))
                g <- drop(Dstar %*% rowSums(r))

                list(g = g, Dstar = Dstar)
        })
        g <- unlist(lapply(gammastar, "[[", "g"))
        gs <- matrix(g, byrow = TRUE, ncol = NCOL(Y) * NCOL(W))
        Dstar <- lapply(gammastar, "[[", "Dstar")
        gamma <- apply(gs, 2, weighted.mean, w = wt)
        
        list(gamma = gamma, gammastar = gs, A = A, wt = wt, Dstar = Dstar)
}

drawSample0 <- function(sg.out, n = 1000) {
        gammastar <- sg.out$gammastar
        Dstar <- sg.out$Dstar
        N <- NROW(gammastar)
        w <- sg.out$wt
        draw <- lapply(seq(n), function(i) {
                use <- sample(1:N, 1, prob = w)
                mvrnorm(1, gammastar[use, ], Dstar[[use]])
        })
        do.call("rbind", draw)
}

drawSample <- function(tlnise.out, n = 1000, Y, V, W = NULL) {
        if(!require(MASS))
                stop("'MASS' package required for 'drawSample'")
        sg.out <- sampleGamma(tlnise.out, V, Y, W)
        gammastar <- sg.out$gammastar
        Dstar <- sg.out$Dstar
        N <- NROW(gammastar)
        w <- sg.out$wt
        draw <- lapply(seq(n), function(i) {
                use <- sample(1:N, 1, prob = w)
                mvrnorm(1, gammastar[use, ], Dstar[[use]])
        })
        do.call("rbind", draw)
}

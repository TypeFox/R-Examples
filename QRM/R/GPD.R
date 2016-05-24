## Copyright (C) 2013 Marius Hofert, Bernhard Pfaff
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


## GPD function (fine for all q in IR, xi in IR, beta > 0)
## vectorized in q
pGPD <- function(q, xi, beta = 1)
{
    stopifnot(beta > 0)
    q. <- q / beta
    res <- if(xi != 0) 1 -(1+xi*q.)^(-1/xi) else pexp(q.)
    res[q. < 0] <- 0
    if(xi < 0) res[q. > -1 / xi] <- 1
    res
}

## quantile of the GPD function (fine for all p in [0,1], xi in IR, beta > 0)
## vectorized in p
qGPD <- function(p, xi, beta = 1){
    stopifnot(beta > 0)
    p. <- pmax(pmin(p, 1), 0) # or not allow by stopifnot(0 <= p, p <= 1)? [rest still fine]
    if(xi == 0) qexp(p., rate=1/beta) else (beta/xi)*((1-p.)^(-xi)-1) # p-quantile for xi != 0
}

rGPD <- function(n, xi, beta = 1) qGPD(runif(n), xi, beta)

dGPD <- function(x, xi, beta = 1, log = FALSE){
    xb <- x / beta
    res <- if(xi == 0){
        log(dexp(xb)) - log(beta)
    } else {
        ind <- if(xi < 0) xb > 0 & (xb < 1 / abs(xi)) else xb > 0
        r <- rep(-Inf, length(x))
        r[ind] <- (-1 / xi-1) * log(1 + xi * xb[ind]) - log(beta)
        r
    }
    if(log) res else exp(res)
}

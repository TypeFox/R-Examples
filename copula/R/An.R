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


##' Nonparametric estimators of the Pickands dependence function
##' Bivariate versions
##'
##' @title Rank-based versions of the bivariate Pickands and CFG estimators
##' @param x data
##' @param w points were to estimate A
##' @param estimator CFG or Pickands
##' @param corrected if TRUE, endpoint corrections applied
##' @return values of estimated A at w
##' @author Ivan Kojadinovic
An.biv <- function(x, w, estimator = c("CFG", "Pickands"), corrected = TRUE) {
    n <- nrow(x)
    m <- length(w)

    ## make pseudo-observations
    u <- pobs(x)
    mlu <- -log(u)

    switch(match.arg(estimator),
           "CFG" =
           .C(biv_ACFG,
              as.integer(n),
              mlu[,1],
              mlu[,2],
              as.double(w),
              as.integer(m),
              as.integer(corrected),
              A = double(m))$A,
           "Pickands" =
           .C(biv_AP,
              as.integer(n),
              mlu[,1],
              mlu[,2],
              as.double(w),
              as.integer(m),
              as.integer(corrected),
              A = double(m))$A,
           stop("invalid 'estimator' : ", estimator))
}

Anfun <- function(x, w, estimator = c("CFG", "Pickands"), corrected = TRUE) {
    .Deprecated("An.biv")
    An.biv(x, w, estimator, corrected=corrected)
}


##' Nonparametric estimators of the Pickands dependence function
##' Mulivariate P, CFG, HT corrected versions
##'
##' @title Rank-based versions of the multivariate Pickands and CFG estimators
##' @param x data
##' @param w points were to estimate A
##' @return values of estimated A at w
##' @author Ivan Kojadinovic
An <- function(x, w) {

    d <- ncol(x)

    if (d < 2)
        stop("The data should be at least of dimension 2")
    if (ncol(w) != d)
        stop("The matrices 'x' and 'w' should have the same number of columns")

    n <- nrow(x)
    m <- nrow(w)

    .C(mult_A,
       as.double(pobs(x)),
       as.integer(n),
       as.integer(d),
       as.double(w),
       as.integer(m),
       P = double(m),
       CFG = double(m),
       HT = double(m))[c("P","CFG","HT")]
}

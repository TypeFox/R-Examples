## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


### Wrappers and auxiliary functions for dealing with elliptical (Gauss, t_nu)
### and Archimedean copulas

##' @title Copula class for the given copula object
##' @param cop copula object
##' @return "ellipCopula" or "outer_nacopula" depending on the given copula object
##' @author Marius Hofert
copClass <- function(cop)
{
    cls <- class(cop)
    if(is(cop, "copula") && (cls=="normalCopula" || cls=="tCopula")) "ellipCopula" # note: there could be other "copula" objects which are not elliptical
    else if(cls=="outer_nacopula") "outer_nacopula" # can be Archimedean or nested Archimedean
    else stop("not yet supported copula object")
}

##' @title Copula family for the given copula object
##' @param cop copula object (either elliptical or (nested) Archimedean)
##' @return family string
##' @author Marius Hofert
copFamily <- function(cop)
{
    cls <- class(cop)
    if(is(cop, "copula")){
        if(cls=="normalCopula") "normal"
        else if(cls=="tCopula") "t"
        else stop("unsupported copula family")
    } else if(cls=="outer_nacopula"){
        cop@copula@name # could be nested or not
    } else stop("not yet supported copula object")
}

##' @title Copula family for the given copula object
##' @param cop copula object (either elliptical or (nested) Archimedean)
##' @return family string
##' @author Marius Hofert
copFamilyClass <- function(family)
{
    if(family == "normal" || family == "t")
	"ellipCopula"
    else if(family %in% .ac.longNames ||
	    family %in% paste0("opower:", .ac.longNames))
	"outer_nacopula" # note: opower not really supported yet
    else stop("family ", family, " not yet supported")
}

##' @title Construct a symmetric matrix with 1s on the diagonal from the given
##'        parameter vector
##' @param param parameter vector
##' @param d number of columns (or rows) of the output matrix
##' @return a symmetric matrix with 1s on the diagonal and the values of param
##'         filled column-wise below the diagonal (= row-wise above the diagonal)
##' @author Marius Hofert
p2P <- function(param, d)
{
    P <- diag(1, nrow=d)
    P[lower.tri(P)] <- param
    P <- P+t(P)
    diag(P) <- rep.int(1, d)
    P
}

##' @title Extract the vector of column-wise below-diagonal entries from a matrix
##' @param P matrix (typically a symmetric matrix as used for elliptical copulas)
##' @return the vector of column-wise below-diagonal entries of P (they are equal
##'         to the row-wise above-diagonal entries in case of a symmetric matrix)
##' @author Marius Hofert
##' Note: This is used "by foot" at several points in the package.
P2p <- function(P) P[lower.tri(P)]

##' @title Construct matrix Sigma from a given elliptical copula
##' @param copula copula
##' @return (d, d) matrix Sigma containing the parameter vector rho
##' @author Marius Hofert
getSigma <- function(copula)
{
    stopifnot(is(copula, "ellipCopula"))
    d <- copula@dimension
    rho <- copula@getRho(copula)
    switch(copula@dispstr,
	   "ex" = {
	       Sigma <- matrix(rho[1], nrow=d, ncol=d)
	       diag(Sigma) <- rep(1, d)
	   },
	   "ar1" = {
	       Sigma <- rho^abs(outer(1:d, 1:d, FUN="-"))
	   },
	   "un" = {
	       Sigma <- p2P(rho, d)
	   },
	   "toep" = {
	       rho <- c(rho, 1)
	       ind <- outer(1:d, 1:d, FUN=function(i, j) abs(i-j))
	       diag(ind) <- length(rho)
	       Sigma <- matrix(rho[ind], nrow=d, ncol=d)
	   },
	   stop("invalid 'dispstr'"))
    Sigma
}

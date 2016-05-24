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


## constructor #################################################################

fgmCopula <- function(param = NA_real_, dim = 2L) {
    if (!is.numeric(dim) || (dim <- as.integer(dim)) < 2)
        stop("dim should be an integer of at least 2")
    npar <- 2^dim - dim - 1
    if(is.na(param)) param <- rep(param, length.out = npar)
    else if (!is.numeric(param) && length(param) != npar)
        stop("wrong parameters")

    ## power set of {1,...,dim} in integer notation
    subsets <-  .C(k_power_set,
                   as.integer(dim),
                   as.integer(dim),
                   subsets = integer(2^dim))$subsets
    ## power set in character vector: {}, {1}, {2}, ..., {1,2}, ..., {1,...,dim}
    subsets.char <-  .C(k_power_set_char,
                        as.integer(dim),
                        as.integer(2^dim),
                        as.integer(subsets),
                        sc = character(2^dim))$sc

    ## expression of the cdf
    cdfExpr <- function(n,sc) {

        expr1 <- "u1"
        for (i in 2:n)
            expr1 <- paste0(expr1, " * u", i)

        expr2 <- "1"
        for (i in (dim + 2):2^dim)
        {
            expr3 <- paste0("alpha",i)
            sub <- substr(sc[i],2,nchar(sc[i])-1)
            for (j in eval(strsplit(sub,",")[[1]]))
                expr3 <- paste0(expr3, " * (1 - u", j, ")")
            expr2 <- paste(expr2,"+",expr3)
        }

        expr <- paste(expr1," * (", expr2, ")")
        parse(text = expr)
    }

    ## expression of the pdf
    pdfExpr <- function(n,sc) {
        expr2 <- "1"
        for (i in (dim + 2):2^dim)
        {
            expr3 <- paste0("alpha",i)
            sub <- substr(sc[i],2,nchar(sc[i])-1)
            for (j in eval(strsplit(sub,",")[[1]]))
                expr3 <- paste0(expr3, " * (1 - 2 * u", j, ")")
            expr2 <- paste(expr2,"+",expr3)
        }
        parse(text = expr2)
    }

    cdf <- cdfExpr(dim,subsets.char)
    pdf <- pdfExpr(dim,subsets.char)

    ## create new object
    new("fgmCopula",
        dimension = dim,
        parameters = param,
        exprdist = c(cdf = cdf, pdf = pdf),
        param.names = paste0("param",subsets.char[(dim+2):2^dim]),
        param.lowbnd = rep(-1, 2^dim - dim - 1),
        param.upbnd = rep(1, 2^dim - dim - 1),
        fullname = "Farlie-Gumbel-Morgenstern copula family")
}


### random number generation ###################################################

##' @title Random number generation for a FGM copula
##' @param n sample size
##' @param copula object of type 'fgmCopula'
##' @param method method (R code based on Remillard (2013, p. 310); C code from
##'        C code from Ivan Kojadinovic)
##' @return (n, d) matrix of copula samples
##' @author Marius Hofert
rfgmCopula <- function(n, copula, method=c("C", "R"))
{
    stopifnot(is(object = copula, "fgmCopula"))
    d <- copula@dimension
    alpha <- copula@parameters
    stopifnot(d >= 2, -1 <= alpha, alpha <= 1)
    method <- match.arg(method)

    switch(method,
           "R" = {
               ## note: this is *only* for the case where S = {1,...,d}
               ##       => the FGM has only one parameter in this case
               ##       see Jaworski, Durante, Haerdle, Rychlik (2009, p. 19)
               warning("random number generation only for homogeneous (one-parameter) case")
               stopifnot(n >= 1, d >= 2, -1 <= alpha, alpha <= 1)
               U <- matrix(runif(n*d), nrow=n, ncol=d)
               B <- alpha * apply(1-2*U[,-d, drop=FALSE], 1, prod)
               C <- sqrt( (1 + B)^2 - 4 * B * U[,d])
               U[,d] <- 2 * U[,d] / (1 + B + C)
               U
           },
           "C" = {
               if (d > 2)
                   warning("random generation in C needs to be properly tested")
               matrix(.C(rfgm,
                         as.integer(d),
                         c(rep.int(0., d+1), alpha),
                         as.integer(n),
                         out = double(n * d))$out, n, d, byrow=TRUE)
           },
           stop("wrong 'method': ", method))
}


### cdf of the copula ##########################################################

pfgmCopula <- function(u, copula) {
    ## if(..) fails when u has NAs ; treatment now via generic pCopula()
    ## if (any(u < 0) || any(u > 1))
    ##     stop("u values should lie between 0 and 1")
    dim <- copula@dimension
    param <- copula@parameters
    cdf <- copula@exprdist$cdf
    for (i in 1:dim) assign(paste0("u", i), u[,i])
    for (i in (dim + 2):2^dim) assign(paste0("alpha", i), param[i - dim - 1])
    eval(cdf)
}

## pdf of the copula ###########################################################

dfgmCopula <- function(u, copula, log=FALSE, ...) {
    ## if(..) fails when u has NAs ; treatment now via generic dCopula()
    ## if (any(u < 0) || any(u > 1))
    ##     stop("u values should lie between 0 and 1")
    dim <- copula@dimension
    param <- copula@parameters
    pdf <- copula@exprdist$pdf
    for (i in 1:dim) assign(paste0("u", i), u[,i])
    for (i in (dim + 2):2^dim) assign(paste0("alpha", i), param[i - dim - 1])
    ## FIXME: improve log-case
    if(log) log(eval(pdf)) else eval(pdf)
}

## Kendall's tau

tauFgmCopula <- function(copula) {
    alpha <- copula@parameters[1]
    2 * alpha / 9
}

## Spearman's rho

rhoFgmCopula <- function(copula) {
    alpha <- copula@parameters[1]
    1 * alpha / 3
}

## calibration via tau

iTauFgmCopula <- function(copula, tau) {
  if (any(tau < -2/9 | tau > 2/9))
    warning("tau is out of the range [-2/9, 2/9]")
  pmax(-1, pmin(1, 9 * tau / 2))
}

## calibration via rho

iRhoFgmCopula <- function(copula, rho) {
  if (any(rho < -1/3 | rho > 1/3))
    warning("rho is out of the range [-1/3, 1/3]")
  pmax(-1, pmin(1, 3 * rho))
}


################################################################################

setMethod("rCopula", signature("numeric", "fgmCopula"), rfgmCopula)

setMethod("pCopula", signature("numeric", "fgmCopula"),pfgmCopula)
setMethod("pCopula", signature("matrix", "fgmCopula"), pfgmCopula)
setMethod("dCopula", signature("numeric", "fgmCopula"),dfgmCopula)
setMethod("dCopula", signature("matrix", "fgmCopula"), dfgmCopula)

setMethod("tau", signature("fgmCopula"), tauFgmCopula)
setMethod("rho", signature("fgmCopula"), rhoFgmCopula)
setMethod("iTau", signature("fgmCopula"), iTauFgmCopula)
setMethod("iRho", signature("fgmCopula"), iRhoFgmCopula)

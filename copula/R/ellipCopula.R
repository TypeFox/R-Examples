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

##' Length of parameter vector for an 'ellipCopula' (tCopula or normCopula)
##'
##' @title Parameter length of an 'ellipCopula'
##' @param dim dimension of copula
##' @param dispstr dispersion string
##' @param df.fixed [for tCopula():] logical indicating if 'df' is fixed or a parameter
##' @return integer
##' @author Martin Maechler
## NOTE: This is related to  validRho() in ./Classes.R
npar.ellip <- function(dim, dispstr, df.fixed = TRUE) {
    (!df.fixed) + # add 1 if 'df' is not fixed
    switch(dispstr, ## also checking for correct 'dispstr'
	   "ar1" =, "ex" = 1 ,
	   "un" = dim * (dim - 1) / 2,
	   "toep" = dim - 1,
	   ## otherwise
	   return("'dispstr' not supported (yet)"))
}

ellipCopula <- function(family, param = NA_real_, dim = 2L, dispstr = "ex",
                        df = 4, ...)
{
  familiesImplemented <- c("normal", "t")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", paste(familiesImplemented, collapse=", ")))
  switch(fam,
	 normalCopula(param, dim = dim, dispstr = dispstr),
	      tCopula(param, dim = dim, dispstr = dispstr, df = df, ...)
         )
}

iTauEllipCopula <- function(copula, tau) sin((tau * pi) / 2)

iRhoEllipCopula <- function(copula, rho) sin(pi * rho / 6) * 2

dTauEllipCopula <- function(copula)  {
  2 / (pi * sqrt(1 - copula@getRho(copula)^2))
}

dTauFunEllipCopula <- function(copula)  {
  function(x) 2 / (pi * sqrt(1 - x^2))
}

dRhoEllipCopula <- function(copula) {
  6 / (pi * sqrt(4 - copula@getRho(copula)^2))
}

dRhoFunEllipCopula <- function(copula) {
  function(x) 6 / (pi * sqrt(4 - x^2))
}

setMethod("iTau", signature("ellipCopula"), iTauEllipCopula)
setMethod("iRho", signature("ellipCopula"), iRhoEllipCopula)

setMethod("dTau", signature("ellipCopula"), dTauEllipCopula)
setMethod("dRho", signature("ellipCopula"), dRhoEllipCopula)

setMethod("dTauFun", signature("ellipCopula"), dTauFunEllipCopula)
setMethod("dRhoFun", signature("ellipCopula"), dRhoFunEllipCopula)

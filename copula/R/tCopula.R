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


## DONE: allow "param = NA" --
##  2) npar.ellip(dim, dispstr, df.fixed) |-->  "dimension" (length) of param
##  --> TODO: use it ---> define a generic  nparam(.) {with methods for all copula}
### TODO:  {also for "normalCopula"}
##  3) validity should check  pos.definiteness for "un"structured (maybe "toeplitz"
##
tCopula <- function(param = NA_real_, dim = 2L, dispstr = "ex",
		    df = 4, df.fixed = FALSE)
{
    dim <- as.integer(dim)
    stopifnot((pdim <- length(param)) >= 1, is.numeric(param))
    if(pdim == 1 && is.na(param)) ## extend it (rho only!)
	pdim <- length(param <- rep(param, length.out = npar.ellip(dim, dispstr)))
    parameters <- param
    param.names <- paste("rho", seq_len(pdim), sep=".")
    param.lowbnd <- rep.int(-1, pdim)
    param.upbnd	 <- rep.int( 1, pdim)
    if (!df.fixed) { ## df is another parameter __at end__
	parameters <- c(parameters, df)
	param.names <- c(param.names, "df")
	param.lowbnd <- c(param.lowbnd, 1e-6)
	param.upbnd <- c(param.upbnd, Inf)
    }

    new("tCopula",
	dispstr = dispstr,
	dimension = dim,
	parameters = parameters,
	df = df,
	df.fixed = df.fixed,
	param.names = param.names,
	param.lowbnd = param.lowbnd,
	param.upbnd = param.upbnd,
	fullname = paste("t copula family", if(df.fixed) paste("df fixed at", df)),
	getRho = function(obj) {
	    par <- obj@parameters
	    if (obj@df.fixed) par else par[-length(par)]
	}
	)
}

## used for tCopula *and* tevCopula :
as.df.fixed <- function(cop, classDef = getClass(class(cop))) {
    stopifnot( ## fast .slotNames()
	      any(names(classDef@slots) == "df.fixed"))
    if (!cop@df.fixed) { ## df is another parameter __at end__
	cop@df.fixed <- TRUE
	ip <- seq_len(length(cop@parameters) - 1L)
	cop@parameters	 <- cop@parameters  [ip]
	cop@param.names	 <- cop@param.names [ip]
	cop@param.lowbnd <- cop@param.lowbnd[ip]
	cop@param.upbnd	 <- cop@param.upbnd [ip]
    }
    cop
}

getdf <- function(object) {
  if (object@df.fixed) object@df
  else ## the last par. is 'df'
      object@parameters[length(object@parameters)]
}

rtCopula <- function(n, copula) {
  df <- getdf(copula)
  pt(rmvt(n, sigma = getSigma(copula), df = df), df = df)
}


ptCopula <- function(u, copula, ...) {
  dim <- copula@dimension
  i.lower <- rep.int(-Inf, dim)
  sigma <- getSigma(copula)
  df <- getdf(copula)
  if(!(df==Inf || df == as.integer(df)))
    stop("'df' is not integer (or Inf); therefore, pCopula() cannot be computed yet")
  ## more checks now  pCopula() *generic*
  apply(u, 1, function(x) if(any(is.na(x))) NA_real_ else
	pmvt(lower = i.lower, upper = qt(x, df = df), sigma = sigma, df = df, ...))
}

dtCopula <- function(u, copula, log = FALSE, ...) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  df <- getdf(copula)
  ## more checks now  dCopula() *generic*
  r <- numeric(nrow(u)) # i.e. 0  by default (i.e. "outside")
  ok <- u.in.01(u)
  x <- qt(u[ok, , drop=FALSE], df)
  ## work in log-scale [less over-/under-flow, then (maybe) transform]:
  r[ok] <- dmvt(x, delta = rep.int(0, dim), sigma = sigma, df = df, log = TRUE) -
      rowSums(dt(x, df = df, log=TRUE))
  if(log) r else exp(r)
}

showTCopula <- function(object) {
  print.copula(object)
  if (object@dimension > 2) cat("dispstr: ", object@dispstr, "\n")
  if (object@df.fixed) cat("df is fixed at", object@df, "\n")
  invisible(object)
}

tailIndexTCopula <- function(copula) {
### McNeil, Frey, Embrechts (2005), p.211
  df <- getdf(copula)
  rho <- copula@getRho(copula)
  upper <- lower <- 2 * pt(- sqrt((df + 1) * (1 - rho) / (1 + rho)), df=df + 1)
  c(upper=upper, lower=lower)
}

tauTCopula <- function(copula) {
  rho <- copula@getRho(copula)
  2 * asin(rho) /pi
}

rhoTCopula <- function(copula) {
  rho <- copula@getRho(copula)
  asin(rho / 2) * 6 / pi
}

setMethod("rCopula", signature("numeric", "tCopula"), rtCopula)

setMethod("pCopula", signature("matrix", "tCopula"), ptCopula)
setMethod("pCopula", signature("numeric", "tCopula"),ptCopula)
setMethod("dCopula", signature("matrix", "tCopula"), dtCopula)
setMethod("dCopula", signature("numeric", "tCopula"),dtCopula)



setMethod("show", signature("tCopula"), showTCopula)

setMethod("tau", signature("tCopula"), tauTCopula)
setMethod("rho", signature("tCopula"), rhoTCopula)
setMethod("tailIndex", signature("tCopula"), tailIndexTCopula)

setMethod("iTau", signature("tCopula"), iTauEllipCopula)
setMethod("iRho", signature("tCopula"), iRhoEllipCopula)

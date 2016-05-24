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


joeCopula <- function(param = NA_real_, dim = 2L,
		      use.indepC = c("message", "TRUE", "FALSE"))
{
  stopifnot(length(param) == 1)
  if(!is.na(param) && param == 1) {
      use.indepC <- match.arg(use.indepC)
      if(!identical(use.indepC, "FALSE")) {
	  if(identical(use.indepC, "message"))
	      message("parameter at boundary ==> returning indepCopula()")
	  return( indepCopula(dim=dim) )
      }
  }
  new("joeCopula",
      dimension = as.integer(dim),
      parameters = param[1],
      param.names = "param",
      param.lowbnd = 1, # 0.238733989880086 for tau >= -1 -- is NOT valid
      param.upbnd = Inf,
      fullname = "Joe copula family; Archimedean copula")
}

setMethod("rCopula", signature("numeric", "joeCopula"),
	  function (n, copula, ...)
	  racopula(n, copJoe, theta=copula@parameters, d = copula@dimension))

setMethod("pCopula", signature("numeric", "joeCopula"),
	  function(u, copula, ...) pacopula(u, copJoe, theta=copula@parameters))
setMethod("pCopula", signature("matrix", "joeCopula"),
	  function(u, copula, ...) .pacopula(u, copJoe, theta=copula@parameters))

setMethod("dCopula" , signature("matrix", "joeCopula"),
	  function (u, copula, log = FALSE, ...)
	  copJoe@dacopula(u, theta=copula@parameters, log=log, ...))
setMethod("dCopula", signature("numeric", "joeCopula"),
	  function (u, copula, log = FALSE, ...)
	  copJoe@dacopula(rbind(u, deparse.level=0L),
			  theta=copula@parameters, log=log, ...))

setMethod("psi", signature("joeCopula"),
	  function(copula, s) copJoe@psi(t=s, theta=copula@parameters))
setMethod("iPsi", signature("joeCopula"),
	  function(copula, u) copJoe@iPsi(u, theta=copula@parameters))
setMethod("diPsi", signature("joeCopula"),
	  function(copula, u, degree=1, log=FALSE, ...)
      {
	  s <- if(log || degree %% 2 == 0) 1. else -1.
	  s* copJoe@absdiPsi(u, theta=copula@parameters, degree=degree, log=log, ...)
      })

setMethod("tau", signature("joeCopula"),
          function(copula) tauJoe(theta=copula@parameters))
setMethod("tailIndex", signature("joeCopula"),
	  function(copula) c(lower=0,
			     upper=copJoe@lambdaU(theta=copula@parameters)))

setMethod("iTau", signature("joeCopula"),
	  function(copula, tau, tol = 1e-7) copJoe@iTau(tau, tol=tol))
                                        # now that tauJoe() is accurate

## "TODO"
## setMethod("rho", signature("joeCopula"), ... ? ...)
## setMethod("iRho", signature("joeCopula"), function(copula, rho) ...)

## "TODO"
## setMethod("dRho", signature("joeCopula"), ...)
## setMethod("dTau", signature("joeCopula"), ...)

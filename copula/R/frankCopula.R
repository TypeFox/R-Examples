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

### Parameter  'a', also called 'alpha' or 'theta'
### must be in (0, Inf)   for completely monotone psi ( <==> tau >= 0 )
### a = 0  <==> independence copula

iPsiFrank <- function(copula, u) .iPsiFrank(u, copula@parameters[1])
.iPsiFrank <- function(u, a) { ## FIXME: copFrank@iPsi() is clearly better for |a| << 1
  - log(expm1(- a * u) / expm1(-a))
}

psiFrank <- function(copula, s) .psiFrank(s, copula@parameters[1])
.psiFrank <- function(s, a) {
  ## -1/a * log(1 + exp(-s) * (exp(-a) - 1))
  stopifnot(length(a) == 1)
  if(a > 0)
    copFrank@psi(s, a) # using log1mexp(.)
  else if(a == 0)
    exp(-s)
  else if(a <= log(.Machine$double.eps))# -36.04
    ## exp(-a) -1  ~= exp(-a)
    -log1pexp(-(s+a))/a
  else ## -36.04 < a < 0 :
    -log1p(exp(-s) * expm1(-a))/a
}

## psiDerFrank <- function(copula, s, n) {
##   eval(psiDerFrank.expr[n + 1], list(s=s, alpha=copula@parameters[1]))
## }

frankCopula <- function(param = NA_real_, dim = 2L,
			use.indepC = c("message", "TRUE", "FALSE"))
{
  stopifnot(length(param) == 1)
  if((dim <- as.integer(dim)) > 2 && !is.na(param) && param < 0)
    stop("param can be negative only for dim = 2")
  if(!is.na(param) && param == 0) {
      use.indepC <- match.arg(use.indepC)
      if(!identical(use.indepC, "FALSE")) {
	  if(identical(use.indepC, "message"))
	      message("parameter at boundary ==> returning indepCopula()")
	  return( indepCopula(dim=dim) )
      }
  }

  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <-   "- log( (exp(- alpha * u1) - 1) / (exp(- alpha) - 1) )"
    for (i in 2:n) {
      cur <- paste0("- log( (exp(- alpha * u", i, ") - 1) / (exp(- alpha) - 1))")
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- paste("-1/alpha * log(1 + exp(-(", expr, ")) * (exp(-alpha) - 1))")
    parse(text = expr)
  }

  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste0("u", i))
    }
    val
  }

  cdf <- cdfExpr(dim)
  pdf <- if (dim <= 6) pdfExpr(cdf, dim) # else NULL
  new("frankCopula",
      dimension = dim,
      parameters = param,
      exprdist = c(cdf = cdf, pdf = pdf),
      param.names = "param",
      param.lowbnd = if(dim == 2) -Inf else 0,
      param.upbnd = Inf,
      fullname = "Frank copula family; Archimedean copula")
}

rfrankBivCopula <- function(n, alpha) {
  U <- runif(n); V <- runif(n)
  ## FIXME : |alpha| << 1  (including alpha == 0)
  ## to fix numerical rounding problems for alpha >35 but not for alpha < -35 :
  a <- -abs(alpha)
  ## reference: Joe (1997, p.147)
  V <- -1/a * log1p(-V * expm1(-a) / (exp(-a * U) * (V - 1) - V))
  cbind(U, if(alpha > 0) 1 - V else V,  deparse.level=0L)
}

rfrankCopula <- function(n, copula) {
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if (abs(alpha) < .Machine$double.eps ^ (1/3))
    return(rCopula(n, indepCopula(dim)))
##   if (abs(alpha) <= .Machine$double.eps^.9)
##     return (matrix(runif(n * dim), nrow = n))
  if (dim == 2) return (rfrankBivCopula(n, alpha))
  ## the frailty is a log series distribution with a = 1 - exp(-alpha)
  fr <- rlogseries(n, -expm1(-alpha))
  fr <- matrix(fr, nrow = n, ncol = dim)
  U <- matrix(runif(dim * n), nrow = n)
  psi(copula, - log(U) / fr)
}


pfrankCopula <- function(copula, u) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  cdf <- copula@exprdist$cdf
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (apply(u, 1, prod))
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  eval(cdf)
}

dfrankCopula <- function(u, copula, log=FALSE, ...) {
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  pdf <- copula@exprdist$pdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  alpha <- copula@parameters[1]
  if(log) stop("'log=TRUE' not yet implemented")
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep(1, nrow(u)))
  val <- eval(pdf)
#  val[apply(u, 1, function(v) any(v <= 0))] <- 0
#  val[apply(u, 1, function(v) any(v >= 1))] <- 0
  val
}

## dfrankCopula.expr <- function(copula, u) {
##   if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
##   s <- apply(iPsiFrank(copula, u), 1, sum)
##   pdf <- psiDerFrank(copula, s, copula@dimension) *
##     apply(iPsiDerFrank(copula, u, 1), 1, prod)
##   pdf
## }

dfrankCopula.pdf <- function(u, copula, log=FALSE) {
  dim <- copula@dimension
  if (dim > 6) stop("Frank copula PDF not implemented for dimension > 6.")
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  for (i in 1:dim) assign(paste0("u", i), u[,i])
  alpha <- copula@parameters[1]
  ## FIXME: improve log-case
  if(log)
    log(c(eval(frankCopula.pdf.algr[dim])))
  else  c(eval(frankCopula.pdf.algr[dim]))
}


tauFrankCopula <- function(copula) .tauFrankCopula(copula@parameters)
.tauFrankCopula <- function(a) { # 'a', also called 'alpha' or 'theta'
  ## For small a (not just a = 0), need Taylor approx:
  ##  1 - 4 / a * (1 - debye1(a)) = a/9 *(1 - (a/10)^2  + O(a^4))  <<- MM, 16.May 2014
  if(abs(a) < 10*sqrt(.Machine$double.eps))
    a/9
  else
    1 - 4 / a * (1 - debye1(a))
}

rhoFrankCopula <- function(copula) .rhoFrankCopula(copula@parameters)
.rhoFrankCopula <- function(a) { # 'a', also called 'alpha' or 'theta'
  ## For small alpha (not just alpha = 0), need Taylor approx:
  ## Genest (1987), "Frank's family of bivariate distributions" (Biometrika, 7, 549--555)
  ##  1 - 12/a * (debye1(a) - debye2(a)) =
  ##  1 + 12/a * (debye2(a) - debye1(a)) = a/6 * (1 - a^2 / 75 + O(a^4)) <<- MM, 16.May 2014
  if(abs(a) < sqrt(75*.Machine$double.eps))
    a/6
  else
    1 + 12/a * (debye2(a) - debye1(a))
}

dTauFrankCopula <- function(copula) .dTauFrankCopula(copula@parameters)
.dTauFrankCopula <- function(a) {
  ## FIXME: use Taylor for small |a|
  (2/a)^2 * (a/expm1(a) + 1 - 1/a * debye1(a))
}

dRhoFrankCopula <- function(copula) .dRhoFrankCopula(copula@parameters)
.dRhoFrankCopula <- function(a) {
  ## FIXME: use Taylor for small |a|
  12 / a^2 * (a/expm1(a) - 3 * debye2(a) + 2 * debye1(a))
}


pMatFrank <- function (u, copula, ...) {
    ## was  pfrankCopula
    stopifnot(!is.null(d <- ncol(u)), d == copula@dimension)
    th <- copula@parameters
    if(d == 2 && !copFrank@paraConstr(th)) # for now, .. to support negative tau
        pfrankCopula(copula, u=u)
    else
        pacopula(u, copFrank, theta=copula@parameters, ...)
}

dMatFrank <- function (u, copula, log = FALSE, checkPar=TRUE, ...) {
    ## was  dfrankCopula.pdf
    stopifnot(!is.null(d <- ncol(u)), d == copula@dimension)
    th <- copula@parameters
    if(d == 2 && th < 0) # for now, copFrank does not yet support negative tau
        dfrankCopula.pdf(u, copula, log=log)
    else
        copFrank@dacopula(u, theta=th, log=log, checkPar=checkPar, ...)
}
setMethod("rCopula", signature("numeric", "frankCopula"), rfrankCopula)

setMethod("pCopula", signature("numeric", "frankCopula"),
	  function (u, copula, ...)
	  pMatFrank(matrix(u, ncol = dim(copula)), copula, ...))
setMethod("pCopula", signature("matrix", "frankCopula"), pMatFrank)

setMethod("dCopula", signature("numeric", "frankCopula"),
	  function (u, copula, log=FALSE, ...)
	  dMatFrank(matrix(u, ncol = dim(copula)), copula, log=log, ...))
setMethod("dCopula", signature("matrix", "frankCopula"), dMatFrank)

setMethod("iPsi", signature("frankCopula"), iPsiFrank)
setMethod("psi",  signature("frankCopula"),  psiFrank)

## setMethod("psiDer", signature("frankCopula"), psiDerFrank)
setMethod("diPsi", signature("frankCopula"),
	  function(copula, u, degree=1, log=FALSE, ...)
      {
	  s <- if(log || degree %% 2 == 0) 1. else -1.
	  s* copFrank@absdiPsi(u, theta=copula@parameters, degree=degree, log=log, ...)
      })



setMethod("tau", signature("frankCopula"), tauFrankCopula)
setMethod("rho", signature("frankCopula"), rhoFrankCopula)
setMethod("tailIndex", signature("frankCopula"), function(copula) c(lower=0, upper=0))

setMethod("iTau", signature("frankCopula"),
	  function(copula, tau, tol = 1e-7) copFrank@iTau(tau, tol=tol))
setMethod("iRho", signature("frankCopula"), iRhoCopula)

setMethod("dRho", signature("frankCopula"), dRhoFrankCopula)
setMethod("dTau", signature("frankCopula"), dTauFrankCopula)

# Copied from flexmirt.R
Tnom.trend <- function(nc) {
  T <- matrix(0,nc,nc-1)
  for (i in 1:nc) {
    T[i,1] <- i-1
  }
  for (k in 2:(nc-1)) {
    for (i in 2:(nc-1)) {
      T[k,i] <- sin(pi*(i-1)*(k-1)/(nc-1))
    }
  }
  return(T[-1,])
}

# Copied from flexmirt.R
Tnom.ida <- function(nc) {
  T <- matrix(0,nc,nc-1)
  T[nc,1] <- nc-1
  T[2:(nc-1),2:(nc-1)] <- diag(nc-2)
  return(T[-1,])
}

# Copied from flexmirt.R
Tnom.idc <- function(nc) {
  T <- matrix(0,nc,nc-1)
  T[2:nc,] <- diag(nc-1)
  return(T[-1,])
}

build.T <- function(outcomes, got, type) {
	if (is.null(got)) stop("rpf.nrm: T matrix specified as NULL")
  if (!is.matrix(got)) {
    if (got == "id") {
	    if (type == 'a') {
		    got <- Tnom.ida(outcomes)
	    } else {
		    got <- Tnom.idc(outcomes)
	    }
    } else if (got == "trend") {
      got <- Tnom.trend(outcomes)
    } else if (got == "random") {
      while (1) {
        side <- outcomes-1
        got <- matrix(rnorm(side*side), side, side)
        invertible <- try(solve(got), silent=TRUE)
        if (!inherits(invertible, "try-error")) break
      }
    } else {
      stop(paste("T matrix", deparse(got), "not recognized"))
    }
  }
  if (all(dim(got) == c(outcomes,outcomes-1))) {
    if (any(got[1,] != 0)) warning("Non-zero T[1,] will be ignored")
    got <- got[-1,]
  }
  if (all(dim(got) != rep(outcomes-1, 2))) {
    stop(paste("T matrix must be of dimension",
               paste(rep(outcomes-1, 2), collapse="x"),
               "not", paste(dim(got), collapse="x")))
  }
  got
}

##' Create a nominal response model
##'
##' This function instantiates a nominal response model.
##'
##' The transformation matrices T.a and T.c are chosen by the analyst
##' and not estimated.  The T matrices must be invertible square
##' matrices of size outcomes-1. As a shortcut, either T matrix
##' can be specified as "trend" for a Fourier basis or as "id" for an
##' identity basis. The response probability function is
##'
##' \deqn{a = T_a \alpha}
##' \deqn{c = T_c \gamma}
##' \deqn{\mathrm P(\mathrm{pick}=k|s,a_k,c_k,\theta) = C\ \frac{1}{1+\exp(-(s \theta a_k + c_k))}}
##'
##' where \eqn{a_k} and \eqn{c_k} are the result of multiplying two vectors
##' of free parameters \eqn{\alpha} and \eqn{\gamma} by fixed matrices \eqn{T_a} and \eqn{T_c}, respectively;
##' \eqn{a_0} and \eqn{c_0} are fixed to 0 for identification;
##' and \eqn{C} is a normalizing factor to ensure that \eqn{\sum_k \mathrm P(\mathrm{pick}=k) = 1}.
##' 
##' @param outcomes The number of choices available
##' @param factors the number of factors
##' @param T.a the T matrix for slope parameters
##' @param T.c the T matrix for intercept parameters
##' @return an item model
##' @references Thissen, D., Cai, L., & Bock, R. D. (2010). The
##' Nominal Categories Item Response Model. In M. L. Nering &
##' R. Ostini (Eds.), \emph{Handbook of Polytomous Item Response
##' Theory Models} (pp. 43--75). Routledge.
##' @examples
##' spec <- rpf.nrm()
##' rpf.prob(spec, rpf.rparam(spec), 0)
##' # typical parameterization for the Generalized Partial Credit Model
##' gpcm <- function(outcomes) rpf.nrm(outcomes, T.c=lower.tri(diag(outcomes-1),TRUE) * -1)
##' spec <- gpcm(4)
##' rpf.prob(spec, rpf.rparam(spec), 0)
rpf.nrm <- function(outcomes=3, factors=1, T.a="trend", T.c="trend") {
  if (outcomes < 3) stop("Minimum number of outcomes is 3")
  T.a <- build.T(outcomes, T.a, 'a')
  T.c <- build.T(outcomes, T.c, 'c')
  id <- rpf.id_of("nominal")
  m <- new("rpf.mdim.nrm",
           outcomes=outcomes,
           factors=factors)
  m@spec <- c(id, outcomes, factors, T.a, T.c, solve(T.a), solve(T.c))
  m
}

getT <- function(m, tx) {
  Tsize <- (m@outcomes-1L)^2
  offset <- 4 + tx * Tsize
  matrix(m@spec[offset:(offset+Tsize-1)], m@outcomes-1L, m@outcomes-1L)
}

setMethod("rpf.rparam", signature(m="rpf.mdim.nrm"),
          function(m, version) {
		  if (m@factors == 0) {
			  return(c(gam=ck <- sort(rnorm(m@outcomes-1))))
		  }
            a <- rlnorm(m@factors, sdlog=.5)
            ak <- abs(rnorm(m@outcomes-1, mean=1, sd=.25))
            ck <- sort(rnorm(m@outcomes-1))
            c(a=a,
              alf=getT(m,2) %*% ak,
              gam=getT(m,3) %*% ck)
          })

setMethod("rpf.modify", signature(m="rpf.mdim.nrm", factors="numeric"),
          function(m, factors) {
              rpf.nrm(m@outcomes, factors, T.a=getT(m,0), T.c=getT(m,1))
          })

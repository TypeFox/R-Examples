
#' Derive variable ranges from linear restrictions
#'
#' Gaussian and/or Fourier-Motzkin elimination is used 
#' to derive upper and lower limits implied by a system of (in)equations.
#' 
#'   
#' @param A \code{[numeric]} Matrix 
#' @param b \code{[numeric]} vector
#' @param neq [\code{numeric}] The first \code{neq} rows in \code{A} and
#'   \code{b} are treated as linear equalities. 
#' @param nleq [\code{numeric}] The \code{nleq} rows after \code{neq} are treated as
#' inequations of the form \code{a.x<=b}. All remaining rows are treated as strict inequations
#' of the form \code{a.x<b}.
#' @param eps \code{[numeric]} Coefficients with absolute value  \code{<= eps} are treated as zero.
#' using Fourier-Motzkin elimination.
#'
#' @export
ranges <- function(A, b, neq=nrow(A), nleq=0, eps=1e-8){
  nvar <- ncol(A)
  R <- array(NA, dim=c(nvar,2),dimnames=list(variable=1:nvar, range=c("lower","upper")))
  I <- 1:nvar
  for ( i in I ){
    L <- eliminate_variables(A=A, b=b, variables=setdiff(I,i), neq=neq, nleq=nleq, eps=eps)
    L <- compact(L$A,L$b,neq=L$neq, nleq=L$nleq, eps=eps)
    if ( nrow(L$A) == 0 ){ # no restrictions
      lim <- c(-Inf,Inf) 
    } else if (L$neq > 0){ # at least one equation
      lim <- c(L$b[1]/L$A[1,1], L$b[1]/L$A[1,1])
    } else { # all inequations
      ilw <- L$A[,1] < 0
      lwr <- infmax(L$b[ilw]/L$A[ilw,1])
      iup <- L$A[,1] > 0
      upr <- infmin(L$b[iup]/L$A[iup,1])
      lim <- c(lwr,upr)
    }
    R[i,] <- lim
  }
  R
}

# if numeric(0) is passed, min and max have desired behaviour: returning
# Inf, -Inf, but with a warning we can avoid.
infmax <- function(x) suppressWarnings(max(x))
infmin <- function(x) suppressWarnings(min(x))

eliminate_variables <- function(A, b, variables, neq=nrow(A),nleq=0,eps=1e-8){
  
  L <- list(A=A, b=b, neq=neq, nleq=nleq, H=NULL, h=0)
  for ( v in variables ){
    L <- eliminate(A=L$A, b=L$b, neq=L$neq, nleq=L$nleq, variable=v, H=L$H, h=L$h)
  }
  L
}



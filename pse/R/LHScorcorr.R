#' Correlation Matrix Correction
#' 
#' Corrects the correlation matrix of a given Latin Hypercube Sample.
#' 
#' This function changes the order in which data is organized in order
#' to force the correlation matrix to a prescribed value. This implementation
#' uses the Hungtington-Lyrintzis algorithm.
#'
#' This is mainly intended for use inside of the \code{\link{LHS}} function.
#'
#' If you intend to use non-zero correlation terms, read Chalom & Prado (2012) for some
#' important theoretical restrictions.
#' 
#' The correlation matrix may be specified by a Pearson or Spearman method. In order to
#' generate the Spearman correlation, the function "rank transforms" the data using the
#' \code{\link[base]{order}} function, and thus works only if there are no ties in the data.
#'
#' @param vars The data.frame or matrix containing the parameters from the "raw" Latin Hypercube Sample.
#'	  Each column corresponds to one variable, and each line to one observation.
#'
#' @param COR The desired correlation matrix. The default is to have 0 correlation.
#'	  You can supply a numeric square matrix with M rows, where M is the number of
#'	  input factors. The *lower* triangular part of the matrix will be used as the 
#'	  desired correlation matrix.
#' @param method A character string, which may be "Spearman" or "Pearson", indicating the 
#' correlation method to be used.
#' @param eps The tolerance for the deviation between the prescribed correlation matrix and the result.
#' @param echo Set to true to display information messages.
#' @param maxIt Maximum number of iterations before giving up. 
#' Set to 0 to use a heuristic based on the size of the hypercube.
#' Set to a negative number to never give up. *CAUTION*, this might result in an infinite loop.
#' @return A data.frame containing the same variables, but with the correlation matrix corrected.
#' @references 
#'  Huntington, D.E. and Lyrintzis, C.S. 1998 Improvements to and limitations of 
#'  Latin hypercube sampling. \emph{Prob. Engng. Mech.} 13(4): 245-253.
#'
#'  Chalom, A. and Prado, P.I.K.L. 2012. Parameter space exploration of ecological models
#'  \emph{arXiv}:1210.6278 [q-bio.QM]
#' @export
#' @import stats
#' @useDynLib pse, corcorr
LHScorcorr <-
	function (vars, COR = 0, method=c("Pearson", "Spearman"), eps = 0.005, echo=FALSE, maxIt = 0) {
    method <- match.arg(method)
		if (! is.matrix(COR)) COR = matrix(0,dim(vars)[2],dim(vars)[2])
		if (maxIt==0) maxIt = 2*sqrt(dim(vars)[1])
    if (method=="Pearson") { # Simply proceeds with the calculation
  		return(internal.LHScorcorr(vars, COR, 2, eps, 1, echo, maxIt))
    } else { # "Rank" transforms the data, calculates the new arrangement, then transforms back
      Rvars <- apply(vars, 2, order)
      if (! all.equal(Rvars, apply(vars, 2, rank)))
        warning("Ties detected on LHScorcorr, result may not be exact")
  		Rvars <- internal.LHScorcorr(Rvars, COR, 2, eps, 1, echo, maxIt)
      for (i in 1:(dim(vars)[2]))
        vars[, i] <- vars[Rvars[,i] , i]
      return(vars)
	}
}

internal.LHScorcorr <- function (vars, COR, l, eps, it, echo, maxIt) {
  N <- dim(vars)[1]; M <- dim(vars)[2]
  my.names <- colnames(vars)
  # Stop condition: 
  if (l == M + 1) {
    return (vars);
  }
  # Skipping this column: correlations are close enough to the prescribed
  if (max(abs(cor(vars)[l,1:(l-1)] - COR[l,1:(l-1)])) < eps) {
    return (internal.LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1, echo=echo, maxIt=maxIt));
  }
  # Skipping this column: maximum iterations for the same variable
  # If maxIt is set to -1, NEVER gives up
  if (it > maxIt | maxIt < 0) { 
    warning("LHScorcorr: correlation does not converge after maximum iterations");
    return (internal.LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1, echo=echo, maxIt=maxIt));
  }
  if (echo==T) cat(paste("Info: Correlation correction being made for l =",l,"/",M, "\n"))
  # Here we start correcting the correlation for var[,l]
  V <- .C(corcorr, vars=as.double(as.matrix(vars)),cor=as.double(COR), N=as.integer(N), M=as.integer(M), l=as.integer(l), FLAGSTOP=as.integer(0))
  vars <- as.data.frame(matrix(V$vars, nrow=N, ncol=M))
  names(vars) <- my.names
  if (V$FLAGSTOP == 1) { # Convergence, going for next
    return(internal.LHScorcorr (vars, COR, l = l + 1, eps = eps, it = 1, echo=echo, maxIt=maxIt));
  } else {
    # Repeat the proccess with the same variable
    internal.LHScorcorr(vars,COR,l=l,eps=eps,it=it+1, echo=echo, maxIt=maxIt)
  }
}

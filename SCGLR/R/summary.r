#' @export
#' @importFrom stats cor
#' @importFrom utils combn
#' @title Summarizing SCGLR fits
#' @description Summary method for class "SCGLR".
#' @method summary SCGLR
#' @param object an object of class "SCGLR", usually a result of a call to \code{\link{scglr}}.
#' @param \dots Not used.
#' @return
#' an object of class "summary.SCGLR".
#'    \item{inertia}{inertia per component.}
#'    \item{deviance}{deviance for each \eqn{Y_k}.}
#'    \item{rho}{squared correlations with numerical covariates.}
#'    \item{rho.pred}{squared correlations with linear predictors.}
#'    \item{coefficients}{contains the coefficients of the regression on the components.}
#'    \item{pvalue}{contains the pvalues of the coefficients of the regression on the components.}
summary.SCGLR <- function(object, ...) {
  rho <- as.data.frame(cor(object$xNumeric,object$compr))
  rho.pred <- as.data.frame(cor(object$lin.pred,object$compr))
  ncomp <- ncol(object$compr)
  if(ncomp>1) {
    cmp_pairs <- combn(ncomp, 2, simplify=FALSE)
  
    best_plane <- function(var) {
      magni <- lapply(cmp_pairs, function(pair) sum(var[pair]^2))
      ind <- which.max(magni)
      list(bp=paste(as.character(cmp_pairs[[ind]]),collapse="/"),val=magni[ind])
    }
    tmp <- unlist(apply(rho,1,best_plane))
    tmp.pred <- unlist(apply(rho.pred,1,best_plane))
    rho <- data.frame(rho^2,best_plane=tmp[seq(1,length(tmp),2)],best_val=as.numeric(tmp[seq(2,length(tmp),2)])) 
    rho <- rho[order(rho$best_val,decreasing=TRUE),]
    rho.pred <- data.frame(rho.pred^2,best_plane=tmp.pred[seq(1,length(tmp.pred),2)],best_val=as.numeric(tmp.pred[seq(2,length(tmp.pred),2)]))
    rho.pred <- rho.pred[order(rho.pred$best_val,decreasing=TRUE),]
  } else {
    rho <- rho^2
    rho <- rho[order(rho[,1],decreasing=TRUE),,drop=FALSE]
    rho.pred <- NA
  }
  coef <- sapply(object$gamma, function(x) x[,1])
  pvalue <- sapply(object$gamma, function(x) x[,4])
  colnames(coef) <- colnames(object$lin.pred)
  
  structure(list(
    call=object$call,
    inertia=object$inertia,
    deviance=object$deviance,
    rho=rho,
    rho.pred=rho.pred,
    coefficients=coef,
    pvalue=pvalue),
  class="summary.SCGLR")
}

#' @export
#' @rdname summary.SCGLR
#' @method print summary.SCGLR
#' @param x an object of class "summary.SCGLR", usually a result of a call to summary.SCGLR.
#' @param digits the number of significant digits to use when printing.
#' @param cutoff print coefficients with pvalue lower than or equal to cutoff (default to 1).
print.summary.SCGLR <- function(x, digits=3, cutoff=1, ...) {

  cat("Squared correlations with numerical covariates (in decreasing order):\n")
  print(x$rho, print.gap=2, digits=digits)
  
  if( is.data.frame(x$rho.pred) ) {
    cat("\nSquared correlations with linear predictors (in decreasing order):\n")  
    print(x$rho.pred ,print.gap=2, digits=digits)
  }
  
  coef <- x$coefficients
  coef[x$pvalue>=cutoff] <- NA
  
  cat("\nCoefficients for dependant variables:\n")
  print(coef, na.print="", digits=digits)
  cat("\n")
  invisible(x)
}

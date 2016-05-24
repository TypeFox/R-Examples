summary.NoncompMAR <- function(object, CI=c(2.5, 97.5),...){

  qoi <- cbind(object$itt, object$cace, object$pc, object$base)
  qoi <- cbind(apply(qoi, 2, mean), apply(qoi, 2, sd),
               apply(qoi, 2, quantile, min(CI)/100),
               apply(qoi, 2, quantile, max(CI)/100))
  colnames(qoi) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                     paste(max(CI), "%", sep=""))
  
  if (!is.null(object$coefficientsC)) {
    coefC <- object$coefficientsC
    coefC <- cbind(apply(coefC, 2, mean), apply(coefC, 2, sd),
                   apply(coefC, 2, quantile, min(CI)/100),
                   apply(coefC, 2, quantile, max(CI)/100))
    colnames(coefC) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                         paste(max(CI), "%", sep=""))
  }
  else
    coefC <- NULL

  if (!is.null(object$coefficientsO)) {
    coefO <- object$coefficientsO
    coefO <- cbind(apply(coefO, 2, mean), apply(coefO, 2, sd),
                   apply(coefO, 2, quantile, min(CI)/100),
                   apply(coefO, 2, quantile, max(CI)/100)) 
    colnames(coefO) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                         paste(max(CI), "%", sep=""))
  }
  else
    coefO <- NULL

  if (!is.null(object$coefficientsS)) {
    coefS <- object$coefficientsS
    coefS <- cbind(apply(coefS, 2, mean), apply(coefS, 2, sd),
                   apply(coefS, 2, quantile, min(CI)/100),
                   apply(coefS, 2, quantile, max(CI)/100)) 
    colnames(coefS) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                         paste(max(CI), "%", sep=""))
  }
  else
    coefS <- NULL

  
  ans <- list(call = object$call, n.obs = length(object$Y), n.draws =
              object$n.draws, qoi.table = qoi, coefC.table = coefC,
              coefO.table = coefO, coefS.table = coefS)  

  if (!is.null(object$thresholds)) {
    tauO <- object$thresholds
    tauO <- cbind(apply(tauO, 2, mean), apply(tauO, 2, sd),
                  apply(tauO, 2, quantile, min(CI)/100),
                  apply(tauO, 2, quantile, max(CI)/100)) 
    colnames(tauO) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                        paste(max(CI), "%", sep=""))
    ans$tauO.table <- tauO
  }

  class(ans) <- "summary.NoncompMAR"
  return(ans)
}

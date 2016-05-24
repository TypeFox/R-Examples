##' @export
summary.pec <-  function(object,
                           times,
                           what,
                           models,
                           digits=3,
                           print=TRUE,
                           ...){
  
  if (missing(models)) models <- names(object$models)
  if (missing(what) || is.null(what)){
    what <- grep(c("Err$"),names(object),value=TRUE)
  }
  if (print==TRUE) cat("\nPrediction error curves\n\n")
  if (print==TRUE) print(object$splitMethod)
  otime <- object$time
  if (missing(times) && (length(times <- otime) > 20)){
    warning("Missing times argument: prediction error curves evaluated at the quantiles of fitted times\n")
    times <- quantile(otime)
  }
  tindex <- prodlim::sindex(jump.times=object$time,eval.times=times)
  out <- lapply(what,function(w){
    if (print==TRUE) cat("\n",w,"\n")
    tmp <- rbind(0, do.call("cbind",object[[w]][models]))[tindex+1,,drop=FALSE]
    tmp <- cbind(time=times,n.risk=c(object$n.risk[1],object$n.risk)[tindex+1],tmp)
    rownames(tmp) <- 1:NROW(tmp)
    if (print==TRUE) prmatrix(round(tmp,digits=digits),...)
    tmp
  })
  names(out) <- what
  if (!is.null(object$multiSplitTest))
    if (print==TRUE) print(object$multiSplitTest)
  if (print==TRUE) cat("\n")
  invisible(out)
}

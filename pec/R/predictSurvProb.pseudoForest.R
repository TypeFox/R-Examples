##' @export
predictSurvProb.pseudoForest <- function(object,
                                         newdata,
                                         times,
                                         digits=8,
                                         ...){
  stopifnot(object$model.type=="survival")
  # {{{ get forests

  # Extract forests
  #  - NOTE: if more than one the forests are independent and only
  #          made for easy extraction to different time points.
  #
  forestList <- object$forest
  if (class(forestList[[1]])!="randomForest")
    stop("Only works for 'randomForest'")

  L <- length(forestList)

  # }}}
  # {{{ predict to given time points
  # find the forest
  ## pos <- prodlim::sindex(jump.times=object$times,eval.times=times)
  pos <- match(times,object$times,nomatch=FALSE)
    if (any(pos==FALSE))
    stop("Requested forests at times ",paste(times[!pos],collapse=", "),"not available. Available are forests at times:",paste(object$times,collapse=", "))
  p <- do.call("cbind",lapply(pos,function(t){
    getForest <- forestList[[t]]
    ## print(names(newdata))
    ## print(str(getForest))
    p.t <- stats::predict(getForest,newdata=newdata)
    p.t <- round(p.t,digits=digits)
  }))
  # }}}
  # {{{ return
  # check dim.
  if (is.null(dim(p))) {
    if (length(p)!=length(times))
       stop("Prediction failed")
  }
  else{
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  }
  p
  # }}}
}

##' @export
predictEventProb.pseudoForest <- function(object,
                                          newdata,
                                          times,
                                          cause,
                                          ...){
    stopifnot(object$model.type=="competing.risks")
    # {{{ cause
    if (missing(cause)) cause <- 1
    # }}}
    # {{{ get forests
    forestList <- object$forest
    ##   if (class(forestList[[1]])!="randomForest")
    ##     stop("Only works for 'randomForest'")
    L <- length(forestList)
    # }}}
    # {{{ predict to given time points
    # find the forest
    pos <- match(times,object$times,nomatch=FALSE)
    if (any(pos==FALSE))
        stop("Requested forests at times ",paste(times[!pos],collapse=", "),"not available. Available are forests at times:",paste(object$times,collapse=", "))
    ##   pos <- prodlim::sindex(jump.times=object$times,eval.times=times)
    p <- do.call("cbind",lapply(pos,function(t){
        getForest <- forestList[[t]]
        if (class(getForest)!="randomForest")
            pseudo.t <- rep(getForest, NROW(newdata))
        else
            pseudo.t <- stats::predict(getForest,newdata=newdata)
        ## pseduo.t <- round(pseudo.t,digits=digits)
        # Pseudo-value correction: return only [0;1]
        no.negative <- pmax(pseudo.t,0)
        under.one <- pmin(no.negative,1)
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

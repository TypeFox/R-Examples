##' @export
summary.Cindex <- function(object,what=NULL,digits=3,print=TRUE,...){
    ccr <- attr(object$response,"model")=="competing.risks"
    if (print==TRUE) cat("\nThe c-index for right censored event times\n\n")
    # {{{ echo models
    if (print==TRUE) cat("Prediction models:\n\n")
    printModels <- sapply(object$models,function(m){
        if (class(m) %in% c("character","call"))
            m
        else
            if (class(try(m$call,silent=TRUE))=="try-error")
                "unknown formula"
            else
                m$call
    })
    if (print==TRUE) print(printModels,quote=FALSE)
    # }}}
    # {{{ echo response
    if (print==TRUE) print(object$response)
    # }}}
    # {{{ echo cens model and splitmethod
    if (print==TRUE){
        if (!is.null(object$cens.model)){
            if (object$cens.model!="none")
                cat("\nCensoring model for IPCW:",object$cens.model,"model",ifelse(object$cens.model=="marginal","(Kaplan-Meier for censoring distribution)",""),"\n")
            else cat("\nno censoring")}
        if (!is.null(object$splitMethod)) print(object$splitMethod)
    }
    # }}}
    # {{{ discover what to print
    if (missing(what) || is.null(what)){
        what <- grep(c("Cindex$"),names(object),value=TRUE)
    }
    # }}}
    # {{{ result table
    out <- lapply(what,function(r){
        out <- do.call("rbind",lapply(1:length(object$models),function(m){
            object[[r]][[m]]
        }))
        if (is.matrix(out)){
            rownames(out) <- names(object$models)
            coln <- paste("time=",round(object$time,1),sep="")
            coln[object$time<1] <- paste("time=",round(object$time[object$time<1],4),sep="")
            colnames(out) <- coln
        }
        out
    })
    names(out) <- what
    # {{{ if only one time point
    if (NCOL(out[[1]])==1){
        if (print==TRUE)    cat("\nEstimated C-index in % at",colnames(out[[1]]),"\n\n")
        outMat <- 100*do.call("cbind",out)
        colnames(outMat) <- what
        if (ccr){
            if (!is.null(object$Pairs))
                outMat <- cbind(outMat,
                                "Pairs (Di=1,Ti<Tj)"=round(unlist(sapply(object$Pairs,function(x)x$A),1)),
                                "Concordant"=round(unlist(sapply(object$Concordant,function(x)x$A),1)),
                                "Pairs (Di=1,Dj=2)"=round(unlist(sapply(object$Pairs,function(x)x$B),1)),
                                "Concordant"=round(unlist(sapply(object$Concordant,function(x)x$B),1)))
        }
        else{
            if (!is.null(object$Pairs))
                outMat <- cbind(outMat,
                                Pairs=round(object$Pairs[[1]],1),
                                Concordant=round(unlist(sapply(object$Concordant,function(x)x),1)))
        }
        if (print==TRUE)    print(outMat,digits)
    }
    # }}}
    # {{{ multiple time points
    else{
        if (print==TRUE)     cat("\nEstimated C-index in %\n\n")
        if (print==TRUE) print(lapply(out,function(x)x*100),digits)
    }
    if(object$splitMethod$name=="BootCv")
        if (print==TRUE) cat("\nAppCindex    : Apparent (training data) performance\nBootCvCindex : Bootstrap crossvalidated performance\n\n")
    # }}}
    invisible(out)
}


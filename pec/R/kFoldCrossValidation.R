kFoldCrossValidation <- function(object,
                                 data,
                                 Y,
                                 status,
                                 event,
                                 times,
                                 cause,
                                 ipcw,
                                 splitMethod,
                                 giveToModel,
                                 predictHandlerFun,
                                 keep,
                                 verbose){
    # {{{ initializing
    B <- splitMethod$B
    N <- splitMethod$N
    k <- splitMethod$k
    NT <- length(times)
    NF <- length(object) 
    ResampleIndex <- splitMethod$index
    if (missing(giveToModel)) extraArgs=NULL
    # }}}
    CrossValErrMat <- lapply(1:B,function(b,extraArgs=extraArgs){
        if (verbose==TRUE) internalTalk(b,B)
        groups <- ResampleIndex[,b,drop=TRUE]
        ## each subject belongs to exactly one group
        ## the prediction `p[i]' is obtained with the reduced data
        if (k==N-1) k <- N
        subjectPred <- lapply(1:k,function(g){
            if (verbose==TRUE) internalTalk(g,k)
            # {{{ training and validation data
            id <- groups==g
            train.k <- data[!id,,drop=FALSE]
            val.k <- data[id,,drop=FALSE]
            # }}}
            # {{{ Building the models in training data
            trainModels <- lapply(1:NF,function(f){
                fit.k <- internalReevalFit(object=object[[f]],data=train.k,step=paste("CV group=",k),silent=FALSE,verbose=verbose)
                ## this was a good idea to reduce the memory usage:
                ## fit.k$call <- object[[f]]$call
                ## fit.k$call <- NULL
                ## however, it does not work with the new version of the survival package
                ## in which the survfit.coxph function checks the response 'y'
                ## fit.k$call$data <- substitute(train.k)
                fit.k
            })
            # }}}
            # {{{ Predicting the validation data
            modelPred <- lapply(1:NF,function(f){
                fit.k <- trainModels[[f]]
                extraArgs <- giveToModel[[f]]
                if (predictHandlerFun == "predictEventProb"){      
                    p.group <- do.call(predictHandlerFun,c(list(object=fit.k,newdata=val.k,times=times,cause=cause),extraArgs))
                }
                else{
                    p.group <- do.call(predictHandlerFun,c(list(object=fit.k,newdata=val.k,times=times),extraArgs))
                }
                if(is.null(dim(p.group))) {
                    p.group <- do.call("rbind",lapply(1:NROW(val.k),function(x){p.group}))
                }
                p.group
            })
            # }}}
            modelPred
        })
        # {{{ Compute prediction error curves for step b
        pec <- lapply(1:NF,function(f){
            pred.b <- do.call("rbind",lapply(subjectPred,function(x)x[[f]]))
            if (splitMethod$internal.name!="loocv"){
                pred.b <- pred.b[order(order(groups)),]
            }
            if (predictHandlerFun=="predictEventProb")
                innerCrossValErr <- .C("pecCR",
                                       pec=double(NT),
                                       as.double(Y),
                                       as.double(status),
                                       as.double(event),
                                       as.double(times),
                                       as.double(pred.b),
                                       as.double(ipcw$IPCW.times),
                                       as.double(ipcw$IPCW.subjectTimes),
                                       as.integer(N),
                                       as.integer(NT),
                                       as.integer(ipcw$dim),
                                       as.integer(is.null(dim(pred.b))),
                                       NAOK=TRUE,
                                       PACKAGE="pec")$pec
            else
                innerCrossValErr <- .C("pec",pec=double(NT),as.double(Y),as.double(status),as.double(times),as.double(pred.b),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$pec
            innerCrossValErr
        })
        names(pec) <- names(object)
        pec
    })
    # }}}
    if (B>1){
        CrossValErr <- lapply(1:NF,function(f){
            rowMeans(do.call("cbind",lapply(CrossValErrMat,function(b)b[[f]])))
        })
    }
    else
        CrossValErr <- CrossValErrMat[[1]]
    out <- list(CrossValErr=CrossValErr)
    if (keep==TRUE && B>1) out$CrossValErrMat <- CrossValErrMat
    out
}

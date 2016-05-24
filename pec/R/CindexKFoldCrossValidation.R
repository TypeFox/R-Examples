CindexKFoldCrossValidation <- function(object,
                                       data,
                                       Y,
                                       status,
                                       event,
                                       tindex,
                                       eval.times,
                                       pred.times,
                                       cause,
                                       weights,
                                       ipcw.refit=FALSE,
                                       ipcw.call,
                                       tiedPredictionsIn,
                                       tiedOutcomeIn,
                                       tiedMatchIn,
                                       splitMethod,
                                       multiSplitTest,
                                       keepResiduals,
                                       testTimes,
                                       confInt,
                                       confLevel,
                                       getFromModel,
                                       giveToModel,
                                       predictHandlerFun,
                                       keepMatrix,
                                       verbose,
                                       savePath,slaveseed){
    # {{{ initializing
  B <- splitMethod$B
  N <- splitMethod$N
  k <- splitMethod$k
  NT <- length(eval.times)
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
            internalTalk(g,k)
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
          p.group <- do.call(predictHandlerFun,c(list(object=fit.k,newdata=val.k,times=eval.times,cause=cause),extraArgs))
        }
        else{
          p.group <- do.call(predictHandlerFun,c(list(object=fit.k,newdata=val.k,times=eval.times),extraArgs))
        }
        if(is.null(dim(p.group))) {
          p.group <- do.call("rbind",lapply(1:NROW(val.k),function(x){p.group}))
        }
        p.group
      })

            # }}}
            modelPred
        })
        ipcw.i <- weights$weight.i
        if (is.null(dim(weights$weight.j))){
            ipcw.j <- weights$weight.j
        }
        else{
            ipcw.j <- weights$weight.j
        }
        # {{{ Compute cindex for step b
        PredCindexStepB  <- lapply(1:NF,function(f){
            pred.b <- do.call("rbind",lapply(subjectPred,function(x)x[[f]]))
            if (splitMethod$internal.name!="loocv"){
                pred.b <- pred.b[order(order(groups)),]
            }
            if (predictHandlerFun=="predictEventProb"){
                Step.b.CindexResult <- .C("ccr",
                                          cindex=double(NT),
                                          concA=double(NT),
                                          pairsA=double(NT),
                                          concB=double(NT),
                                          pairsB=double(NT),
                                          as.integer(tindex),
                                          as.double(Y),
                                          as.integer(status),
                                          as.integer(event),
                                          as.double(eval.times),
                                          as.double(ipcw.i),
                                          as.double(ipcw.j),
                                          as.double(pred.b),
                                          as.integer(N),
                                          as.integer(NT),
                                          as.integer(tiedPredictionsIn),
                                          as.integer(tiedOutcomeIn),
                                          as.integer(tiedMatchIn),
                                          as.integer(!is.null(dim(ipcw.j))),
                                          NAOK=TRUE,
                                          package="pec")
                Step.b.Cindex <- Step.b.CindexResult$cindex
                Step.b.PairsA <- Step.b.CindexResult$pairsA
                Step.b.ConcordantA <- Step.b.CindexResult$concA
                Step.b.PairsB <- Step.b.CindexResult$pairsB
                Step.b.ConcordantB <- Step.b.CindexResult$concB
                list(Cindex.b=Step.b.Cindex,Pairs.b=list(A=Step.b.PairsA,B=Step.b.PairsB),Concordant.b=list(A=Step.b.ConcordantA,B=Step.b.ConcordantB))
            }
            else{
                cindexOut <- .C("cindex",
                                cindex=double(NT),
                                conc=double(NT),
                                pairs=double(NT),
                                as.integer(tindex),
                                as.double(Y),
                                as.integer(status),
                                as.double(eval.times),
                                as.double(ipcw.i),
                                as.double(ipcw.j),
                                as.double(pred.b),
                                as.integer(N),
                                as.integer(NT),
                                as.integer(tiedPredictionsIn),
                                as.integer(tiedOutcomeIn),
                                as.integer(tiedMatchIn),
                                as.integer(!is.null(dim(ipcw.j))),
                                NAOK=TRUE,
                                package="pec")            
                Cindex.b <- cindexOut$cindex
                Pairs.b <- cindexOut$pairs 
                Concordant.b <- cindexOut$conc
                list(Cindex.b=Cindex.b,Pairs.b=Pairs.b,Concordant.b=Concordant.b)
            }
        })
        names(PredCindexStepB) <- names(object)
        PredCindexStepB
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
  if (keepMatrix==TRUE && B>1) out$CrossValErrMat <- CrossValErrMat
  out
}

CindexBootstrapCrossValidation <- function(object,
                                           data,
                                           Y,
                                           status,
                                           event,
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
  M <- splitMethod$M
  NT <- length(eval.times)
  NF <- length(object) 
  ResampleIndex <- splitMethod$index
  # }}}
  step <- function(b,seed){
    if (verbose==TRUE) internalTalk(b,B)
    # {{{ training and validation data
    vindex.b <- match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0
    Y.b <- Y[vindex.b]
    tindex.b <- match(Y.b,unique(Y.b))
    val.b <- data[vindex.b,,drop=FALSE]
    ## browser()
    train.b <- data[ResampleIndex[,b],,drop=FALSE]
    ## if (b==1) print(train.b$days)
    ## if (b==1) print(val.b$days)
    NV=sum(vindex.b)                    # NROW(val.b)
    # }}}
    # {{{ IPCW
    ## if (ipcw.refit==TRUE){
    ## ipcw.call.b.i <- ipcw.call$weight.i
    ## ipcw.call.b.j <- ipcw.call$weight.j
    ## ipcw.call.b.i$data <- val.b
    ## ipcw.call.b.j$data <- val.b
    ## ipcw.call.b.i$subjectTimes <- Y.b
    ## ipcw.call.b.j$subjectTimes <- Y.b
    ## ipcw.b.i <- do.call("ipcw",ipcw.call.b.i)$IPCW.subjectTimes
    ## ipcw.b.j <- do.call("ipcw",ipcw.call.b.j)$IPCW.times
    ## }
    ## else{
    ipcw.b.i <- weights$weight.i[vindex.b]
    if (is.null(dim(weights$weight.j))){
      ipcw.b.j <- weights$weight.j
    }
    else{
      ipcw.b.j <- weights$weight.j[vindex.b,]
    }
    ## }
    # }}}
    # {{{ Building the models in training data
    if (!is.null(seed)) {
      set.seed(seed)
      ## message("seed:",seed)
    }
    trainModels <- lapply(1:NF,function(f){
      fit.b <- internalReevalFit(object=object[[f]],
                                 data=train.b,
                                 step=b,
                                 silent=FALSE,
                                 verbose=verbose)
      ## fit.b$call <- object[[f]]$call
      fit.b
    })
    # }}}
    # {{{ Saving the models?
    if (!is.null(savePath)){
      nix <- lapply(1:NF,function(f){
        fit.b <- trainModels[[f]]
        ## print(object.size(fit.b))
        fit.b$formula <- NULL
        ## print(environment(fit.b$formula))
        save(fit.b,file=paste(paste(savePath,"/",names(object)[f],"-bootstrap-",b,sep=""),".rda",sep=""))
      })
    }
    # }}}
    # {{{ Extracting parameters?
    if (!is.null(getFromModel)){
      ModelParameters <- lapply(1:NF,function(f){
        getParms <- getFromModel[[f]]
        if (is.null(getParms)) trainModels[[f]][getParms] else NULL
      })
    }
    # }}}
    # {{{ Check fits
    fitFailed <- lapply(trainModels,function(fit.b) (is.null(fit.b)))
    # }}}
    # {{{ Predicting the validation data
    predVal <- lapply(1:NF,function(f){
      fit.b <- trainModels[[f]]
      extraArgs <- giveToModel[[f]]
      if (predictHandlerFun %in% c("predictEventProb","predictLifeYearsLost")){
        try2predict <- try(pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=val.b,times=pred.times,cause=cause),extraArgs)))
      }
      else{
        try2predict <- try(pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=val.b,times=pred.times),extraArgs)))
      }
      ## browser()
      ## print(pred.b[1:5])
      if (inherits(try2predict,"try-error")==TRUE){
        if (verbose==TRUE) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
        NULL}
      else{
        pred.b
      }
    })
    # }}}
    # {{{ Compute cindex  for step b
    if (multiSplitTest==TRUE){
      stop("not yet defined: residual test for cindex")
      Residuals <- lapply(predVal,function(pred.b){
        if (is.null(pred.b))
          NA
        else{
          if (predictHandlerFun %in% c("predictEventProb","predictLifeYearsLost")){
            1
            ## matrix(.C("pecResidualsCR",pec=double(NT),resid=double(NT*NV),as.double(Y[vindex.b]),as.double(status[vindex.b]),as.double(event[vindex.b]),as.double(times),as.double(pred.b),as.double(ipcwTimes.b),as.double(IPCW.subjectTimes.b),as.integer(NV),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$resid,ncol=NT,byrow=FALSE)
          }
          else{
            1
            ## matrix(.C("pecResiduals",pec=double(NT),resid=double(NT*NV),as.double(Y[vindex.b]),as.double(status[vindex.b]),as.double(times),as.double(pred.b),as.double(ipcwTimes.b),as.double(IPCW.subjectTimes.b),as.integer(NV),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$resid,ncol=NT,byrow=FALSE)
          }
        }
      })
      names(Residuals) <- names(object)
      ## PredCindexStepB=lapply(Residuals,function(x){colMeans(x)})
      PredCindexStepB=1
    }
    else{
      PredCindexStepB <- lapply(predVal,function(pred.b){
        if (is.null(pred.b))
          NA
        else{
          if (predictHandlerFun %in% c("predictEventProb","predictLifeYearsLost")){
            Step.b.CindexResult <- .C("ccr",cindex=double(NT),concA=double(NT),pairsA=double(NT),concB=double(NT),pairsB=double(NT),as.integer(tindex.b),as.double(Y.b),as.integer(status[vindex.b]),as.integer(event[vindex.b]),as.double(eval.times),as.double(ipcw.b.i),as.double(ipcw.b.j),as.double(pred.b),as.integer(sum(vindex.b)),as.integer(NT),as.integer(tiedPredictionsIn),as.integer(tiedOutcomeIn),as.integer(tiedMatchIn),as.integer(!is.null(dim(ipcw.b.j))),NAOK=TRUE,package="pec")
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
                            as.integer(tindex.b),
                            as.double(Y.b),
                            as.integer(status[vindex.b]),
                            as.double(eval.times),
                            as.double(ipcw.b.i),
                            as.double(ipcw.b.j),
                            as.double(pred.b),
                            as.integer(sum(vindex.b)),
                            as.integer(NT),
                            as.integer(tiedPredictionsIn),
                            as.integer(tiedOutcomeIn),
                            as.integer(tiedMatchIn),
                            as.integer(!is.null(dim(ipcw.b.j))),
                            NAOK=TRUE,
                            package="pec")            
            Cindex.b <- cindexOut$cindex
            Pairs.b <- cindexOut$pairs 
            Concordant.b <- cindexOut$conc
            list(Cindex.b=Cindex.b,Pairs.b=Pairs.b,Concordant.b=Concordant.b)
          }
        }
      })
    }
    # }}}
    # {{{ van de Wiel's test
    ## if (multiSplitTest==TRUE){
    ## testedResid <- testResiduals(Residuals,times=times,testTimes=testTimes,rangeInt=testIBS,confInt=confInt,confLevel=confLevel)
    ## }
    # }}}
    # {{{ looping output
    ##     if (multiSplitTest==TRUE)
    ##       loopOut=list(PredCindexStepB=PredCindexStepB,testedResid=testedResid)
    ##     else
    loopOut=list(PredCindexStepB=PredCindexStepB)
    ##     if (keepResiduals==TRUE)  
    ##       loopOut=c(loopOut,list(Residuals=lapply(Residuals,function(R){
    ##         R[,prodlim::sindex(eval.times=testTimes,jump.times=times)]
    ##       })))

    if (!is.null(getFromModel)){
      loopOut=c(loopOut,list(ModelParameters=ModelParameters))
    }
    loopOut
  }
  ## })
  b <- 1
  ## if (require(foreach)){
  if (missing(slaveseed)||is.null(slaveseed))
    slaveseed <- sample(1:1000000,size=B,replace=FALSE)
  Looping <- foreach::foreach (b= 1:B) %dopar% step(b,slaveseed[[b]])
  ## }
  ## else{
  ## Looping <- lapply(1:B,function(b){step(b,seed=NULL)})
  ## }
  # }}}
  # {{{ output
  ## 
  ## 
  ##    1. a list of NF matrices each with B (rows) and NT columns
  ##       the prediction error curves
  ## 
  ## if (verbose==TRUE && B>1) cat("\n")
  BootstrapCrossValCindexMat <- lapply(1:NF,function(f){
    ## matrix with NT columns and b rows
    do.call("rbind",lapply(Looping,function(b){
      c.b <- b$PredCindexStepB[[f]]$Cindex.b
      c.b
      ## pairs.b <- b$PredCindexStepB[[f]]$Pairs.b
      ## conc.b <- b$PredCindexStepB[[f]]$Concordant.b
    }))
  })
  ## 
  ##    2. a list of NF average out-of-bag prediction error curves
  ##       with length NT
  ##
  BootstrapCrossValCindex <- lapply(BootstrapCrossValCindexMat,colMeans)
  out <- list(BootstrapCrossValCindex=BootstrapCrossValCindex)
  ## 
  ##   3. the results of B residual tests 
  ##
  ##   print(str(Looping))
  if (multiSplitTest==TRUE){
    out$testedResid <- lapply(Looping,function(x)x$testedResid)
  }
  ## 
  ##   4. model parameters
  ##
  if (!is.null(getFromModel)){
    out$ModelParameters <- lapply(1:NF,function(f){
      lapply(Looping,function(x)x$ModelParameters[[f]])
    })
  }
  ## 
  ##   5. bootstrap crossvalidation results
  ##
  if (keepMatrix==TRUE)
    out$BootstrapCrossValCindexMat <- BootstrapCrossValCindexMat
  ## 
  ##   6. residuals
  ##
  ##   if (keepResiduals==TRUE){
  ##     out$Residuals <- lapply(1:NF,function(f){
  ##       bootResiduals <- lapply(Looping,function(b){
  ##         b$Residuals[[f]]
  ##       })
  ##       names(bootResiduals) <- paste("testSample",1:B,sep=".")
  ##       bootResiduals
  ##     })
  ##     names(out$Residuals) <- names(object)
  ##   }
  out
  # }}}
}

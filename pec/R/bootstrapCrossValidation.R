bootstrapCrossValidation <- function(object,
                                     data,
                                     Y,
                                     status,
                                     event,
                                     times,
                                     cause,
                                     ipcw,
                                     ipcw.refit=FALSE,
                                     ipcw.call,
                                     splitMethod,
                                     multiSplitTest,
                                     keepResiduals,
                                     testIBS,
                                     testTimes,
                                     newdata,
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
  NT <- length(times)
  NF <- length(object) 
  ResampleIndex <- splitMethod$index
  # }}}
  step <- function(b,seed){
    if (verbose==TRUE) internalTalk(b,B)
    # {{{ training and validation data
    vindex.b <- match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0
    val.b <- data[vindex.b,,drop=FALSE]
    train.b <- data[ResampleIndex[,b],,drop=FALSE]
    ## print(c(NROW(train.b), NROW(val.b)))
    NV=sum(vindex.b)                    # NROW(val.b)
    # }}}
    # {{{ IPCW
    if (ipcw.refit==TRUE){
      ipcw.call.b <- ipcw.call
      ipcw.call.b$data <- val.b
      ipcw.call.b$subjectTimes <- Y[vindex.b]
      ipcw.b <- do.call("ipcw",ipcw.call.b)
      ipcwTimes.b <- ipcw.b$IPCW.times
      IPCW.subjectTimes.b <- ipcw.b$IPCW.subjectTimes
    }
    else{
      IPCW.subjectTimes.b <- ipcw$IPCW.subjectTimes[vindex.b]
      if (ipcw$dim==1)
        ipcwTimes.b <- ipcw$IPCW.times[vindex.b,]
      else
        ipcwTimes.b <- ipcw$IPCW.times
    }
    
    # }}}
    # {{{ Building the models in training data
    if (!is.null(seed)) {
        set.seed(seed)
        ## if (verbose) message("seed:",seed)
    }
    trainModels <- lapply(1:NF,function(f){
        fit.b <- internalReevalFit(object=object[[f]],data=train.b,step=b,silent=FALSE,verbose=verbose)
        ## this was a good idea to reduce the memory usage:
        ## fit.b$call <- object[[f]]$call
        ## fit.b$call <- NULL
        ## however, it does not work with the new version of the survival package
        ## in which the survfit.coxph function checks the response 'y'
        ## next try
        ## print("before")
        ## print(object.size(fit.b))
        ## print("after")
        ## browser()
        ## fit.b$call$data <- substitute(train.b)
        ## print(object.size(fit.b))
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
            print(trainModels[[f]][getParms])
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
        if (predictHandlerFun == "predictEventProb"){
            try2predict <- try(pred.b <- do.call(predictHandlerFun,
                                                 c(list(object=fit.b,newdata=val.b,times=times,cause=cause),extraArgs)))
        }
        else{
            try2predict <- try(pred.b <- do.call(predictHandlerFun,
                                                 c(list(object=fit.b,newdata=val.b,times=times),extraArgs)))
        }
        if (inherits(try2predict,"try-error")==TRUE){
            if (verbose==TRUE) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
            NULL}
        else{
            pred.b
        }
    })
    # }}}
    # {{{ Compute prediction error curves for step b
    if (multiSplitTest==TRUE){
        Residuals <- lapply(predVal,function(pred.b){
            if (is.null(pred.b))
                NA
            else{
                if (predictHandlerFun == "predictEventProb"){
                    matrix(.C("pecResidualsCR",
                              pec=double(NT),
                              resid=double(NT*NV),
                              as.double(Y[vindex.b]),
                              as.double(status[vindex.b]),
                              as.double(event[vindex.b]),
                              as.double(times),
                              as.double(pred.b),
                              as.double(ipcwTimes.b),
                              as.double(IPCW.subjectTimes.b),
                              as.integer(NV),
                              as.integer(NT),
                              as.integer(ipcw$dim),
                              as.integer(is.null(dim(pred.b))),
                              NAOK=TRUE,
                              PACKAGE="pec")$resid,ncol=NT,byrow=FALSE)
                }
                else{
                    matrix(.C("pecResiduals",
                              pec=double(NT),
                              resid=double(NT*NV),
                              as.double(Y[vindex.b]),
                              as.double(status[vindex.b]),
                              as.double(times),
                              as.double(pred.b),
                              as.double(ipcwTimes.b),
                              as.double(IPCW.subjectTimes.b),
                              as.integer(NV),
                              as.integer(NT),
                              as.integer(ipcw$dim),
                              as.integer(is.null(dim(pred.b))),
                              NAOK=TRUE,
                              PACKAGE="pec")$resid,ncol=NT,byrow=FALSE)
                }
            }
        })
        names(Residuals) <- names(object)
        PredErrStepB=lapply(Residuals,function(x){colMeans(x)})
    }
    else{
        PredErrStepB <- lapply(predVal,function(pred.b){
            if (is.null(pred.b))
                NA
            else{
                if (predictHandlerFun=="predictEventProb")
                    pecOut <- .C("pecCR",
                                 pec=double(NT),
                                 as.double(Y[vindex.b]),
                                 as.double(status[vindex.b]),
                                 as.double(event[vindex.b]),
                                 as.double(times),
                                 as.double(pred.b),
                                 as.double(ipcwTimes.b),
                                 as.double(IPCW.subjectTimes.b),
                                 as.integer(NV),
                                 as.integer(NT),
                                 as.integer(ipcw$dim),
                                 as.integer(is.null(dim(pred.b))),
                                 NAOK=TRUE,
                                 PACKAGE="pec")$pec
                else
                    pecOut <- .C("pec",pec=double(NT),as.double(Y[vindex.b]),as.double(status[vindex.b]),as.double(times),as.double(pred.b),as.double(ipcwTimes.b),as.double(IPCW.subjectTimes.b),as.integer(NV),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$pec
            }
        })
    }
    # }}}
    # {{{ van de Wiel's test
    if (multiSplitTest==TRUE){
        testedResid <- testResiduals(Residuals,times=times,testTimes=testTimes,rangeInt=testIBS,confInt=confInt,confLevel=confLevel)
    }
    # }}}
    # {{{ looping output 
    if (multiSplitTest==TRUE)
        loopOut=list(PredErrStepB=PredErrStepB,testedResid=testedResid)
    else
        loopOut=list(PredErrStepB=PredErrStepB)
    if (keepResiduals==TRUE)  
        loopOut=c(loopOut,list(Residuals=lapply(Residuals,function(R){
            R[,prodlim::sindex(eval.times=testTimes,jump.times=times)]
        })))
    if (!is.null(getFromModel)){
        loopOut=c(loopOut,list(ModelParameters=ModelParameters))
    }
    loopOut
}
  b <- 1
  ## if (require(foreach)){
  if (missing(slaveseed)||is.null(slaveseed))
      slaveseed <- sample(1:1000000,size=B,replace=FALSE)
  
  Looping <- foreach::foreach (b= 1:B)  %dopar% step(b,slaveseed[[b]])
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
  if (verbose==TRUE && B>1) cat("\n")
  BootstrapCrossValErrMat <- lapply(1:NF,function(f){
      ## matrix with NT columns and b rows
      do.call("rbind",lapply(Looping,function(b){
          b$PredErrStepB[[f]]
      }))
  })
  ## 
  ##    2. a list of NF average out-of-bag prediction error curves
  ##       with length NT
  ## 
  BootstrapCrossValErr <- lapply(BootstrapCrossValErrMat,colMeans)
  ##   function(x){
  ##     if (na.accept>0) colMeans(x,na.rm=sum(is.na(b))<na.accept)
  ##     else
  ##     colMeans(x)
  ##   })
  out <- list(BootstrapCrossValErr=BootstrapCrossValErr)
  ## 
  ##   3. the results of B residual tests 
  ##   
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
      out$BootstrapCrossValErrMat <- BootstrapCrossValErrMat
  ## 
  ##   6. residuals
  ##
  if (keepResiduals==TRUE){
      out$Residuals <- lapply(1:NF,function(f){
          bootResiduals <- lapply(Looping,function(b){
              b$Residuals[[f]]
          })
          names(bootResiduals) <- paste("testSample",1:B,sep=".")
          bootResiduals
      })
      names(out$Residuals) <- names(object)
  }
  out
  # }}}
}

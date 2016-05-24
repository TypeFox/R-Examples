pseudoPec <- function(object,
                      formula,
                      data,
                      times,
                      cause,
                      start,
                      maxtime,
                      exact=TRUE,
                      exactness=100,
                      fillChar=NA,
                      cens.model="cox",
                      ipcw.refit=FALSE,
                      splitMethod="none",
                      B,
                      M,
                      reference=TRUE,
                      model.args=NULL,
                      model.parms=NULL,
                      keep.index=FALSE,
                      keep.matrix=FALSE,
                      keep.models=FALSE,
                      keep.residuals=FALSE,
                      keep.pvalues=FALSE,
                      noinf.permute=FALSE,
                      multiSplitTest=FALSE,
                      testIBS,
                      testTimes,
                      confInt=FALSE,
                      confLevel=0.95,
                      verbose=TRUE,
                      savePath=NULL,
                      ...){
  
  # {{{ checking integrity some arguments
  theCall=match.call()
  if (match("replan",names(theCall),nomatch=FALSE))
    stop("Argument name 'replan' has been replaced by 'splitMethod'.")
  
  if (!missing(testIBS) && (!(is.logical(testIBS) || (length(testIBS)==2 && is.numeric(testIBS)))))
      stop("Argument testIBS can be TRUE/FALSE or a vector of two numeric values.")
  if (missing(testIBS)) testIBS <- FALSE
  if (keep.residuals && missing(testTimes))
      stop("To keep.residuals please spseudoPecify testTimes.")
  if (missing(splitMethod) && multiSplitTest==TRUE){
      stop("Need data splitting to compute van de Wiel's test")
  }
  if (missing(M) && multiSplitTest)
      M <- NA
  # }}}
  # {{{ check and convert object
  if (class(object)[1]!="list") {
      object <- list(object)
  }
  # }}}
  # {{{ formula
  if (missing(formula)){
      if (length(grep("~",as.character(object[[1]]$call$formula)))==0){
          stop(paste("Argument formula is missing and first model has no usable formula:",as.character(object[[1]]$call$formula)))
      } else{
          ftry <- try(formula <- eval(object[[1]]$call$formula),silent=TRUE)
          if ((class(ftry)=="try-error") || match("formula",class(formula),nomatch=0)==0)
              stop("Argument formula is missing and first model has no usable formula.")
          else if (verbose)
              warning("Formula missing. Using formula from first model")
      }
  }  
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[1]=="~")
      ||
      (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
    stop("Invalid spseudoPecification of formula. Perhaps forgotten right hand side?\nNote that any subsetting, ie data$var or data[,\"var\"], is invalid for this function.")
  }
  else{
    if (!(formula.names[2] %in% c("Surv","Hist")))
      survp <- FALSE
    else
      survp <- TRUE
  }
  # }}}
  # {{{ data
  if (missing(data)){
    data <- eval(object[[1]]$call$data)
    if (match("data.frame",class(data),nomatch=0)==0)
      stop("Argument data is missing.")
    else
      if (verbose)
        warning("Argument data is missing. I use the data from the call to the first model instead.")
  }
  # }}}
  # {{{ response
  m <- model.frame(formula,data,na.action=na.fail)
  response <- model.response(m)
  if (match("Surv",class(response),nomatch=0)!=0){
    attr(response,"model") <- "survival"
    attr(response,"cens.type") <- "rightCensored"
    model.type <- "survival"
  }
  model.type <- attr(response,"model")
  if (model.type=="competing.risks"){
    predictHandlerFun <- "predictEventProb"
    if (missing(cause))
      cause <- attr(response,"state")[1]
  }
  else{
    if (survp==FALSE && NCOL(response)!=1) stop("Response must be one-dimensional.")
    if (survp==TRUE && NCOL(response)!=2) stop("Survival response must have two columns: time and status.")
    predictHandlerFun <- "predictSurvProb"
  }
  # }}}
  # {{{ prediction models
  if (reference==TRUE) {
    ProdLimform <- update(formula,".~1")
    ProdLimfit <- prodlim::prodlim(ProdLimform,data)
    ProdLimfit$call$data <- NULL
    ProdLimfit$formula <- NULL
    ProdLimfit$call$formula=ProdLimform
    if (model.type=="competing.risks")
      object <- c(list(AalenJohannsen=ProdLimfit),object)
    else
      object <- c(list(KaplanMeier=ProdLimfit),object)
  }
  if (is.null(names(object))){
    names(object) <- sapply(object,function(o)class(o)[1])
  }
  else{
    names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
  }
  names(object) <- make.names(names(object),unique=TRUE)
  NF <- length(object) 

  # }}}  
  # {{{ sort the data 

  if (survp){
    neworder <- order(response[,"time"],-response[,"status"])
    if (predictHandlerFun=="predictEventProb"){    
      event <- prodlim::getEvent(response,mode="character")
      event <- event[neworder]
    }
    response <- response[neworder,,drop=FALSE]
    Y <- response[,"time"]
    status <- response[,"status"]
  }
  else{
    cens.model <- "none"
    neworder <- order(response)
    Y <- response[neworder]
    status <- rep(1,length(Y))
  }
  ## for competing risks find the cause of interest.
  if (predictHandlerFun=="predictEventProb"){
    availableCauses <- unique(event)
    if (!match(cause, availableCauses,nomatch=FALSE))
      stop("Cause ",cause," is not among the available causes: ",paste(availableCauses,collapse=", "))
    event <- event==cause
  }

  data <- data[neworder,]
  
  unique.Y <- unique(Y)
  N <- length(Y)
  NU <- length(unique.Y)

  # }}}
  # {{{ splitMethod
  ##   N <- NROW(data)
  splitMethod <- resolvesplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
  B <- splitMethod$B
  
  ResampleIndex <- splitMethod$index
  k <- splitMethod$k
  do.resample <- !(is.null(ResampleIndex))
  if (keep.matrix==TRUE & !do.resample){
    warning("Argument keep.matrix set to FALSE, since no resampling/crossvalidation is requested.")
    keep.matrix <- FALSE
  }

  # }}}      
  # {{{ find maxtime, start, and jumptimes in the range of the response 
  if (missing(maxtime) || is.null(maxtime))
    maxtime <- unique.Y[NU]

  if (missing(start))
    if (survp==TRUE)
      start <- 0  ## survival times are positive
    else
      start <- min(unique.Y) 
  
  if (missing(times)){
    if (exact==TRUE)
      times <- unique(c(start,unique.Y))
    else
      times <- seq(start,maxtime,(maxtime - start)/exactness)
  }
  else{
    if (exact==TRUE) 
      times <- sort(c(start,unique(times),unique.Y))
    else
      times <- sort(unique(c(start,times)))
  }
  times <- times[times<=maxtime]
  NT <-  length(times)
  
  # }}}
  # {{{ compute jackknife pseudo values
  if (predictHandlerFun=="predictEventProb"){
    ## aj <- prodlim(Hist(time,event)~1,data=data.frame(time=response[,"time"],status=response[,"event"]))
    YY <- prodlim::jackknife.competing.risks(ProdLimfit,times=times,cause=cause)
    ## FIXME:
    YY[,1] <- 0
  }else{
    ## km <- prodlim(Hist(time,status)~1,data=data.frame(time=response[,"time"],status=response[,"status"]))
    YY <- prodlim::jackknife.survival(ProdLimfit,times=times)
  }
  # }}}
  # {{{ checking the models for compatibility with resampling
  if (do.resample){
    cm <- checkModels(object=object,model.args=model.args,model.parms=model.parms,splitMethod=splitMethod$internal.name)
    model.args <- cm$model.args
    model.parms <- cm$model.parms
  }
  # }}}
  # {{{ ---------------------------Apparent error---------------------------
  AppErr <- lapply(1:NF,function(f){
    fit <- object[[f]]
    extraArgs <- model.args[[f]]

    if (predictHandlerFun=="predictEventProb"){
      pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times,cause=cause),extraArgs))
      if (class(fit)[[1]]%in% c("matrix","numeric")) pred <- pred[neworder,]
    }
    else{
      pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times),extraArgs))
      if (class(fit)[[1]]%in% c("matrix","numeric")) pred <- pred[neworder,]
    }
    if (NCOL(pred)==1) pred <- matrix(rep(pred,N),byrow=TRUE,ncol=NT,nrow=N)
    pec <- colMeans(pred^2) + colMeans((1-2*pred) * YY)
    pec
  })
  names(AppErr) <- names(object)
  # }}}
  # {{{------------------------No information error------------------------

  if (splitMethod$internal.name %in% c("Boot632plus")){
    if (noinf.permute==FALSE){
      NoInfErr <- lapply(1:NF,function(f){
        fit <- object[[f]]
        extraArgs <- model.args[[f]]
        pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times),extraArgs))
        extraArgs <- model.args[[f]]
        if (predictHandlerFun=="predictEventProb")
          .C("pseudoPec_noinfCR",pseudoPec=double(NT),as.double(YY),as.double(event),as.double(times),as.double(pred),as.integer(N),as.integer(NT),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pseudoPec
        else
          .C("pseudoPec_noinf",pseudoPec=double(NT),as.double(YY),as.double(times),as.double(pred),as.integer(N),as.integer(NT),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pseudoPec
      })
      names(NoInfErr) <- names(object)
    }else{
      NoInfErrList <- lapply(1:B,function(b){
        if (verbose==TRUE){
          internalTalk(b,B,sign=".")
        }
        responseNames <- colnames(response)
        noinf.b <- data[sample(1:NROW(data),replace=FALSE),-match(responseNames,names(data))]
        noinf.b[,responseNames] <- response
        noinfPredErr <- lapply(1:NF,function(f){
          fit.b <- internalReevalFit(object=object[[f]],data=noinf.b,step=b,silent=FALSE,verbose=verbose)
          ## fit.b$call <- object[[f]]$call
          extraArgs <- model.args[[f]]
          pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=noinf.b,times=times),extraArgs))
          if (predictHandlerFun=="predictEventProb")
            .C("pseudoPecCR",pseudoPec=double(NT),as.double(YY),as.double(event),as.double(times),as.double(pred.b),as.integer(N),as.integer(NT),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$pseudoPec
          else
            .C("pseudoPec",pseudoPec=double(NT),as.double(YY),as.double(times),as.double(pred.b),as.integer(N),as.integer(NT),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$pseudoPec
        })
        noinfPredErr
      })
      NoInfErrMat <- lapply(1:NF,function(f){
        do.call("rbind",lapply(NoInfErrList,function(x){
          x[[f]]
        }))})
      NoInfErr <- lapply(NoInfErrMat,colMeans)
      names(NoInfErr) <- names(object)
    }
  }

  # }}}
  # {{{--------------k-fold and leave-one-out CrossValidation-----------------------
  if (splitMethod$internal.name %in% c("crossval","loocv")){
    kCV <- pseudo.kFoldCrossValidation(object=object,data=data,Y=Y,status=status,event=event,times=times,cause=cause,splitMethod=splitMethod,giveToModel=model.args,predictHandlerFun=predictHandlerFun,keep=keep.matrix,verbose=verbose)
    CrossValErr <- kCV$CrossValErr
    if (keep.matrix && B>1)
      CrossValErrMat <- kCV$CrossValErrMat
  }
  # }}}
  # {{{ ----------------------BootstrapCrossValidation----------------------

  if (splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632")){
    if (missing(testTimes)){
      testTimes <- NULL
    }
    if (missing(testIBS)){
      testIBS <- NULL
    }
    BootCv <- pseudo.bootstrapCrossValidation(object=object,data=data,Y=Y,status=status,event=event,times=times,cause=cause,splitMethod=splitMethod,multiSplitTest=multiSplitTest,testIBS=testIBS,testTimes=testTimes,confInt=confInt,confLevel=confLevel,getFromModel=model.parms,giveToModel=model.args,predictHandlerFun=predictHandlerFun,keepMatrix=keep.matrix,keepResiduals=keep.residuals,verbose=verbose,savePath=savePath)
    BootstrapCrossValErr <- BootCv$BootstrapCrossValErr
    Residuals <- BootCv$Residuals
    names(BootstrapCrossValErr) <- names(object)
    if (multiSplitTest==TRUE){
      comparisons <- allComparisons(names(object))
      multiSplitTest <- list(testIBS=testIBS,B=B,M=M,N=N,testTimes=testTimes)
      multiSplitTest$Comparisons <- lapply(1:length(comparisons),function(cc){
        if (length(testTimes)>0){
          allPairwisePvaluesTimes <- do.call("rbind",lapply(BootCv$testedResid,function(b){
            b$pValue[[cc]]}))
          out <- list(pValueTimes=apply(allPairwisePvaluesTimes,2,median))
          if (keep.pvalues==TRUE){
            out$allPairwisePvaluesTimes <- allPairwisePvaluesTimes
          }
        }
        else out <- NULL
        if(length(testIBS)>0)
          allPairwisePvaluesIBS <- sapply(BootCv$testedResid,function(b){
            b$IBSpValue[[cc]]})
        out$pValueIBS <- median(allPairwisePvaluesIBS)
        if (keep.pvalues==TRUE){
          out$allPairwisePvaluesIBS <- allPairwisePvaluesIBS}
        out
      })
      names(multiSplitTest$Comparisons) <- names(comparisons)
      ## multiSplitTest$splitMethod <- splitMethod
      class(multiSplitTest) <- "multiSplitTest"
    }
    ## upperLimits <- lapply(BootCv$testedResid,function(x){x[,1:length(testTimes)]})
    ##     if (testIBS==TRUE){
    ##       wtestIBSpValues <- do.call("cbind",apply(BootCv$testedResid,function(x){x[,length(testTimes)+1]}))
    ##     }
    ## wtestIBSupper <- BootCv$testedResid$wtestIBSupper
    ##   }
    if (keep.matrix==TRUE){
      BootstrapCrossValErrMat <- BootCv$BootstrapCrossValErrMat
      names(BootstrapCrossValErr) <- names(object)
    }
  }

  # }}}
  # {{{ Bootstrap .632

  if (splitMethod$internal.name=="Boot632"){
    B632Err <- lapply(1:NF,function(f){
      B632Err <- .368 * AppErr[[f]] + .632 * BootstrapCrossValErr[[f]]
    })
    names(B632Err) <- names(object)
  }

  # }}}    
  # {{{ Bootstrap .632+
  if (splitMethod$internal.name=="Boot632plus"){
    B632plusErr <- lapply(1:NF,function(f){
      Err1 <- pmin(BootstrapCrossValErr[[f]],NoInfErr[[f]])
      overfit <- (Err1 - AppErr[[f]]) / (NoInfErr[[f]] - AppErr[[f]])
      overfit[!(Err1>AppErr[[f]])] <- 0
      w <- .632 / (1 - .368 * overfit)
      B632plusErr <- (1-w) * AppErr[[f]]  + w * Err1
      B632plusErr
      ## w[NoInfErr<=BootstrapCrossValErr] <- 1
      ## B632plus.error <- (1-w) * AppErr  + w * BootstrapCrossValErr
    })
    names(B632plusErr) <- names(object)
  }
  # }}}
  # {{{ prepare output

  out <- switch(splitMethod$internal.name,
                "noPlan"=list("AppErr"=AppErr),
                "Boot632plus"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr,"NoInfErr"=NoInfErr,"Boot632plusErr"=B632plusErr),
                "Boot632"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr,"Boot632Err"=B632Err),
                "BootCv"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr),
                "loocv"=list("AppErr"=AppErr,"loocvErr"=CrossValErr),
                "crossval"=list("AppErr"=AppErr,"crossvalErr"=CrossValErr),
                "noinf"=list("AppErr"=AppErr,"NoInfErr"=NoInfErr))
  observed.maxtime <- sapply(out,function(x){
    lapply(x,function(y){times[length(y)-sum(is.na(y))-1]})
  })
  minmaxtime <- min(unlist(observed.maxtime))
  if (multiSplitTest==TRUE){
    out <- c(out,list(multiSplitTest=multiSplitTest))
  }
  if (keep.residuals==TRUE){
    out <- c(out,list(Residuals=Residuals))
  }
  if (keep.matrix==TRUE && splitMethod$internal.name!="no"){
    ##       if (splitMethod$internal.name=="plain") out <- c(out,"BootErrMat"=BootErrMat)
    ##       else
    if (splitMethod$internal.name %in% c("crossval","loocv")){
      if (B>1)
        out <- c(out,list("CrossValErrMat"=CrossValErrMat))
    }
    else{
      if (splitMethod$internal.name!="noinf")
        out <- c(out,list("BootstrapCrossValErrMat"=BootstrapCrossValErrMat))
    }
  }
  if (!is.na(fillChar))
    out <- lapply(out,function(o){
      o[is.na(o)] <- fillChar
      o
    })
  if (!is.null(model.parms)) out <- c(out,list("ModelParameters"=BootCv$ModelParameters))
  ## if (na.accept>0) out <- c(out,list("failed"=failed))
  out
  n.risk <- N - prodlim::sindex(Y,times)
  if (!keep.index) splitMethod$index <- NULL

  # }}}    
  # {{{ put out
  
  if(keep.models==TRUE){
      outmodels <- object
  } else{
        outmodels <- names(object)
        names(outmodels) <- names(object)
    }
  out <- c(out,
           list(call=theCall,
                time=times,
                n.risk=n.risk,
                models=outmodels,
                maxtime=maxtime,
                observed.maxtime=observed.maxtime,
                minmaxtime=minmaxtime,
                reference=reference,
                start=min(times),
                cens.model=cens.model,
                exact=exact,
                cens.model="pseudovalues",
                splitMethod=splitMethod))
  class(out) <- "pec"
  out

  # }}}
}



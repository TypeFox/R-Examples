ConfScore <- function(object,
                      data,
                      newdata,
                      splitMethod="bootcv",
                      B,
                      M,
                      verbose=TRUE){
  # {{{ Find splits
  stopifnot(B>.368*NROW(data))
  splitMethod <- MgSplitMethods(splitMethod=splitMethod,B=B,N=NROW(data),M=M)
  ResampleIndex <- splitMethod$index
  # }}}
  # {{{ Who to predict
  if (missing(newdata)){
    predTestSet <- TRUE
    N <- NROW(data)
    
  }
  else{
    predTestSet <- FALSE
    N <- NROW(newdata)
  }
  # }}}
  # {{{ fit models in a loop over training sets  
  predictedValues <- lapply(1:B,function(b){
    MgTalk(b,B)
    d <- data[ResampleIndex[,b],,drop=FALSE]
    predList <- lapply(object,function(model){
      fit.b <- MgRefit(model,data=d,silent=FALSE)
      fit.b$call <- NULL
      gc()
      if (predTestSet==TRUE){
        vindex.b <- match(1:NROW(data),ResampleIndex[,b],nomatch=0)==0
        newdata <- data[vindex.b,,drop=FALSE]
      }
      try2predict <- try(pred.b <- do.call("predictStatusProb",
                                           c(list(object=fit.b,newdata=newdata))))
      if (inherits(try2predict,"try-error")==TRUE){
        if (verbose==TRUE) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
        NULL}
      else{
        if (predTestSet==TRUE){
          names(pred.b) <- (1:NROW(data))[vindex.b]
        }
        pred.b
      }
    })
    predList
  })

  # }}}
  # {{{ collect the predictions
  ##
  ## the output will have an entry for each model 
  ## each contains three things:
  ##
  ##     1. score: a vector with confidence scores (1-sd) for newdata
  ##
  ##     2. meanPred: a vector with mean predictions for newdata
  ## 
  ##     3. B: only if newdata is missing, for each subject in data,
  ##           the number of times where the subject was not in the training set
  ##
  out <- list(splitMethod=splitMethod,predTestSet=predTestSet)
  out$models <- lapply(1:length(object),function(m){
    predSubject <- lapply(1:N,function(i){
      if (predTestSet==TRUE){
        inTest <- is.na(apply(ResampleIndex,2,function(index){match(i,unique(index),nomatch=NA)}))
      }
      else
        inTest <- rep(TRUE,B)
      predTrain <- sapply(predictedValues[inTest],function(x){
        x[[m]][as.character(i)]
      })
      predTrain
    })
    ## browser()
    score <- sapply(predSubject,function(i){1-sd(i)})
    meanPred <- sapply(predSubject,function(i){mean(i)})
    B <- sapply(predSubject,function(i){length(i)})
    if (predTestSet==TRUE)
      origPred <- predictStatusProb(object[[m]],newdata=data)
    else
      origPred <- predictStatusProb(object[[m]],newdata=newdata)
    list(B=B,score=score,meanPred=meanPred,origPred=origPred)
  })
  # }}}
  gc()
  cat("\n")
  names(out$models) <- names(object)
  class(out) <- "confScore"
  out
}

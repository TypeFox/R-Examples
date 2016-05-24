ConfScoreSurv <- function(object,
                          data,
                          newdata,
                          times,
                          splitMethod="BootCv",
                          B,
                          M,
                          verbose=TRUE){
  NF <- length(object)
  NT <- length(times)
  N <- NROW(data)
  # {{{ Find splits

  ## require(pec)
  splitMethod <- resolvesplitMethod(splitMethod=splitMethod,B=B,N=NROW(data),M=M)
  ResampleIndex <- splitMethod$index

  # }}}
  # {{{ Who to predict
  if (missing(newdata)){
    predTestSet <- TRUE
    NTEST <- NROW(data)
    stopifnot(M<N)
  }
  else{
    predTestSet <- FALSE
    NTEST <- NROW(newdata)
    M <- NTEST
  }
  # }}}
  # {{{ fit models in a loop over training sets
  ##   predictedValues <- lapply(1:B,function(b){
  predictedValues <- lapply(names(object),function(m){
    array(dim=c(B,NT,N))
  })
  names(predictedValues) <- names(object)
  for(b in 1:B){
    internalTalk(b,B)
    d <- data[ResampleIndex[,b],,drop=FALSE]
    for (j in 1:NF){
      model <- object[[j]]
      fit.b <- internalReevalFit(model,data=d,silent=FALSE)
      ## fit.b$call <- NULL
      ## gc()
      if (predTestSet==TRUE){
        vindex.b <- match(1:NTEST,unique(ResampleIndex[,b]),nomatch=0)==0
        newdata <- data[vindex.b,,drop=FALSE]
      }
      try2predict <- try(pred.b <- predictSurvProb(object=fit.b,
                                                   newdata=newdata,
                                                   times=times))
      if (inherits(try2predict,"try-error")==TRUE){
        if (verbose==TRUE)
          warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
      } else{
        if (predTestSet==TRUE){
          for (i in 1:NROW(pred.b)){
            who <- (1:NTEST)[vindex.b]
            predictedValues[[j]][b,,who[i]] <- pred.b[i,]
          }
        }
        else{
          for (i in (1:NTEST))
            predictedValues[[j]][b,,i] <- pred.b[i,]
        }
      }
    }
    rm(d,fit.b,pred.b)
    gc()
  }
  # }}}
  # {{{ collect the predictions
  out <- list(times=times,splitMethod=splitMethod,predTestSet=predTestSet)
  out$models <- lapply(1:NF,function(j){
    pmodel <- predictedValues[[j]]
    score <- t(apply(pmodel,3,function(x){
      1-apply(x,2,sd,na.rm=TRUE)
    }))
    meanPred <- t(apply(pmodel,3,function(x){
      colMeans(x,na.rm=TRUE)
    }))
    B <- apply(pmodel,3,function(x){
      sum(!is.na(x[,1]))
    })
    list(B=B,score=score,meanPred=meanPred)
  })
  # {{{ return a confScoreSurv object
  gc()
  cat("\n")
  names(out$models) <- names(object)
  class(out) <- "confScoreSurv"
  out
  # }}}
}
  
    ##
    ## the output will have an entry for each model 
    ## each contains three things:
    ##
    ##     1. score: a matrix with confidence scores (1-sd), each row corresponds
    ##               to a subject in newdata and each column to a time point
  ##
    ##     2. meanPred: a matrix of average training predictions, each row corresponds
    ##                  to a subject in newdata and each column to a time point
    ## 
    ##     3. B: only if newdata is missing, this for each subject in data,
    ##           the number of times where the subject was not in the training set
    ##
    ##     out <- list(times=times,splitMethod=splitMethod,predTestSet=predTestSet)
    ##       out$models <- lapply(1:length(object),function(m){
    ##         predSubject <- lapply(1:NTEST,function(i){
    ##           if (predTestSet==TRUE){
    ##             inTest <- is.na(apply(ResampleIndex,2,function(index){match(i,unique(index),nomatch=NA)}))
    ##           }
    ##           else
    ##             inTest <- rep(TRUE,B)
    ##           predTrain <- do.call("rbind",lapply(predictedValues[inTest],function(x){
    ##             x[[m]][as.character(i),]
    ##           }))
    ##           predTrain
    ##         })
    ##         score <- do.call("rbind",lapply(predSubject,function(i){
    ##           if (is.null(dim(i)))
    ##             rep(NA,length(times))
  ##           else
    ##             1-apply(i,2,sd)
    ##         }))
    ##         meanPrediction <- do.call("rbind",lapply(predSubject,function(i){
    ##           if (is.null(dim(i)))
    ##             rep(NA,length(times))
    ##           else
    ##             colMeans(i)
    ##         }))
    ##         B <- sapply(predSubject,function(i){
    ##           NROW(i)
    ##         })
    ##         list(B=B,score=score,meanPred=meanPrediction)
    ##       })
    # }}}

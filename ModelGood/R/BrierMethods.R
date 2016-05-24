# {{{ UseMethod
##' @export
Brier <- function(object,...){
  UseMethod("Brier",object=object)
}
# }}}
# {{{ Brier.default, Brier.glm, etc

#' @method Brier default
#' @S3method Brier default
Brier.default <- function(object,y,cbRatio=1,...){
  N <- length(object)
  if(length(y)!=N) stop("Arguments must have the same length")
  if(length(unique(y))!=2) stop("y must be binary")
  Y <- as.integer(as.character(factor(y,labels=c("0","1"))))
  resid <- Y-object
  resid[Y==1] <- cbRatio[1] * resid[Y==1]
  bs <- mean(resid^2)
  bs
}

#' @method Roc formula
#' @S3method Roc formula
Brier.formula <- function(object,formula,data,...){
    Brier.list(object=list("logistic.regression"=object),formula=object,data,...)
}

#' @S3method Brier glm
Brier.glm <- function(object,formula,data,...){
    stopifnot(object$family$family=="binomial")
    Brier.list(object=list(object),formula,data,...)
}

#' @method Brier lrm
#' @S3method Brier lrm
Brier.lrm <- function(object,formula,data,...){
    Brier.list(object=list(object),formula,data,...)
}

#' @method Brier ElasticNet
#' @S3method Brier ElasticNet
Brier.ElasticNet <- function(object,formula,data,...){
  Brier.list(list(object),formula,data,...)
}

#' @method Brier rpart
#' @S3method Brier rpart
Brier.rpart <- function(object,formula,data,...){
    Brier.list(object=list(object),formula,data,...)
}

#' @method Brier randomForest
#' @S3method Brier randomForest
Brier.randomForest <- function(object,formula,data,...){
    Brier.list(object=list(object),formula,data,...)
}

# }}}
# {{{ Brier.list
#' @method Brier list
#' @S3method Brier list
Brier.list <- function(object,
                       formula,
                       data,
                       splitMethod="noSplitMethod",
                       noinf.method=c("simulate"),
                       simulate="reeval",
                       cbRatio=1,
                       B,
                       M,
                       model.args=NULL,
                       model.parms=NULL,
                       keepModels=FALSE,
                       keepSampleIndex=FALSE,
                       keepCrossValRes=FALSE,
                       na.accept=0,
                       verbose=FALSE,
                       ...){
    # }}}
    theCall=match.call()
    if (match("replan",names(theCall),nomatch=FALSE))
        stop("Argument name 'replan' has been replaced by 'splitMethod'.")
    # {{{ models
    NF <- length(object) 
    if (is.null(names(object)))names(object) <- sapply(object,function(o)class(o)[1])
    else{names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])}
    names(object) <- make.names(names(object),unique=TRUE)
    # }}}
    # {{{ formula
    if (missing(formula)){
        formula <- eval(object[[1]]$call$formula)
        if (class(formula)!="formula")
            stop("Argument formula is missing.")
        else if (verbose)
            warning("Argument formula is missing. Try to use the formula from the call to the first model instead.")
    }
    # }}}
    # {{{ data
    if (missing(data)){
        data <- eval(object[[1]]$call$data)
        if (class(data)!="data.frame")
            stop("Argument data is missing.")
        else  if (verbose)
            warning("Argument data is missing. I have (ab)used the data from the call\n of the first model instead.")
    }
    # }}}
    # {{{ response
    m <- model.frame(formula,data,na.action=na.fail)
    Y <- model.response(m)
    if (is.factor(Y) && (length(levels(Y))==2) || length(unique(Y))==2) {
        Y <- factor(Y)
        Y <- as.numeric(Y==levels(Y)[2])
    }
    N <- length(Y)
    # }}}
    # {{{ SplitMethod
    SplitMethod <- MgSplitMethods(splitMethod=splitMethod,B=B,N=N,M=M,k=k)
    B <- SplitMethod$B
    CrossvalIndex <- SplitMethod$index
    if (!keepSampleIndex) SplitMethod$index <- NULL
    k <- SplitMethod$k
    do.crossval <- !(is.null(CrossvalIndex))
    ## if (missing(keepCrossValRes)) keepCrossValRes <- do.crossval
    # }}}
    # {{{ checking the models for compatibility with cross-validation
    if (do.crossval){
        cm <- MgCheck(object=object,model.args=model.args,model.parms=model.parms,SplitMethod=SplitMethod)
        model.args <- cm$model.args
        model.parms <- cm$model.parms
    }
    # }}}
    # {{{ computation of the Brier in a loop over the models 
    list.out <- lapply(1:NF,function(f){
        if (verbose==TRUE) cat("\n",names(object)[f],"\n")
        fit <- object[[f]]
        extract <- model.parms[[f]]
        # }}}
        # {{{ apparent Brier (use the same data for fitting and validation)
        pred <- do.call("predictStatusProb",c(list(object=fit,newdata=data),model.args[[f]]))
        AppBS <- Brier.default(object=pred,y=Y,cbRatio=cbRatio)
        # }}}
        # {{{ No information error  
        if (SplitMethod$internal.name %in% c("boot632plus","noinf")){
            if (noinf.method=="simulate"){
                if (verbose==TRUE)
                    cat("\nSimulate no information performance\n")
                NoInfBSList <- lapply(1:B,function(b){
                    if (verbose==TRUE) MgTalk(b,B)
                    data.b <- data
                    ## permute the response variable
                    responseName <- all.vars(formula)[1]
                    data.b[,responseName] <- sample(factor(Y),replace=FALSE)
                    fit.b <- MgRefit(object=fit,data=data.b,step=b,silent=na.accept>0,verbose=verbose)
                    if (is.null(fit.b)){
                        failed <- "fit"
                        NoInfBS <- NA
                    }
                    else{
                        try2predict <- try(pred.b <- do.call("predictStatusProb",c(list(object=fit.b,newdata=data.b),model.args[[f]])),silent=TRUE)
                        if (inherits(try2predict,"try-error")==TRUE){
                            failed <- "prediction"
                            NoInfBS <- NA
                        }
                        else{
                            failed <- NA
                            NoInfBS <- Brier.default(object=pred.b,y=data.b[,responseName],cbRatio=cbRatio)
                        }
                    }
                    NoInfBS
                })
                if (verbose==TRUE) cat("\n")

        # }}}
        # {{{ averaging the B noinf BS curves
        NoInfBS <- mean(unlist(NoInfBSList))
    }
      else{ ## constant
          NoInfBS <- .C("brier_noinf",bs=double(1),as.double(Y),as.double(pred),as.integer(N),NAOK=TRUE,PACKAGE="ModelGood")$bs
      }
  }
    if (SplitMethod$internal.name %in% c("boot632plus","bootcv","boot632")){
        # }}}
        # {{{ Bootcv aka BootstrapCrossValidation
        if (verbose==TRUE)
            cat("\nBootstrap cross-validation performance\n")
        compute.BootcvBSList <- lapply(1:B,function(b){
            if (verbose==TRUE) MgTalk(b,B)
            vindex.b <- match(1:N,CrossvalIndex[,b],nomatch=0)==0
            val.b <- data[vindex.b,,drop=FALSE]
            train.b <- data[CrossvalIndex[,b],,drop=FALSE]
            fit.b <- MgRefit(object=fit,data=train.b,step=b,silent=na.accept>0,verbose=verbose)
            if (!is.null(extract)) fit.parms <- fit.b[extract]
            else fit.parms <- NULL
            if (is.null(fit.b)){
                failed <- "fit"
                innerBootcvBS <- NA
            }
            else{
                try2predict <- try(pred.b <- do.call("predictStatusProb",c(list(object=fit.b,newdata=val.b),model.args[[f]])),silent=na.accept>0)
                if (inherits(try2predict,"try-error")==TRUE){
                    if (verbose==TRUE) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
                    failed <- "prediction"
                    innerBootcvBS <- NA
                }
                else{
                    failed <- NA
                    innerBootcvBS <- Brier.default(object=pred.b,y=Y[vindex.b],cbRatio=cbRatio)
                }
            }
            list("innerBootcvBS"=innerBootcvBS,"fit.parms"=fit.parms,"failed"=failed)
        })
        if (verbose==TRUE) cat("\n")
        if (!is.null(extract)) fitParms <- lapply(compute.BootcvBSList,function(x)x$fit.parms)
        failed <- na.omit(sapply(compute.BootcvBSList,function(x)x$failed))
        BootcvBSList <- lapply(compute.BootcvBSList,function(x)x$innerBootcvBS)
        # }}}
        # {{{ averaging the B bootcv Brier scores
        BCVBS <- mean(unlist(BootcvBSList))
    }
    # }}}
    # {{{ Bootstrap .632
    if (SplitMethod$internal.name=="boot632"){
        B632BS <- .368 * AppBS + .632 * BCVBS
    }
    # }}}
    # {{{ Bootstrap .632+
    if (SplitMethod$internal.name=="boot632plus"){
        R632PlusBS <- MgFormule632(App=AppBS,BCV=BCVBS,NoInf=NoInfBS,SmallerBetter=TRUE)
    }
    # }}}
    # {{{ output
    out <- switch(SplitMethod$internal.name,
                  "noSplitMethod"=list("BS"=AppBS),
                  "boot632"=list("AppBS"=AppBS,"BootcvBS"=BCVBS,"BS"=B632BS),
                  "boot632plus"=list("AppBS"=AppBS,"BootcvBS"=BCVBS,"NoInfBS"=NoInfBS,"weight"=R632PlusBS$weight,"overfit"=R632PlusBS$overfit,"BS"=R632PlusBS$B632Plus),
                  "bootcv"=list("AppBS"=AppBS,"BS"=BCVBS),
                  "noinf"=list("AppBS"=AppBS,"PredBS"=NoInfBS))
    if (keepCrossValRes==TRUE && SplitMethod$internal.name!="noSplitMethod"){
        if (SplitMethod$internal.name!="noinf")
            out <- c(out,list("BootcvBSList"=BootcvBSList))
    }
    if (!is.null(extract)) out <- c(out,list("fitParms"=fitParms))
    if (na.accept>0) out <- c(out,list("failed"=failed))
    out
})
  names.lout <- names(list.out[[1]])
  out <- lapply(names.lout,function(w){
      e <- lapply(list.out,function(x){x[[w]]})
      names(e) <- names(object)
      e
  })
  names(out) <- names.lout

  if(keepModels==TRUE)
    outmodels <- object
  else if (keepModels=="Call"){
    outmodels <- lapply(object,function(o)o$call)
    names(outmodels) <- names(object)
  }
  else{
    outmodels <- names(object)
    names(outmodels) <- names(object)
  }
  out <- c(out,
           list(call=match.call(),
                models=outmodels,
                method=SplitMethod))
  if (verbose==TRUE) cat("\n")
  ##   bs <- Brier(Y,P,N)$pred.error
  ##   pNull <- sum(Y)/N
  ##   nullbs <- sum((Y-pNull)^2)/N
  ## list(Response=table(Y),Brier=bs,R2=1-bs/nullbs)
  # }}}
  class(out) <- "Brier"
  out
}

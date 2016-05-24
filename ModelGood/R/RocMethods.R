# {{{ Roc.list
#' Comparing prediction models with Receiver operating characteristics and
#' Brier scores
#' 
#' Evaluation of the performance of risk prediction models with binary status
#' response variable (case/control or similar). Roc curves are either based on a
#' single continuous marker, or on the probability prediction of an event. 
#' Probability predictions are extracted from a given (statistical) model, such as logistic
#' regression, or algorithm, such as random forest. The area under the curve
#' and the Brier score is used to summarize and compare the performance.
#'
#' All functions work on a list of models to ease comparison. 
#' 
#' Bootstrap-crossvalidation techniques are implemented to estimate the
#' generalization performance of the model(s), i.e., the performance which can
#' be expected in new subjects.
#'
#' By default, when crossvalidation is involved, the ROC curve is
#' approximated on a grid of either sensitivities or specificities and
#' not computed at all unique changepoints of the crossvalidated ROC curves,
#' see Fawcett, T. (2006). The (density of the) grid can be controlled with the
#' argument: RocAverageGrid 
#' 
#' Missing data in the response or in the marker/predicted risk cause a failure.
#' 
#' For each R object which potentially can predict a probability for an event,
#' there should be a corresponding \code{predictStatusProb} method:
#' 
#' For example, to assess a prediction model which evaluates to a
#' \code{myclass} object one defines a function called
#' \code{predictStatusProb.myclass} with arguments
#' \code{object,newdata,...}. For example, the
#' function predictStatusProb.lrm looks like this:
#'
#' predictStatusProb.lrm <- function(object,newdata,...){
#'  p <- as.numeric(predict(object,newdata=newdata,type='fitted'))
#'  class(p) <- 'predictStatusProb'
#'  p
#'}
#' 
#' Currently implemented are \code{predictStatusProb} methods for the following
#' R-functions: \describe{
#' \item{}{\code{numeric} (marker values are passed on)}
#' \item{}{\code{formula} (single predictor: extracted from newdata and passed on,
#' multiple predictors: projected to score by logistic regression)}
#' \item{}{\code{glm} (from \code{library(stats)}}
#' \item{}{\code{lrm} (from \code{library(Design)}}
#' \item{}{\code{rpart} (from \code{library(rpart)})}
#' \item{}{\code{BinaryTree} (from \code{library(party)})}
#' \item{}{\code{ElasticNet} (a wrapper for glmnet from \code{library(glmnet)})}
#' \item{}{\code{randomForest} from \code{library(randomForest)}}
#' \item{}{\code{rfsrc} from \code{library(randomForestSRC)}}
#' }
#' 
#' @aliases Brier Brier.list Brier.glm Brier.lrm Brier.randomForest Brier.rpart
#' Roc Roc.list Roc.glm Roc.lrm Roc.randomForest Roc.rpart
#' @usage
#' \method{Roc}{list} (object, formula, data, splitMethod='noSplitMethod',
#' noinf.method=c('simulate'), simulate='reeval', B, M, breaks, cbRatio=1,
#' RocAverageMethod='vertical',
#' RocAverageGrid=switch(RocAverageMethod, 'vertical'=seq(0,1,.01),
#' 'horizontal'=seq(1,0,-.01)), model.args=NULL, model.parms=NULL,
#' keepModels=FALSE, keepSampleIndex=FALSE, keepCrossValRes=FALSE,
#' keepNoInfSimu, slaveseed, cores=1, na.accept=0, verbose=FALSE, ...)
#' 
#' @param object A named list of R objects that represent predictive markers, prediction models,
#' or prediction algorithms. The function \link{predictStatusProb} is called on the R objects
#' to extract the predicted risk (see details). 

#' For  cross-validation (e.g. when \code{splitMethod} is 'bootcv') all the
#' R objects in this list must include a \code{call} which can be evaluated in
#' a learning subset of the data.
#' @param formula A formula whose left hand side is used to identify the binary
#' outcome variable in \code{data}. If missing, use the formula of the (first) model
#' in object. 
#' @param data A data set in which to validate the prediction models.
#' If missing, the function tries to extract the data from the call of
#' the (first) model in object. 
#' 
#' The data set needs to have the same structure, variable names,
#' factor levels, etc., as the data in which the models
#' were trained. If the subjects in data were not used to train the
#' models given in \code{object}, this leads to an external validation
#' situation.
#' 
#' However, note that if one of the elements in \code{object} is a formula
#' then it is evaluated in this data set. 
#' 
#' @param splitMethod Method for estimating the generalization error.
#' 
#' \code{none}:Assess the models in the data given by \code{data}.
#' If this data set coincides with the train data where the models
#' were fitted this yields an apparent (or re-substitution) estimate
#' of performance. Otherwise, this leads to an external validation
#' situation.
#' 
#' \code{bootCV}: Internal bootstrap cross validation. The prediction models are trained
#' on \code{B} bootstrap samples of the \code{data}. Bootstrap samples
#' are either drawn with replacement from \code{data} (same size), or without
#' replacement of the size \code{M} where \code{M} is a number smaller
#' than \code{NROW(data)}. The model performance parameters (Roc, Brier, AUC)
#' are estimated with the observations that are NOT in the current
#' bootstrap sample.
#' 
#' \code{boot632}: Linear combination of the apparent performance and the
#' BootCV performance using the constant weight .632 (see Efron & Tibshirani, 1997).
#' 
#' \code{boot632plus}: Linear combination of apparent performance and Bootcv
#' using weights dependent on how the models perform in permuted data (see
#' Efron & Tibshirani, 1997).
#' 
#' \code{noinf}: Assess the models trained in permutations of \code{data}.
#' @param noinf.method Experimental: For .632+ method the way to obtain no-information performance. This can either be 'simulate' or 'none'.
#' @param simulate Experimental: For .632+ method. If \code{'reeval'} then the models are re-build
#' in the current permuted data for computing the no-information Roc curve.
#' @param B Number of repetitions for internal crossvalidation. The meaning depends on the argument
#' \code{splitMethod}: When \code{splitMethod in
#' c('Bootcv','Boot632','Boot632plus')} it is the number of bootstrap samples, default is 100.
#' Otherwise it is ignored.
#' @param M The size of the bootstrap samples for cross-validation without
#' replacement.
#' @param breaks Break points for computing the Roc curve. Defaults to
#' \code{seq(0,1,.01)} when crossvalidation is applied, i.e., when  \code{splitMethod in
#' c('Bootcv','Boot632','Boot632plus')}. Otherwise use all unique values of the
#' predictive marker.
#' 
#' @param cbRatio Experimental. Cost/benefit ratio. Default
#' value is 1, meaning that misclassified cases are as bad as misclassified
#' controls.
#' @param RocAverageMethod Method for averaging ROC curves across data splits.
#' If \code{'horizontal'} average crossvalidated specificities for fixed sensitivity values,
#' specified in \code{RocAverageGrid}, otherwise, if  \code{'vertical'},
#' average crossvalidated specificities for fixed sensitivity values.
#' See Fawcett, T. (2006) for details.
#' @param RocAverageGrid Grid points for the averaging of Roc curves.
#' A sequence of values at which to compute averages across the ROC curves
#' obtained for different data splits during crossvalidation.
#' @param model.args List of extra arguments that can be passed to the
#' \code{predictStatusProb} methods. The list must have an entry for each entry
#' in \code{object}.
#' @param model.parms List with exactly one entry for each entry in
#' \code{object}.  Each entry names parts of the value of the fitted models
#' that should be extracted and added to the output (see value).
#' @param keepModels If \code{FALSE} keep only the names of the elements of
#' object.  If \code{'Call'} then keep the call of the elements of object.
#' Else, add the object as it is to the output.
#' @param keepSampleIndex Logical. If \code{FALSE} remove the cross-validation
#' index (which tells who was in the learn and who in the validation set) from
#' the output list which otherwise is included in the method part of the output
#' list.
#' @param keepCrossValRes Logical. If \code{TRUE} add all \code{B}
#' crossvalidation results to the output (see value). Defaults to \code{TRUE}.
#' @param keepNoInfSimu Logical. If \code{TRUE} add the \code{B} results in
#' permuted data (for no-information performance) to the output (see value).
#' Defaults to \code{FALSE}.
#' @param slaveseed Vector of seeds, as long as \code{B}, to be given to the
#' slaves in parallel computing to control the models build in crossvalidation loop.
#' @param cores Number of cores for parallel computing.
#' Passed as the value of the argument \code{mc.cores}
#' when calling \code{\link{mclapply}}.
#' @param na.accept For 'Bootcv' estimate of performance. The
#' maximal number of bootstrap samples in which the training the models may fail
#' This should usually be a small number relative to \code{B}.
#' @param verbose if \code{TRUE} the procedure is reporting details of the
#' progress, e.g. it prints the current step in cross-validation procedures.
#' @param ... Used to pass arguments to submodules.
#' @return Object of class \code{Roc} or class \code{Brier}. 
#' 
#' Depending on the  \code{splitMethod} the object includes the following components:
#' 
#' \item{Roc, Brier, Auc}{A list of Roc curve(s), Brier scores (BS), and areas under the curves (Auc),
#' one for each element of argument \code{object},
#' estimated according to \code{splitMethod}.}
#' 
#' \item{weight}{The weight used to linear combine the \code{AppRoc} and
#' the \code{BootcvRoc} Only available if \code{splitMethod} is one of 'Boot632', or 'Boot632plus'.
#' }
#'
#' \item{overfit}{ Estimated \code{overfit} of the model(s). Only if
#' \code{splitMethod} is one of 'Boot632', or 'Boot632plus'.  }
#'
#' \item{call}{The call that produced the object}
#'
#' \item{models}{See keepModels}
#' 
#' \item{method}{Summary of the splitMethod used.}
#' @author Thomas Gerds \email{tag@@biostat.ku.dk}
#' @references Fawcett, T. (2006). An introduction to ROC analysis. Pattern
#' Recognition Letters, 27, 861-874.
#' 
#' Gerds, Cai & Schumacher (2008). The Performance of Risk Prediction Models.
#' Biometrical Journal, Vol 50, 4, 457-479.
#' 
#' Efron, Tibshirani (1997) Journal of the American Statistical Association 92,
#' 548--560 Improvement On Cross-Validation: The .632+ Bootstrap Method.
#' 
#' Wehberg, S and Schumacher, M (2004) A comparison of nonparametric error rate
#' estimation methods in classification problems. Biometrical Journal, Vol 46,
#' 35--47
#' @keywords models
##' @examples
##' 
##' ## Generate some data with binary response Y
##' ## depending on X1 and X2 and X1*X2
##' set.seed(40)
##' N <- 40
##' X1 <- rnorm(N)
##' X2 <- abs(rnorm(N,4))
##' X3 <- rbinom(N,1,.4)
##' expit <- function(x) exp(x)/(1+exp(x))
##' lp <- expit(-2 + X1 + X2 + X3 - X3*X2)
##' Y <- factor(rbinom(N,1,lp))
##' dat <- data.frame(Y=Y,X1=X1,X2=X2)
##'
##' # single markers, one by one
##' r1 <- Roc(Y~X1,data=dat)
##' plot(r1,col=1)
##' r2 <- Roc(Y~X2,data=dat)
##' lines(r2,col=2)
##'
##' # or, directly multiple in one
##' r12 <- Roc(list(Y~X1,Y~X2),data=dat)
##' plot(r12)
##' 
##' ## compare logistic regression
##' lm1 <- glm(Y~X1,data=dat,family="binomial")
##' lm2 <- glm(Y~X1+X2,data=dat,family="binomial")
##' r1=Roc(list(LR.X1=lm1,LR.X1.X2=lm2))
##' summary(r1)
##' Brier(list(lm1,lm2))
##' 
##' # machine learning
##' library(randomForest)
##' dat$Y=factor(dat$Y)
##' rf <- randomForest(Y~X2,data=dat)
##' rocCV=Roc(list(RandomForest=rf,LogisticRegression=lm2),
##'     data=dat,
##'     splitMethod="bootcv",
##'     B=3,
##'     cbRatio=1)
##' plot(rocCV)
##'
##' # compute .632+ estimate of Brier score
##' bs <- Brier(list(LR.X1=lm1,LR.X2=lm2),
##'     data=dat,
##'     splitMethod="boot632+",
##'     B=3)
##' bs
##' #'
#' @export
# {{{ UseMethod
Roc <- function(object,...){
    UseMethod("Roc",object=object)
}
# }}}
#' @S3method Roc list
Roc.list <- function(object,
                     formula,
                     data,
                     splitMethod="noSplitMethod",
                     noinf.method=c("simulate"),
                     simulate="reeval",
                     B,
                     M,
                     breaks,
                     cbRatio=1,
                     RocAverageMethod="vertical",
                     RocAverageGrid=switch(RocAverageMethod,"vertical"=seq(0,1,.01),"horizontal"=seq(1,0,-.01)),
                     model.args=NULL,
                     model.parms=NULL,
                     keepModels=FALSE,
                     keepSampleIndex=FALSE,
                     keepCrossValRes=FALSE,
                     keepNoInfSimu,
                     slaveseed,
                     cores=1,
                     na.accept=0,
                     verbose=FALSE,
                     ...){
    # }}}
    theCall=match.call()
    if (match("replan",names(theCall),nomatch=FALSE))
        stop("Argument name 'replan' has been replaced by 'splitMethod'.")
    # {{{ models
    NF <- length(object) 
    if (is.null(names(object))){
        names(object) <- sapply(object,function(o){
            cl <- class(o)[1]
            if (cl=="formula"){
                ff <- update.formula(o,"NULL~.")
                vv <- all.vars(ff)
                switch(length(vv),
                       "1"=vv,
                       "2"=paste(vv,collapse="+"),
                       paste("Formula",length(vv),"vars",sep="."))
            } else
                cl
        })
    }
    else{
        names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o){
            cl <- class(o)[1]
            if (cl=="formula"){
                ff <- update.formula(o,"NULL~.")
                vv <- all.vars(ff)
                switch(length(vv),
                       "1"=vv,
                       "2"=paste(vv,collapse="+"),
                       paste("Formula",length(vv),"vars",sep="."))
            } else
                cl
        })}
    object.names = names(object)
    names(object) <- make.names(names(object),unique=TRUE)
    # }}}
    # {{{ formula
    brier.obj <- sapply(object,function(x){
        if ((class(x)[1]=="formula" && (length(all.vars(x))<=2))
            ||
            (class(x)[1]=="numeric" && (any(x<0)|any(x>1))))
            FALSE
        else
            TRUE
    })
    if (missing(formula)){
        if (match("formula",class(object[[1]]),nomatch=FALSE))
            formula <- object[[1]]
        else
            formula <- eval(object[[1]]$call$formula)
        if (class(formula)!="formula")
            stop("Argument formula is missing.")
        else
            if (verbose)
                warning("Argument formula is missing. I use the formula from the call to the first model instead.")
    }
    # }}}
    # {{{ data
    if (missing(data)){
        data <- eval(object[[1]]$call$data)
        if (class(data)!="data.frame")
            stop("Argument data is missing.")
        else
            if (verbose)
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
    # {{{ break points for the ROC
    if (missing(breaks))
        if (splitMethod=="noSplitMethod")
            breaks <- NULL
        else
            breaks <- seq(0,1,.01)
    
    # }}}
    # {{{ SplitMethod
    SplitMethod <- MgSplitMethods(splitMethod=splitMethod,B=B,N=N,M=M,k=k)
    if (SplitMethod$internal.name=="crossval")
        stop("K-fold cross-validation not available, but you can get similar results with bootstrap subsampling: set (1) splitMethod='bootcv' and (2) M = NROW(mydata)-round(0.1*NROW(mydata))")
    B <- SplitMethod$B
    CrossvalIndex <- SplitMethod$index
    if (!keepSampleIndex) SplitMethod$index <- NULL
    k <- SplitMethod$k
    do.crossval <- !(is.null(CrossvalIndex))
    if (missing(keepCrossValRes)) keepCrossValRes <- do.crossval
    if (missing(keepNoInfSimu)) keepNoInfSimu <- FALSE
    if (missing(slaveseed)||is.null(slaveseed))
        slaveseed <- sample(1:1000000,size=B,replace=FALSE)
    # }}}
    # {{{ checking the models for compatibility with cross-validation
    if (do.crossval){
        cm <- MgCheck(object=object,model.args=model.args,model.parms=model.parms,SplitMethod=SplitMethod,verbose=verbose)
        model.args <- cm$model.args
        model.parms <- cm$model.parms
    }
    # }}}
    # {{{ computation of ROC curves in a loop over the models 
    list.out <- lapply(1:NF,function(f){
        if (verbose && NF>1) message("\n",names(object)[f],"\n")
        fit <- object[[f]]
        extract <- model.parms[[f]]
        # }}}
        # {{{ apparent ROC (use the same data for fitting and validation)
        pred <- do.call("predictStatusProb",c(list(object=fit,newdata=data),model.args[[f]]))
        if (is.null(breaks))
            breaks.f <- sort(unique(pred))
        else
            breaks.f <- breaks
        AppRoc <- Roc.default(object=pred,y=Y,breaks=breaks.f,cbRatio=cbRatio)
        if (length(unique(pred))==1) AppAuc <- 0.5
        AppAuc <- Auc.default(object=AppRoc$Sensitivity,
                              Spec=AppRoc$Specificity)
        if (brier.obj[[f]])
            AppBS <- Brier.default(object=pred,y=Y,cbRatio=cbRatio)
        else
            AppBS <- NA
        # }}}
        # {{{ No information error  
        if (SplitMethod$internal.name %in% c("boot632plus","noinf")){
            if (noinf.method=="simulate"){
                if (verbose)
                    cat("\nSimulate no information performance\n")
                compute.NoInfRocList <- lapply(1:B,function(runb){
                    if (verbose) MgTalk(runb,B)
                    data.index <- data
                    ## permute the response variable
                    responseName <- all.vars(formula)[1]
                    data.index[,responseName] <- sample(factor(Y),replace=FALSE)
                    if (simulate=="reeval"){
                        fit.index <- MgRefit(object=fit,data=data.index,step=runb,silent=na.accept>0,verbose=verbose)
                    }
                    ## evaluate the model in data with reeallocated responses
                    else
                        fit.index <- fit
                    pred.index <- do.call("predictStatusProb",
                                          c(list(object=fit.index,newdata=data.index),
                                            model.args[[f]]))
                    innerNoInfRoc <- Roc.default(object=pred.index,y=data.index[,responseName],breaks=breaks.f,cbRatio=cbRatio)
                    innerNoInfBS <- Brier.default(object=pred.index,y=data.index[,responseName],cbRatio=cbRatio)
                    list("innerNoInfRoc"=innerNoInfRoc,"innerNoInfBS"=innerNoInfBS)
                })
                if (verbose) cat("\n")
                NoInfRocList <- lapply(compute.NoInfRocList,function(x)x$innerNoInfRoc)
                NoInfRoc <- avRoc(list=NoInfRocList,grid=RocAverageGrid,method=RocAverageMethod)
                NoInfBS <- mean(sapply(compute.NoInfRocList,function(x){x$innerNoInfBS}))
                NoInfAuc <- mean(sapply(NoInfRocList,function(nil){Auc.default(object=nil$Sensitivity,nil$Specificity)}))
            }
            else{         
                NoInfRoc <- list(Sensitivity=c(breaks.f,0),Specificity=c(1-breaks.f,1))
                NoInfAuc <- 0.5
                NoInfBS <- .C("brier_noinf",bs=double(1),as.double(Y),as.double(pred),as.integer(N),NAOK=TRUE,PACKAGE="ModelGood")$bs
            }
        }
        # }}}
        # {{{ Bootcv aka BootstrapCrossValidation
        if (SplitMethod$internal.name %in% c("boot632plus","bootcv","boot632")){
            if (verbose)
                cat("\nBootstrap cross-validation performance\n")

            compute.step <- function(runb,seed){
                if (verbose) MgTalk(runb,B)
                vindex.index <- match(1:N,CrossvalIndex[,runb],nomatch=0)==0
                val.index <- data[vindex.index,,drop=FALSE]
                train.index <- data[CrossvalIndex[,runb],,drop=FALSE]
                if (!is.null(seed)) {
                    set.seed(seed)
                }
                fit.index <- MgRefit(object=fit,data=train.index,step=runb,silent=na.accept>0,verbose=verbose)
                if (!is.null(extract)) {
                    fit.parms.index <- fit.index[extract]
                    names(fit.parms.index) <- paste(extract,paste("sample",runb,sep="."),sep=":")
                }
                else fit.parms.index <- NULL
                if (is.null(fit.index)){
                    failed <- "fit"
                    innerBootcvRoc <- list(Sensitivity=NA,Specificity=NA)
                    innerBCVBS <- NA
                }
                else{
                    try2predict <- try(pred.index <- do.call("predictStatusProb",c(list(object=fit.index,newdata=val.index),model.args[[f]])),silent=na.accept>0)
                    if (inherits(try2predict,"try-error")==TRUE){
                        if (verbose) warning(paste("During bootstrapping: prediction for model ",class(fit.index)," failed in step ",runb),immediate.=TRUE)
                        failed <- "prediction"
                        innerBootcvRoc <- list(Sensitivity=NA,Specificity=NA)
                        innerBCVBS <- NA
                    }
                    else{
                        failed <- NA
                        innerBootcvRoc <- Roc.default(y=Y[vindex.index],pred.index,breaks=breaks.f,cbRatio=cbRatio)
                        innerBCVBS <- Brier.default(object=pred.index,y=Y[vindex.index],cbRatio=cbRatio)
                    }
                }
                list("innerBootcvRoc"=innerBootcvRoc,
                     "fit.parms"=fit.parms.index,
                     "failed"=failed,
                     "innerBCVBS"=innerBCVBS)
            }
            ## if (require(foreach)){
            b <- 1
            ## require(foreach)
            ## compute.BootcvRocList <- foreach::foreach (b = 1:B) %dopar% compute.step(runb=b,seed=slaveseed[[b]])
            compute.BootcvRocList <- parallel::mclapply(1:B,function(b){
                compute.step(runb=b,seed=slaveseed[[b]])
            },mc.cores=cores)
            ## }
            ## else{
            ## compute.BootcvRocList <- lapply(1:B,compute.step)
            ## }
            if (verbose) cat("\n")
            if (!is.null(extract)) fitParms <- sapply(compute.BootcvRocList,function(x)x$fit.parms)
            failed <- na.omit(sapply(compute.BootcvRocList,function(x)x$failed))
            BootcvRocList <- lapply(compute.BootcvRocList,function(x)x$innerBootcvRoc)
            BootcvRoc <- avRoc(BootcvRocList,grid=RocAverageGrid,method=RocAverageMethod)
            BCVSens <- BootcvRoc$Sensitivity
            BCVSpec <- BootcvRoc$Specificity
            BCVAuc <- mean(sapply(BootcvRocList,function(ool){Auc.default(object=ool$Sensitivity,ool$Specificity)}))
            BCVBS <- mean(sapply(compute.BootcvRocList,function(x){x$innerBCVBS}))
        }

        # }}}
        # {{{ Bootstrap .632

        if (SplitMethod$internal.name=="boot632"){
            B632Roc <- list(Sensitivity=.368 * AppRoc$Sensitivity + .632 * BCVSens,
                            Specificity=.368 * AppRoc$Specificity + .632 * BCVSpec)
            B632BS <- .368 * AppBS + .632 * BCVBS
            B632Auc <- .368 * AppAuc + .632 * BCVAuc
        }
        # }}}
        # {{{ Bootstrap .632+
        if (SplitMethod$internal.name=="boot632plus"){
            ## first we have to prepare the averaging
            if (RocAverageMethod=="vertical"){
                AppSens <- c(approx(AppRoc$Sensitivity,AppRoc$Specificity,xout=RocAverageGrid,yleft=0,yright=1,ties=median)$y,0)
                AppRoc <- list(Sensitivity=AppSens,Specificity=c(RocAverageGrid,1))
                AppAuc <- Auc.default(object=AppRoc$Sensitivity,Spec=AppRoc$Specificity)
                
                R632Plus <- MgFormule632(App=AppSens,BCV=BCVSens,NoInf=NoInfRoc$Sensitivity,SmallerBetter=FALSE)
                B632plusRoc <- list(Sensitivity=R632Plus$B632Plus,Specificity=c(RocAverageGrid,1))
            }
            else{ ## RocAverageMethod horizontal
                AppSpec <- c(approx(AppRoc$Sensitivity,AppRoc$Specificity,xout=RocAverageGrid,yleft=1,yright=0,ties=median)$y,1)
                AppRoc <- list(Sensitivity=c(RocAverageGrid,0),Specificity=AppSpec)
                R632Plus <- MgFormule632(App=AppSpec,BCV=BCVSpec,NoInf=NoInfRoc$Specificity,SmallerBetter=FALSE)
                B632plusRoc <- list(Sensitivity=c(RocAverageGrid,0),
                                    Specificity=R632Plus$B632Plus)
            }
            ## add the real .632+ AUC
            B632PlusAuc <- MgFormule632(App=AppAuc,BCV=BCVAuc,NoInf=NoInfAuc,SmallerBetter=FALSE)
            B632PlusBS <- MgFormule632(App=AppBS,BCV=BCVBS,NoInf=NoInfBS,SmallerBetter=TRUE)
        }
        # }}}
        # {{{ preparing the output
        out <- switch(SplitMethod$internal.name,
                      "noSplitMethod"=list("Roc"=AppRoc),
                      ## "plain"=list("Roc"=BootRoc,"AppRoc"=AppRoc),
                      "boot632"=list("Roc"=B632Roc,"AppRoc"=AppRoc,"BootcvRoc"=BootcvRoc),
                      "boot632plus"=list("Roc"=B632plusRoc,
                          "AppRoc"=AppRoc,
                          "BootcvRoc"=BootcvRoc,
                          "NoInfRoc"=NoInfRoc,
                          "weight"=R632Plus$weight,
                          "overfit"=R632Plus$overfit),
                      "bootcv"=list("Roc"=BootcvRoc,"AppRoc"=AppRoc),
                      "noinf"=list("AppRoc"=AppRoc,"NoInfRoc"=NoInfRoc))
        
        out$Auc <- switch(SplitMethod$internal.name,
                          "noSplitMethod"=list("Auc"=AppAuc),
                          ## "plain"=list("Auc"=BootAuc,"AppAuc"=AppAuc),
                          "boot632"= list("Auc"=B632Auc,"AppAuc"=AppAuc,"AucBCV"=BCVAuc),
                          "boot632plus"=list("Auc"=B632PlusAuc$B632Plus,
                              "AppAuc"=AppAuc,
                              "BootcvAuc"=BCVAuc,
                              "NoInfAuc"=NoInfAuc,
                              "weight"=B632PlusAuc$weight,
                              "overfit"=B632PlusAuc$overfit),
                          "bootcv"=list("Auc"=BCVAuc,"AppAuc"=AppAuc),
                          "noinf"=list("NoInfAuc"=NoInfAuc,"AppAuc"=AppAuc))
        out$Brier <- switch(SplitMethod$internal.name,
                            "noSplitMethod"=list("BS"=AppBS),
                            ## "plain"=list("BS"=BootBS,"AppBS"=AppBS),
                            "boot632"= list("BS"=B632BS,"AppBS"=AppBS,"BSBCV"=BCVBS),
                            "boot632plus"=list("BS"=B632PlusBS$B632Plus,
                                "AppBS"=AppBS,
                                "BootcvBS"=BCVBS,
                                "NoInfBS"=NoInfBS,
                                "weight"=B632PlusBS$weight,
                                "overfit"=B632PlusBS$overfit),
                            "bootcv"=list("BS"=BCVBS,"AppBS"=AppBS),
                            "noinf"=list("AppBS"=AppBS,"NoInfBS"=NoInfBS))
        out$breaks <- breaks.f
        ##     if (keepCrossValRes==TRUE && SplitMethod$internal.name!="noSplitMethod"){
        if (keepCrossValRes==TRUE && class(try(is.null(BootcvRocList),silent=TRUE))!="try-error"){
            if (SplitMethod$internal.name!="noinf")
                out <- c(out,list("BootcvRocList"=BootcvRocList))
        }
        if (keepNoInfSimu==TRUE && SplitMethod$internal.name!="noSplitMethod"){
            if (SplitMethod$internal.name!="noinf")
                out <- c(out,list("NoInfRocList"=NoInfRocList))
        }
        if (!is.null(extract)) out <- c(out,list("fitParms"=fitParms))
        if (na.accept>0) out <- c(out,list("failed"=failed))
        out
    })
    ## it might be that the first model has no extracted parameters
    ## but one of the other has
    if (length(model.parms)>0)
        names.lout <- c(names(list.out[[1]]),"fitParms")
    else
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
                  Response=Y,
                  models=outmodels,
                  model.names=object.names,
                  method=SplitMethod,
                  ## breaks=breaks,
                  cbRatio=cbRatio))
    if (verbose) cat("\n")
    # }}}
    class(out) <- "Roc"
    out
}

# {{{ Roc.default, Roc.glm,etc.
#' @method Roc default
#' @S3method Roc default
Roc.default <- function(object,
                        y,
                        breaks,
                        cbRatio=1,
                        ## pv=FALSE,
                        ## confint=FALSE,
                        ## confint.method="exact",
                        keep.tables=FALSE,
                        keep.breaks=FALSE,...){
    N <- length(object)
    if(length(y)!=N) stop("Arguments must have the same length")
    if(length(unique(y))!=2) stop("y must be binary")
    Disease <- as.integer(as.character(factor(y,labels=c("0","1"))))
    count.DiseasePos <- sum(Disease==1)
    count.DiseaseNeg <- sum(Disease==0)
    ## print(breaks)
    if (missing(breaks) || is.null(breaks))
        breaks <- sort(unique(object))
    else
        breaks <- sort(unique(breaks))
    if (length(breaks)>1 & !is.factor(object)){
        eval.times <- breaks- min(diff(breaks))/2
        count.TestPos <- N-prodlim::sindex(jump.times=object,eval.times=eval.times)
        count.TestNeg <- N-count.TestPos
        a <- count.DiseasePos-prodlim::sindex(jump.times=object[Disease==1],eval.times=eval.times)
        b <- count.DiseaseNeg-prodlim::sindex(jump.times=object[Disease==0],eval.times=eval.times)
    }
    else{
        tabx <- table(object)
        count.TestPos <- tabx
        count.TestNeg <- tabx
        tabxpos <- table(object[Disease==1])
        tabxneg <- table(object[Disease==0])
        a <- count.DiseasePos-tabxpos
        b <- count.DiseaseNeg-tabxneg
    }
    c <- count.DiseasePos-a
    d <- count.DiseaseNeg-b
    sens <- cbRatio*a/count.DiseasePos
    spec <- d/count.DiseaseNeg
    out <- list("Sensitivity"=c(sens,0),"Specificity"=c(spec,1))
    ##   if (confint==TRUE){
    ##     tmp.sens <- binconf(x=a,n=count.DiseasePos,method=confint.method)
    ##     sens <- tmp.sens[,"PointEst"]
    ##     ci.sens <- tmp.sens[,c("Lower","Upper")]
    ##     tmp.spec <- binconf(x=d,n=count.DiseaseNeg,method=confint.method)
    ##     spec <- tmp.spec[,"PointEst"]
    ##     ci.spec <- tmp.spec[,c("Lower","Upper")]
    ##     out <- list("Sensitivity"=c(sens,0),"Specificity"=c(spec,1),
    ##                 "CI.Sens"=rbind(ci.sens,c(0,0)),"CI.Spec"=rbind(ci.spec,c(1,1)))
    ##   }
    ##   if (pv==TRUE){
    ##     if (confint==TRUE){
    ##       tmp.ppv <- binconf(x=a,n=count.TestPos,method=confint.method)
    ##       ppv <- tmp.ppv[,"PointEst"]
    ##       ci.ppv <- tmp.ppv[,c("Lower","Upper")]
    ##       tmp.npv <- binconf(x=d,n=count.TestNeg,method=confint.method)
    ##       npv <- tmp.npv[,"PointEst"]
    ##       ci.npv <- tmp.npv[,c("Lower","Upper")]
    ##       out <- c(out,list("PPV"=c(ppv,1),"NPV"=c(npv,npv[length(npv)]),
    ##                         list("CI.PPV"=rbind(ci.ppv,c(1,1)),
    ##                              "CI.NPV"=rbind(ci.npv,ci.npv[NCOL(ci.npv),]))))
    ##     }
    ##     else{
    ##   ppv <- a/count.TestPos
    ##   npv <- d/count.TestNeg
    ## }
    ##   out <- c(out,list("PPV"=c(ppv,1),"NPV"=c(npv,npv[length(npv)])))
    ## }
    if (keep.breaks==TRUE)
        out <- c(out,list(breaks=breaks))
    if (keep.tables==TRUE){
        out <- c(out,list(tables=data.frame(a,b,c,d)))
    }
    ##   if (confint==TRUE)
    ##     out$confint.method <- confint.method
    ## class(out) <- "Roc"
    out
}

#' @method Roc formula
#' @S3method Roc formula
Roc.formula <- function(object,formula,data,...){
    ff <- update.formula(object,"NULL~.")
    vv <- all.vars(ff)
    obj <- list(object)
    names(obj) <- switch(length(vv),
                         "1"=vv,
                         "2"=paste(vv,collapse="+"),
                         paste("Formula",length(vv),"vars",sep="."))
    Roc.list(object=obj,formula=object,data,...)
}

#' @S3method Roc numeric
Roc.numeric <- function(object,formula,data,...){
    Roc.list(object=list(object),formula,data,...)
}

#' @S3method Roc integer
Roc.integer <- function(object,formula,data,...){
    Roc.list(object=list(object),formula,data,...)
}

#' @method Roc glm
#' @S3method Roc glm
Roc.glm <- function(object,formula,data,...){
    stopifnot(object$family$family=="binomial")
    Roc.list(object=list(object),formula,data,...)
}

#' @method Roc lrm
#' @S3method Roc lrm
Roc.lrm <- function(object,formula,data,...){
    Roc.list(object=list(object),formula,data,...)
}

#' @method Roc ElasticNet
#' @S3method Roc ElasticNet
Roc.ElasticNet <- function(object,formula,data,...){
  Roc.list(list(object),formula,data,...)
}

#' @method Roc rpart
#' @S3method Roc rpart
Roc.rpart <- function(object,formula,data,...){
    Roc.list(object=list(object),formula,data,...)
}

#' @method Roc randomForest
#' @S3method Roc randomForest
Roc.randomForest <- function(object,formula,data,...){
    Roc.list(object=list(object),formula,data,...)
}

# }}}
# {{{ average Roc curves
avRoc <- function(list,grid,method="vertical"){
    if (missing(grid))
        grid <- switch(method,"vertical"=seq(0,1,.01),"horizontal"=seq(1,0,-.01))
    if (method=="vertical"){
        meanSens <- rowMeans(do.call("cbind",lapply(list,function(Roc){
            approx(x=Roc$Specificity,y=Roc$Sensitivity,xout=grid,ties=median,yleft=0,yright=1)$y
        })))
        meanRoc <- list(Sensitivity=c(meanSens,0),Specificity=c(grid,1))
    }
    else
        if (method=="horizontal"){
            meanSpec <- rowMeans(do.call("cbind",lapply(list,function(Roc){
                approx(x=Roc$Sensitivity,y=Roc$Specificity,xout=grid,ties=median,yleft=0,yright=1)$y
            })))
            meanRoc <- list(Sensitivity=c(grid,0),Specificity=c(meanSpec,1))
        }
    meanRoc
}
# }}}


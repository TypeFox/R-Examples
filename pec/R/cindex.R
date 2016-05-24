#' Concordance index for right censored survival time data
#' 
#' In survival analysis, a pair of patients is called concordant if the risk of
#' the event predicted by a model is lower for the patient who experiences the
#' event at a later timepoint. The concordance probability (C-index) is the
#' frequency of concordant pairs among all pairs of subjects. It can be used to
#' measure and compare the discriminative power of a risk prediction models.
#' The function provides an inverse of the probability of censoring weigthed
#' estimate of the concordance probability to adjust for right censoring.
#' Cross-validation based on bootstrap resampling or bootstrap subsampling can
#' be applied to assess and compare the discriminative power of various
#' regression modelling strategies on the same set of data.
#' 
#' Pairs with identical observed times, where one is uncensored and one is
#' censored, are always considered usuable (independent of the value of
#' \code{tiedOutcomeIn}), as it can be assumed that the event occurs at a later
#' timepoint for the censored observation.
#' 
#' For uncensored response the result equals the one obtained with the
#' functions \code{rcorr.cens} and \code{rcorrcens} from the \code{Hmisc}
#' package (see examples).
#' 
#' @aliases cindex 
#' @param object A named list of prediction models, where allowed entries are
#' (1) R-objects for which a \link{predictSurvProb} method exists (see
#' details), (2) a \code{call} that evaluates to such an R-object (see
#' examples), (3) a matrix with predicted probabilities having as many rows as
#' \code{data} and as many columns as \code{times}. For cross-validation all
#' objects in this list must include their \code{call}.
#' @param formula A survival formula. The left hand side is used to finde the
#' status response variable in \code{data}. For right censored data, the right
#' hand side of the formula is used to specify conditional censoring models.
#' For example, set \code{Surv(time,status)~x1+x2} and \code{cens.model="cox"}.
#' Then the weights are based on a Cox regression model for the censoring times
#' with predictors x1 and x2.  Note that the usual coding is assumed:
#' \code{status=0} for censored times and that each variable name that appears
#' in \code{formula} must be the column name in \code{data}. If there are no
#' covariates, i.e. \code{formula=Surv(time,status)~1} the \code{cens.model} is
#' coerced to \code{"marginal"} and the Kaplan-Meier estimator for the
#' censoring times is used to calculate the weights.  If \code{formula} is
#' \code{missing}, try to extract a formula from the first element in object.
#' @param data A data frame in which to validate the prediction models and to
#' fit the censoring model.  If \code{data} is missing, try to extract a data
#' set from the first element in object.
#' @param eval.times A vector of timepoints for evaluating the discriminative
#' ability. At each timepoint the c-index is computed using only those pairs
#' where one of the event times is known to be earlier than this timepoint. If
#' \code{eval.times} is \code{missing} or \code{Inf} then the largest
#' uncensored event time is used.
#' @param pred.times A vector of timepoints for evaluating the prediction
#' models. This should either be exactly one timepoint used for all
#' \code{eval.times}, or be as long as \code{eval.times}, in which case the
#' predicted order of risk for the jth entry of \code{eval.times} is based on
#' the jth entry of \code{pred.times} corresponding
#' @param cause For competing risks, the event of interest. Defaults to the
#' first state of the response, which is obtained by evaluating the left hand
#' side of \code{formula} in \code{data}.
#' @param lyl If \code{TRUE} rank subjects accoring to predicted
#' life-years-lost (See Andersen due to this cause instead of predicted risk.
#' @param cens.model Method for estimating inverse probability of censoring
#' weigths:
#' 
#' \code{cox}: A semi-parametric Cox proportional hazard model is fitted to the
#' censoring times
#' 
#' \code{marginal}: The Kaplan-Meier estimator for the censoring times
#' 
#' \code{nonpar}: Nonparametric extension of the Kaplan-Meier for the censoring
#' times using symmetric nearest neighborhoods -- available for arbitrary many
#' strata variables on the right hand side of argument \code{formula} but at
#' most one continuous variable. See the documentation of the functions
#' \code{prodlim} and \code{neighborhood} from the prodlim package.
#' 
#' \code{aalen}: The nonparametric Aalen additive model fitted to the censoring
#' times. Requires the timereg package maintained by Thomas Scheike.
#' @param ipcw.refit If \code{TRUE} the inverse probability of censoring
#' weigths are estimated separately in each training set during
#' cross-validation.
#' @param ipcw.args List of arguments passed to function specified by argument \code{cens.model}.
#' @param ipcw.limit Value between 0 and 1 (but no equal to 0!) used to cut for
#' small weights in order to stabilize the estimate at late times were few
#' individuals are observed.
#' @param tiedPredictionsIn If \code{FALSE} pairs with identical predictions
#' are excluded, unless also the event times are identical and uncensored and
#' \code{tiedMatchIn} is set to \code{TRUE}.
#' @param tiedOutcomeIn If \code{TRUE} pairs with identical and uncensored
#' event times are excluded, unless also the predictions are identical and
#' \code{tiedMatchIn} is set to \code{TRUE}.
#' @param tiedMatchIn If \code{TRUE} then pairs with identical predictions and
#' identical and uncensored event times are counted as concordant pairs.
#' @param splitMethod Defines the internal validation design:
#' 
#' \code{none/noPlan}: Assess the models in the give \code{data}, usually
#' either in the same data where they are fitted, or in independent test data.
#' 
#' \code{BootCv}: Bootstrap cross validation. The prediction models are trained
#' on \code{B} bootstrap samples, that are either drawn with replacement of the
#' same size as the original data or without replacement from \code{data} of
#' the size \code{M}.  The models are assessed in the observations that are NOT
#' in the bootstrap sample.
#' 
#' \code{Boot632}: Linear combination of AppCindex and OutOfBagCindex using the
#' constant weight .632.
#' 
#' @param B Number of bootstrap samples. The default depends on argument
#' \code{splitMethod}.  When \code{splitMethod} in c("BootCv","Boot632") the
#' default is 100.  For \code{splitMethod="none"} \code{B} is the number of
#' bootstrap simulations e.g. to obtain bootstrap confidence limits -- default
#' is 0.
#' @param M The size of the bootstrap samples for resampling without
#' replacement. Ignored for resampling with replacement.
#' @param model.args List of extra arguments that can be passed to the
#' \code{predictSurvProb} methods. The list must have an entry for each entry
#' in \code{object}.
#' @param model.parms Experimental.  List of with exactly one entry for each
#' entry in \code{object}.  Each entry names parts of the value of the fitted
#' models that should be extracted and added to the value.
#' @param keep.index Logical. If \code{FALSE} remove the bootstrap or
#' cross-validation index from the output list which otherwise is included in
#' the method part of the output list.
#' @param keep.matrix Logical. If \code{TRUE} add all \code{B} prediction error
#' curves from bootstrapping or cross-validation to the output.
#' @param keep.models Logical. If \code{TRUE} keep the models in object. Since
#' fitted models can be large objects the default is \code{FALSE}.
#' @param keep.residuals Experimental.
#' @param keep.pvalues Experimental.
#' @param keep.weights Experimental.
#' @param multiSplitTest Experimental.
#' @param testTimes A vector of time points for testing differences between
#' models in the time-point specific Brier scores.
#' @param confInt Experimental.
#' @param confLevel Experimental.
#' @param verbose if \code{TRUE} report details of the progress, e.g. count the
#' steps in cross-validation.
#' @param savePath Place in your filesystem (directory) where training models
#' fitted during cross-validation are saved. If \code{missing} training models
#' are not saved.
#' @param slaveseed Vector of seeds, as long as \code{B}, to be given to the
#' slaves in parallel computing.
#' @param na.action Passed immediately to model.frame. Defaults to na.fail. If
#' set otherwise most prediction models will not work.
#' @param ... Not used.
#' @return Estimates of the C-index.
#' @author Thomas A Gerds \email{tag@@biostat.ku.dk}
#' @references
#' 
#' TA Gerds, MW Kattan, M Schumacher, and C Yu. Estimating a time-dependent
#' concordance index for survival prediction models with covariate dependent
#' censoring. Statistics in Medicine, Ahead of print:to appear, 2013. DOI =
#' 10.1002/sim.5681
#' 
#' Wolbers, M and Koller, MT and Witteman, JCM and Gerds, TA (2013) Concordance
#' for prognostic models with competing risks Research report 13/3. Department
#' of Biostatistics, University of Copenhagen
#' 
#' Andersen, PK (2012) A note on the decomposition of number of life years lost
#' according to causes of death Research report 12/2. Department of
#' Biostatistics, University of Copenhagen
#' @keywords survival
#' @examples
#' 
#'  # simulate data based on Weibull regression  
#' library(prodlim)
#'  set.seed(13)
#'  dat <- SimSurv(100)
#'  # fit three different Cox models and a random survival forest
#'  # note: low number of trees for the purpose of illustration 
#'  library(survival)
#'  library(randomForestSRC)
#'  cox12 <- coxph(Surv(time,status)~X1+X2,data=dat)
#'  cox1 <- coxph(Surv(time,status)~X1,data=dat)
#'  cox2 <- coxph(Surv(time,status)~X2,data=dat)
#'  rsf1 <- rfsrc(Surv(time,status)~X1+X2,data=dat,ntree=15,forest=TRUE)
#'  #
#'  # compute the apparent estimate of the C-index at different time points
#'  #
#' ApparrentCindex  <- cindex(list("Cox X1"=cox1,
#' 		       "Cox X2"=cox2,
#' 		       "Cox X1+X2"=cox12,
#'                        "RSF"=rsf1),
#' 		  formula=Surv(time,status)~X1+X2,
#' 		  data=dat,
#' 		  eval.times=seq(5,500,50))
#'   print(ApparrentCindex)
#'   plot(ApparrentCindex)
#'  #
#'  # compute the bootstrap-crossvalidation estimate of
#'  # the C-index at different time points
#'  #
#' set.seed(142)
#' bcvCindex  <- cindex(list("Cox X1"=cox1,
#' 		       "Cox X2"=cox2,
#' 		       "Cox X1+X2"=cox12,
#'                        "RSF"=rsf1),
#' 		  formula=Surv(time,status)~X1+X2,
#' 		  data=dat,
#'                   splitMethod="bootcv",
#'                   B=5,
#' 		  eval.times=seq(5,500,50))
#'   print(bcvCindex)
#'   plot(bcvCindex)
#'  # for uncensored data the results are the same
#'  # as those obtained with the function rcorr.cens from Hmisc
#' library(Hmisc)
#' set.seed(16)
#' dat <- SimSurv(30,cens=FALSE)
#' fit12 <- coxph(Surv(time,status)~X1+X2,data=dat)
#' fit1 <- coxph(Surv(time,status)~X1,data=dat)
#' fit2 <- coxph(Surv(time,status)~X2,data=dat)
#' Cpec <- cindex(list("Cox X1+X2"=fit12,"Cox X1"=fit1,"Cox X2"=fit2),
#' 	       formula=Surv(time,status)~1,
#' 	       data=dat,
#' 	       eval.times=Inf)
#' p1 <- predictSurvProb(fit1,newdata=dat,times=100)
#' p2 <- predictSurvProb(fit2,newdata=dat,times=100)
#' p12 <- predictSurvProb(fit12,newdata=dat,times=100)
#' harrelC1 <- rcorr.cens(p1,with(dat,Surv(time,status)))
#' harrelC2 <- rcorr.cens(p2,with(dat,Surv(time,status)))
#' harrelC12 <- rcorr.cens(p12,with(dat,Surv(time,status)))
#' harrelC1[["C Index"]]==Cpec$AppCindex[["Cox.X1"]]
#' harrelC2[["C Index"]]==Cpec$AppCindex[["Cox.X2"]]
#' harrelC12[["C Index"]]==Cpec$AppCindex[["Cox.X1.X2"]]
#'  #
#'  # competing risks 
#'  #
#' library(riskRegression)
#' library(prodlim)
#' set.seed(30)
#' dcr.learn <- SimCompRisk(30)
#' dcr.val <- SimCompRisk(30)
#' cindex(CSC(Hist(time,event)~X1+X2,data=dcr.learn),data=dcr.val)
#'
#' @export 
# {{{ header cindex.list
cindex <- function(object,
                   formula,
                   data,
                   eval.times,
                   pred.times,
                   cause,
                   lyl=FALSE,
                   cens.model="marginal",
                   ipcw.refit=FALSE,
                   ipcw.args=NULL,
                   ipcw.limit,
                   tiedPredictionsIn=TRUE,
                   tiedOutcomeIn=TRUE,
                   tiedMatchIn=TRUE,
                   splitMethod="noPlan",
                   B,
                   M,
                   model.args=NULL,
                   model.parms=NULL,
                   keep.models=FALSE,
                   keep.residuals=FALSE,
                   keep.pvalues=FALSE,
                   keep.weights=FALSE,
                   keep.index=FALSE,
                   keep.matrix=FALSE,
                   multiSplitTest=FALSE,
                   testTimes,
                   confInt=FALSE,
                   confLevel=0.95,
                   verbose=TRUE,
                   savePath=NULL,
                   slaveseed=NULL,
                   na.action=na.fail,
                   ...){

  # }}}
# {{{ checking integrity some arguments
  theCall=match.call()
  if (match("replan",names(theCall),nomatch=FALSE))
    stop("Argument name 'replan' has been replaced by 'splitMethod'.")
  
  if (keep.residuals && missing(testTimes))
      stop("To keep.residuals please specify testTimes.")
  if (missing(splitMethod) && multiSplitTest==TRUE){
      stop("Need data splitting to compute van de Wiel's test")
  }
  if (missing(M) && multiSplitTest) M <- NA
  stopifnot(as.numeric(tiedPredictionsIn) %in% c(0,1))
  stopifnot(as.numeric(tiedOutcomeIn) %in% c(0,1))
  stopifnot(as.numeric(tiedMatchIn) %in% c(0,1))
  # }}}
  # {{{ check and convert object
  if (class(object)[1]!="list") {
      object <- list(object)
  }
  # }}}
  # {{{  formula
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
      stop("Invalid specification of formula. Perhaps forgotten right hand side?\nNote that any subsetting, ie data$var or data[,\"var\"], is invalid for this function.")
  }
  else{
      if (!(formula.names[2] %in% c("Surv","Hist")))
          survp <- FALSE
      else
          survp <- TRUE
  }
  # }}}
  # {{{  data
  if (missing(data)){
      if (match("call",names(object[[1]]),nomatch=0)==0||is.null(object[[1]]$call$data)){
          stop("Data missing and cannot borrow data from the first object :(")
      }
      data <- eval(object[[1]]$call$data)
      if (match("data.frame",class(data),nomatch=0)==0)
          stop("Argument data is missing.")
      else
          if (verbose)
              warning("Argument data is missing. I use the data from the call to the first model instead.")
  }

# }}}
# {{{  censoring model
  
  cens.model <- match.arg(cens.model,c("cox","marginal","nonpar","aalen","none"))
  
  # }}}
  # {{{ response
  histformula <- formula
  if (histformula[[2]][[1]]==as.name("Surv")){
      histformula[[2]][[1]] <- as.name("Hist")
  }
  m <- model.frame(histformula,data,na.action=na.action)
  response <- model.response(m)
  if (match("Surv",class(response),nomatch=0)!=0){
      attr(response,"model") <- "survival"
      attr(response,"cens.type") <- "rightCensored"
      model.type <- "survival"
  }
  censType <- attr(response,"cens.type")
  model.type <- attr(response,"model")
  if (model.type=="competing.risks"){
      if (lyl==TRUE)
          predictHandlerFun <- "predictLifeYearsLost"
      else
          predictHandlerFun <- "predictEventProb"
      if (missing(cause))
          cause <- attr(response,"state")[1]
  }
  else{
      if (survp==FALSE && NCOL(response)!=1) stop("Response must be one-dimensional.")
      if (survp==TRUE && NCOL(response)!=2) stop("Survival response must have two columns: time and status.")
      predictHandlerFun <- "predictSurvProb"
  }
  if (model.type=="competing.risks")
      if (verbose==TRUE) message("Cindex for competing risks")
  # }}}
  # {{{ prediction models
  NF <- length(object) 
  if (is.null(names(object))){
      names(object) <- sapply(object,function(o)class(o)[1])
      names(object) <- make.names(names(object),unique=TRUE)
  }
  else{ # fix missing names
      if (any(names(object)=="")){
          names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
          names(object) <- make.names(names(object),unique=TRUE)
      }else{
           # leave names as they were given
       }
  }
  # }}}
  # {{{ sort the data 

  if (survp){
    neworder <- order(response[,"time"],-response[,"status"])
    if (model.type=="competing.risks"){
      event <- prodlim::getEvent(response,mode="character")
      event <- event[neworder]
    }
    response <- response[neworder,,drop=FALSE]
    Y <- response[,"time"]
    if (censType=="uncensored"){
      status <- rep(1,length(Y))
      cens.model <- "none"
    }
    else{
      status <- response[,"status"]
    }
  }
  else{
    cens.model <- "none"
    neworder <- order(response)
    Y <- response[neworder]
    status <- rep(1,length(Y))
  }
  ## for competing risks find the cause of interest.
  if (model.type=="competing.risks"){
    availableCauses <- unique(event)
    if (!match(cause, availableCauses,nomatch=FALSE))
      stop("Cause ",cause," is not among the available causes: ",paste(availableCauses,collapse=", "))
    event <- event==cause
  }
  else{
    event <- NULL
  }
  data <- data[neworder,]
  unique.Y <- unique(Y)
  N <- length(Y)
  NU <- length(unique.Y)

  # }}}
  # {{{  splitMethod
  splitMethod <- resolvesplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
  if (splitMethod$internal.name %in% c("Boot632plus")) stop(".632+ method not implemented for c-index.")
  B <- splitMethod$B
  ResampleIndex <- splitMethod$index
  k <- splitMethod$k
  do.resample <- !(is.null(ResampleIndex))
  if (keep.matrix==TRUE & !do.resample){
    warning("Argument keep.matrix set to FALSE, since no resampling/crossvalidation is requested.")
    keep.matrix <- FALSE
  }
  # }}}
  # {{{  define the prediction time(s) and the evaluation time(s)
  maxtime <- unique.Y[NU]
  if (missing(eval.times) || is.infinite(eval.times)){
    ## eval.times <- max(Y) ## maybe less efficient 
    eval.times <- max(Y[status==1])
  }
  else{
    tooLate <- sum(eval.times>maxtime)
    if (tooLate>0){
      if (verbose)
          warning(tooLate," eval.times beyond the maximal evaluation time: ",ifelse(maxtime>1,round(maxtime,1),round(maxtime,3)))
      ## eval.times <- c(eval.times[eval.times<maxtime],maxtime)
    }
  }
  if (missing(pred.times))
    ## pred.times <- median(unique.Y)
    pred.times <- eval.times
  ## FIXME: if the model changes the risk order over time, then we need to care about pred.times
  NT <-  length(eval.times)
  tindex <- match(Y,unique.Y)
  # }}}
  # {{{  IPCW (all equal to 1 without censoring)
  if((cens.model %in% c("aalen","cox","nonpar"))){
    if (all(as.numeric(status)==1) || sum(status)==N){
      message("No censored observations: cens.model coerced to \"none\".")
      cens.model <- "none"
    }
    if (length(attr(terms(formula),"factors"))==0){
      if (verbose==TRUE)
        message("No covariates  specified: cens.model coerced to \"marginal\".\n")
      cens.model <- "marginal"}
  }
  if (model.type=="competing.risks"){
    iFormula <- as.formula(paste("Surv(itime,istatus)","~",as.character(formula)[[3]]))
    iData <- data
    iData$itime <- response[,"time"]
    iData$istatus <- response[,"status"]
    weight.i <- ipcw(formula=iFormula,data=iData,method=cens.model,times=NULL,subjectTimes=Y,subjectTimesLag=1,what="IPCW.subjectTimes")$IPCW.subjectTimes
    weight.j <- ipcw(formula=iFormula,data=iData,method=cens.model,times=unique.Y,subjectTimes=NULL,subjectTimesLag=0,what="IPCW.times")$IPCW.times
    ipcw.call <- NULL
  }
  else{
    #  weights for T_i<=T_j
    #  FIXED: the correct weights are G(T_i|X_j) and G(T_i-|X_i)
    weight.i <- ipcw(formula=formula,data=data,method=cens.model,times=NULL,subjectTimes=Y,subjectTimesLag=1,what="IPCW.subjectTimes")$IPCW.subjectTimes
    weight.j <- ipcw(formula=formula,data=data,method=cens.model,times=unique.Y,subjectTimes=NULL,subjectTimesLag=0,what="IPCW.times")$IPCW.times
  }
  if (ipcw.refit==TRUE)
    stop("pec: internal refitting of censoring distribution not (not yet) supported for competing risks")
  ## if (ipcw.refit==TRUE && splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632"))
  ## ipcw.call <- list(weight.i=list(formula=formula,data=data,method=cens.model,times=NULL,subjectTimes=Y,subjectTimesLag=1,what="IPCW.subjectTimes"),
  ## weight.j=list(formula=formula,data=data,method=cens.model,times=unique.Y,subjectTimes=NULL,subjectTimesLag=0,what="IPCW.times"))
  ## else
  ipcw.call <- NULL
  ## print(weight.j)
  # truncate the weights
  if (!missing(ipcw.limit) && ipcw.limit!=0){
    pfit <- prodlim::prodlim(update(formula,".~1"),data=data)
    limit.i <- 1/(1+c(0,cumsum(pfit$n.lost))[1+prodlim::sindex(jump.times=pfit$time,eval.times=Y-min(diff(Y))/2)])
    limit.times <- 1/(1+c(0,cumsum(pfit$n.lost))[1+prodlim::sindex(jump.times=pfit$time,eval.times=unique.Y)])
    if (ipcw.limit<1 && ipcw.limit>0){
      limit.i <- pmax(ipcw.limit,limit.i)
      limit.times <- pmax(ipcw.limit,limit.times)
    }
    weight.i <- pmax(weight.i,limit.i)
    if (is.null(dim(weight.j))){
      weight.j <- pmax(weight.j,limit.times)
    }
    else{
      weight.j <- t(apply(weight.j,1,function(wj){
        pmax(limit.times,wj)
      }))
    }
  }
  weights <- list(weight.i=weight.i,weight.j=weight.j)
  # }}}
  # {{{  checking the models for compatibility with resampling
  if (do.resample){
    cm <- checkModels(object=object,model.args=model.args,model.parms=model.parms,splitMethod=splitMethod$internal.name)
    model.args <- cm$model.args
    model.parms <- cm$model.parms
  }
  # }}}
  # {{{ -------------------Apparent or test sample cindex----------------------
  
  AppCindexList <- lapply(1:NF,function(f){
    fit <- object[[f]]
    extraArgs <- model.args[[f]]
    if (model.type=="competing.risks"){
      pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=pred.times,cause=cause),extraArgs))
      if (class(fit)[[1]]%in% c("matrix","numeric")) pred <- pred[neworder,]
      if (length(pred.times)==1 && length(pred.times)<length(eval.times))
        pred <- rep(pred,length(eval.times))
      ## stop("Not yet defined: cindex for competing risks")
      AppCindexResult <- .C("ccr",cindex=double(NT),concA=double(NT),pairsA=double(NT),concB=double(NT),pairsB=double(NT),as.integer(tindex),as.double(Y),as.integer(status),as.integer(event),as.double(eval.times),as.double(weight.i),as.double(weight.j),as.double(pred),as.integer(N),as.integer(NT),as.integer(tiedPredictionsIn),as.integer(tiedOutcomeIn),as.integer(tiedMatchIn),as.integer(!is.null(dim(weight.j))),NAOK=TRUE,package="pec")
      AppCindex <- AppCindexResult$cindex
      AppPairsA <- AppCindexResult$pairsA
      AppConcordantA <- AppCindexResult$concA
      AppPairsB <- AppCindexResult$pairsB
      AppConcordantB <- AppCindexResult$concB
      list(AppCindex=AppCindex,
           AppPairs=list(A=AppPairsA,B=AppPairsB),
           AppConcordant=list(A=AppConcordantA,B=AppConcordantB))
    }
    else{
      pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=pred.times),extraArgs))
      if (class(fit)[[1]]%in% c("matrix","numeric")) pred <- pred[neworder,]
      if (length(pred.times)==1 && length(pred.times)<length(eval.times))
        pred <- rep(pred,length(eval.times))
      AppCindexResult <- .C("cindex",cindex=double(NT),conc=double(NT),pairs=double(NT),as.integer(tindex),as.double(Y),as.integer(status),as.double(eval.times),as.double(weight.i),as.double(weight.j),as.double(pred),as.integer(N),as.integer(NT),as.integer(tiedPredictionsIn),as.integer(tiedOutcomeIn),as.integer(tiedMatchIn),as.integer(!is.null(dim(weight.j))),NAOK=TRUE,package="pec")
      AppCindex <- AppCindexResult$cindex
      AppPairs <- AppCindexResult$pairs
      AppConcordant <- AppCindexResult$conc
      list(AppCindex=AppCindex,
           AppPairs=AppPairs,
           AppConcordant=AppConcordant)
    }

  })
  AppCindex <- lapply(AppCindexList,function(x){
    x$AppCindex
  })
  if (predictHandlerFun=="predictSurvProb"){
    AppPairs <- lapply(AppCindexList,function(x){
      2*x$AppPairs
    })
    AppConcordant <- lapply(AppCindexList,function(x){
      2*x$AppConcordant
    })
  }else{
        AppPairs <- lapply(AppCindexList,function(x){
      x$AppPairs
    })
    AppConcordant <- lapply(AppCindexList,function(x){
      x$AppConcordant
    })
  }
  names(AppCindex) <- names(object)
  names(AppPairs) <- names(object)
  names(AppConcordant) <- names(object)
  
  # }}}
  # {{{--------------k-fold and leave-one-out CrossValidation-----------------------
  
  if (splitMethod$internal.name %in% c("crossval","loocv")){
      kCV <- CindexKFoldCrossValidation(object=object,
                                        data=data,
                                        Y=Y,
                                        status=status,
                                        event=event,
                                        tindex=tindex,
                                        eval.times=eval.times,
                                        pred.times=pred.times,
                                        cause=cause,
                                        weights=weights,
                                        ipcw.refit=ipcw.refit,
                                        ipcw.call=ipcw.call,
                                        tiedPredictionsIn=tiedPredictionsIn,
                                        tiedOutcomeIn=tiedOutcomeIn,
                                        tiedMatchIn=tiedMatchIn,
                                        splitMethod=splitMethod,
                                        multiSplitTest=multiSplitTest,
                                        keepResiduals=keep.residuals,
                                        testTimes=testTimes,
                                        confInt=confInt,
                                        confLevel=confLevel,
                                        getFromModel=model.parms,
                                        giveToModel=model.args,
                                        predictHandlerFun=predictHandlerFun,
                                        keepMatrix=keep.matrix,
                                        verbose=verbose,
                                        savePath=savePath,
                                        slaveseed=slaveseed)
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
    BootCv <- CindexBootstrapCrossValidation(object=object,data=data,Y=Y,status=status,event=event,eval.times=eval.times,pred.times=pred.times,cause=cause,weights=weights,ipcw.refit=ipcw.refit,ipcw.call=ipcw.call,splitMethod=splitMethod,multiSplitTest=multiSplitTest,testTimes=testTimes,confInt=confInt,confLevel=confLevel,getFromModel=model.parms,giveToModel=model.args,predictHandlerFun=predictHandlerFun,tiedPredictionsIn=tiedPredictionsIn,tiedOutcomeIn=tiedOutcomeIn,tiedMatchIn=tiedMatchIn,keepMatrix=keep.matrix,keepResiduals=keep.residuals,verbose=verbose,savePath=savePath,slaveseed=slaveseed)
    BootstrapCrossValCindex <- BootCv$BootstrapCrossValCindex
    Residuals <- BootCv$Residuals
    names(BootstrapCrossValCindex) <- names(object)
    if (multiSplitTest==TRUE){
      comparisons <- allComparisons(names(object))
      multiSplitTest <- list(B=B,M=M,N=N,testTimes=testTimes)
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
        ## if (keep.pvalues==TRUE){
        ## out$allPairwisePvaluesIBS <- allPairwisePvaluesIBS}
        out
      })
      names(multiSplitTest$Comparisons) <- names(comparisons)
      class(multiSplitTest) <- "multiSplitTest"
    }
    if (keep.matrix==TRUE){
      BootstrapCrossValCindexMat <- BootCv$BootstrapCrossValCindexMat
      names(BootstrapCrossValCindex) <- names(object)
    }
  }

  # }}}
  # {{{ Bootstrap .632
  if (splitMethod$internal.name=="Boot632"){
      B632Cindex <- lapply(1:NF,function(f){
          .368 * AppCindex[[f]] + .632 * BootstrapCrossValCindex[[f]]
      })
      names(B632Cindex) <- names(object)
  }
  # }}}    
  # {{{ prepare output
  out <- switch(splitMethod$internal.name,
                "noPlan"=list("AppCindex"=AppCindex,
                    "Pairs"=AppPairs,
                    "Concordant"=AppConcordant),
                "Boot632"=list("AppCindex"=AppCindex,
                    ## "Pairs"=AppPairs,
                    ## "Concordant"=AppConcordant,
                    "BootCvCindex"= BootstrapCrossValCindex,
                    "Boot632Cindex"=B632Cindex),
                "BootCv"=list("AppCindex"=AppCindex,
                    ## "Pairs"=AppPairs,
                    ## "Concordant"=AppConcordant,
                    "BootCvCindex"=BootstrapCrossValCindex
                    ## "BootCvConcordant"=BootCvConcordant,
                    ## "BootCvPairs"=BootCvPairs
                    ))
  observed.maxtime <- sapply(out,function(x){
      lapply(x,function(y){
          eval.times[length(y)-sum(is.na(y))]
      })
  })
  minmaxtime <- min(unlist(observed.maxtime))
  
  
  if (multiSplitTest==TRUE){
      out <- c(out,list(multiSplitTest=multiSplitTest))
  }
  if (keep.residuals==TRUE){
      out <- c(out,list(Residuals=Residuals))
  }
  if (keep.matrix==TRUE && splitMethod$internal.name!="noPlan"){
      if (splitMethod$internal.name %in% c("crossval","loocv")){
          if (B>1)
              out <- c(out,list("BootstrapCrossValCindexMat"=BootstrapCrossValCindexMat))
      }
      else{
          if (splitMethod$internal.name!="noinf")
              out <- c(out,list("BootstrapCrossValCindexMat"=BootstrapCrossValCindexMat))
      }
  }
  if (!is.null(model.parms)) out <- c(out,list("ModelParameters"=BootCv$ModelParameters))
  if (!keep.index) splitMethod$index <- NULL
  if(keep.models==TRUE){
      outmodels <- object
  } else{
        outmodels <- names(object)
        names(outmodels) <- names(object)
    }
  out <- c(out,list(call=theCall,
                    time=eval.times,
                    pred.time=pred.times,
                    response=model.response(m),
                    models=outmodels,
                    splitMethod=splitMethod,
                    weights=weights,
                    cens.model=cens.model,
                    minmaxtime=minmaxtime,
                    maxtime=maxtime))
  ## if (verbose==TRUE && do.resample==TRUE) cat("\n")
  # }}}
  class(out) <- "Cindex"
  out
}
    
  

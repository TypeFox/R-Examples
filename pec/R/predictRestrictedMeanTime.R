# methods for survival regression
# --------------------------------------------------------------------
#' Predicting restricted mean time 
#' 
#' Function to extract predicted mean times from various modeling
#' approaches. 
#' 
#' The function predictRestrictedMeanTime is a generic function, meaning that it
#' invokes a different function dependent on the 'class' of the
#' first argument. 
#' 
#' See also \code{\link{predictSurvProb}}.
#' 
#' @aliases predictRestrictedMeanTime predictRestrictedMeanTime.aalen
#' predictRestrictedMeanTime.riskRegression predictRestrictedMeanTime.cox.aalen
#' predictRestrictedMeanTime.coxph predictRestrictedMeanTime.cph predictRestrictedMeanTime.default
#' predictRestrictedMeanTime.rfsrc predictRestrictedMeanTime.matrix predictRestrictedMeanTime.pecCtree
#' predictRestrictedMeanTime.prodlim predictRestrictedMeanTime.psm
#' predictRestrictedMeanTime.selectCox predictRestrictedMeanTime.survfit 
#' predictRestrictedMeanTime.pecRpart
#' @usage
#' \method{predictRestrictedMeanTime}{aalen}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{riskRegression}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{cox.aalen}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{cph}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{coxph}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{matrix}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{selectCox}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{prodlim}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{psm}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{survfit}(object,newdata,times,...)
#' \method{predictRestrictedMeanTime}{pecRpart}(object,newdata,times,...)
#' #' \method{predictRestrictedMeanTime}{pecCtree}(object,newdata,times,...)
#' @param object A fitted model from which to extract predicted survival
#' probabilities
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted survival probabilities.
#' @param times A vector of times in the range of the response variable, e.g.
#' times when the response is a survival object, at which to return the
#' survival probabilities.
#' @param \dots Additional arguments that are passed on to the current method.
#' @return A matrix with as many rows as \code{NROW(newdata)} and as many
#' columns as \code{length(times)}. Each entry should be a probability and in
#' rows the values should be decreasing.
#' @note In order to assess the predictive performance of a new survival model
#' a specific \code{predictRestrictedMeanTime} S3 method has to be written. For examples,
#' see the bodies of the existing methods.
#' 
#' The performance of the assessment procedure, in particular for resampling
#' where the model is repeatedly evaluated, will be improved by supressing in
#' the call to the model all the computations that are not needed for
#' probability prediction. For example, \code{se.fit=FALSE} can be set in the
#' call to \code{cph}.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{predict}},\code{\link{survfit}}
#' @references Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23. URL
#' http://www.jstatsoft.org/v50/i11/.
#' @keywords survival
##' @examples
##' 
##' # generate some survival data
##' library(prodlim)
##' set.seed(100)
##' d <- SimSurv(100)
##' # then fit a Cox model
##' library(rms)
##' coxmodel <- cph(Surv(time,status)~X1+X2,data=d,surv=TRUE)
##' 
##' # predicted survival probabilities can be extracted
##' # at selected time-points:
##' ttt <- quantile(d$time)
##' # for selected predictor values:
##' ndat <- data.frame(X1=c(0.25,0.25,-0.05,0.05),X2=c(0,1,0,1))
##' # as follows
##' predictRestrictedMeanTime(coxmodel,newdata=ndat,times=ttt)
##' 
##' # stratified cox model
##' sfit <- coxph(Surv(time,status)~strata(X1)+X2,data=d,y=TRUE)
##' predictRestrictedMeanTime(sfit,newdata=d[1:3,],times=c(1,3,5,10))
##' 
##' ## simulate some learning and some validation data
##' learndat <- SimSurv(100)
##' valdat <- SimSurv(100)
##' ## use the learning data to fit a Cox model
##' library(survival)
##' fitCox <- coxph(Surv(time,status)~X1+X2,data=learndat)
##' ## suppose we want to predict the survival probabilities for all patients
##' ## in the validation data at the following time points:
##' ## 0, 12, 24, 36, 48, 60
##' psurv <- predictRestrictedMeanTime(fitCox,newdata=valdat,times=seq(0,60,12))
##' ## This is a matrix with survival probabilities
##' ## one column for each of the 5 time points
##' ## one row for each validation set individual
##' 
##' # the same can be done e.g. for a randomSurvivalForest model
##' library(randomForestSRC)
##' rsfmodel <- rfsrc(Surv(time,status)~X1+X2,data=d)
##' predictRestrictedMeanTime(rsfmodel,newdata=ndat,times=ttt)
 
#' @export 
predictRestrictedMeanTime <- function(object,newdata,times,...){
    UseMethod("predictRestrictedMeanTime",object)
}

##' @export
predictRestrictedMeanTime.default <- function(object,newdata,times,...){
  stop("No method for evaluating predicted probabilities from objects in class: ",class(object),call.=FALSE)
}


##' @export
predictRestrictedMeanTime.numeric <- function(object,newdata,times,...){
  if (NROW(object) != NROW(newdata) || NCOL(object) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(object)," x ",NCOL(object),"\n\n",sep=""))
  object
}

##' @export
predictRestrictedMeanTime.matrix <- function(object,newdata,times,...){
    if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(object)," x ",NCOL(object),"\n\n",sep=""))
        ## stop(paste("Prediction matrix has wrong dimensionss: ",NROW(object)," rows and ",NCOL(object)," columns.\n But requested are predicted probabilities for ",NROW(newdata), " subjects (rows) in newdata and ",NCOL(newdata)," time points (columns)",sep=""))
    }
    object
}

##' @export
predictRestrictedMeanTime.aalen <- function(object,newdata,times,...){
    ## require(timereg)
    time.coef <- data.frame(object$cum)
    ntime <- nrow(time.coef)
    objecttime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    covanames <- names(time.coef)[-(1:2)]
    notfound <- match(covanames,names(newdata),nomatch=0)==0
    if (any(notfound))
        stop("\nThe following predictor variables:\n\n",
             paste(covanames[notfound],collapse=","),
             "\n\nwere not found in newdata, which only provides the following variables:\n\n",
             paste(names(newdata),collapse=","),
             "\n\n")
    time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
    nobs <- nrow(newdata)
    hazard <- .C("survest_cox_aalen",
                 timehazard=double(ntime*nobs),
                 as.double(unlist(time.coef[,-1])),
                 as.double(unlist(time.vars)),
                 as.integer(ntimevars+1),
                 as.integer(nobs),
                 as.integer(ntime),PACKAGE="pec")$timehazard
    hazard <- matrix(hazard,ncol=ntime,nrow=nobs,dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
    surv <- pmin(exp(-hazard),1)
    if (missing(times)) times <- sort(unique(objecttime))
    p <- surv[,prodlim::sindex(jump.times=objecttime,eval.times=times)]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictRestrictedMeanTime.cox.aalen <- function(object,newdata,times,...){
    #  require(timereg)
    ##  The time-constant effects first
    const <- c(object$gamma)
    names(const) <- substr(dimnames(object$gamma)[[1]],6,nchar(dimnames(object$gamma)[[1]])-1)
    constant.part <- t(newdata[,names(const)])*const
    constant.part <- exp(colSums(constant.part))
    ##  Then extract the time-varying effects
    time.coef <- data.frame(object$cum)
    ntime <- nrow(time.coef)
    objecttime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
    nobs <- nrow(newdata)
    time.part <- .C("survest_cox_aalen",timehazard=double(ntime*nobs),as.double(unlist(time.coef[,-1])),as.double(unlist(time.vars)),as.integer(ntimevars+1),as.integer(nobs),as.integer(ntime),PACKAGE="pec")$timehazard
    time.part <- matrix(time.part,ncol=ntime,nrow=nobs)
    ## dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
    surv <- pmin(exp(-time.part*constant.part),1)
    if (missing(times)) times <- sort(unique(objecttime))
    p <- surv[,prodlim::sindex(objecttime,times)]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictRestrictedMeanTime.pecRpart <- function(object,newdata,times,...){
    newdata$rpartFactor <- factor(predict(object$rpart,newdata=newdata),
                                  levels=object$levels)
    p <- predictRestrictedMeanTime(object$survfit,newdata=newdata,times=times)
    p
}
    
##' @export
predictRestrictedMeanTime.coxph <- function(object,newdata,times,...){
    ## baselineHazard.coxph(object,times)
    ## require(survival)
    ## new feature of the survival package requires that the
    ## original data are included
    ## survival.survfit.coxph <- getFromNamespace("survfit.coxph",ns="survival")
    ## survival.summary.survfit <- getFromNamespace("summary.survfit",ns="survival")
    ## b <- function(x){browser()}
    ## b()
    eTimes <- unique(sort(object$y))
    pos <- prodlim::sindex(jump.times=eTimes,eval.times=times)
    surv <- predictSurvProb(object,newdata=newdata,times=eTimes)
    rmt <- matrix(unlist(lapply(1:length(pos), function(j) {
                                            pos.j <- 1:(pos[j]+1)
                                            p <- cbind(1,surv)[,pos.j,drop=FALSE]
                                            time.diff <- diff(c(0, eTimes)[pos.j])
                                            apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
                                        })), ncol = length(pos))
    if ((miss.time <- (length(times) - NCOL(rmt)))>0)
        rmt <- cbind(rmt,matrix(rep(NA,miss.time*NROW(rmt)),nrow=NROW(rmt)))
    if (NROW(rmt) != NROW(newdata) || NCOL(rmt) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",
                   NROW(newdata),
                   " x ",
                   length(times),
                   "\nProvided prediction matrix: ",
                   NROW(rmt),
                   " x ",
                   NCOL(rmt),
                   "\n\n",
                   sep=""))
    rmt
}

##' @export
predictRestrictedMeanTime.coxph.penal <- function(object,newdata,times,...){
  ## require(survival)
  frailhistory <- object$history$'frailty(cluster)'$history
  frailVar <- frailhistory[NROW(frailhistory),1]
  ## survfit.object <- survival.survfit.coxph(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
  linearPred <- predict(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
  basehaz <- basehaz(object)
  bhTimes <- basehaz[,2]
  bhValues <- basehaz[,1]
  survPred <- do.call("rbind",lapply(1:NROW(newdata),function(i){
    (1+frailVar*bhValues*exp(linearPred[i]))^{-1/frailVar}
  }))
  where <- prodlim::sindex(jump.times=bhTimes,eval.times=times)
  p <- cbind(1,survPred)[,where+1]
  if ((miss.time <- (length(times) - NCOL(p)))>0)
    p <- cbind(p,matrix(rep(NA,miss.time*NROW(p)),nrow=NROW(p)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}



##' @export
predictRestrictedMeanTime.cph <- function(object,newdata,times,...){
    if (!match("surv",names(object),nomatch=0)) stop("Argument missing: set surv=TRUE in the call to cph!")
    p <- rms::survest(object,times=times,newdata=newdata,se.fit=FALSE,what="survival")$surv
    if (is.null(dim(p))) p <- matrix(p,nrow=NROW(newdata))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictRestrictedMeanTime.selectCox <- function(object,newdata,times,...){
    predictRestrictedMeanTime(object[[1]],newdata=newdata,times=times,...)
}

##' @export
predictRestrictedMeanTime.prodlim <- function(object,newdata,times,...){
    ## require(prodlim)
    p <- predict(object=object,
                 type="surv",
                 newdata=newdata,
                 times=times,
                 mode="matrix",
                 level.chaos=1)
    if (NROW(newdata)==1 && class(p)=="list"){
        p <- unlist(p)
    }
    if (is.null(dim(p)) && NROW(newdata)>=1){
        ## if the model has no covariates
        ## then all cases get the same prediction
        ## in this exceptional case we return a vector
        ## p[is.na(p)] <- 0
        p <- as.vector(p)
        if (length(p)!=length(times))
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    else{
        if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
        ## stop("Prediction failed")
    }
    rownames(p) <- NULL
    p
}

predict.survfit <- function(object,newdata,times,bytimes=TRUE,fill="last",...){
    if (length(class(object))!=1 || class(object)!="survfit" || object$typ !="right")
        stop("Predictions only available \nfor class 'survfit', possibly stratified Kaplan-Meier fits.\n For class 'cph' Cox models see survest.cph.")
    if (missing(newdata))
        npat <- 1
    else
        if (is.data.frame(newdata))
            npat <- nrow(newdata)
        else stop("If argument `newdata' is supplied it must be a dataframe." )
    ntimes <- length(times)
    sfit <- summary(object,times=times)
    if (is.na(fill))
        Fill <- function(x,len){x[1:len]}
    else if (fill=="last")
        Fill <- function(x,len){
            y <- x[1:len]
            y[is.na(y)] <- x[length(x)]
            y}
    else stop("Argument fill must be the string 'last' or NA.")
    if (is.null(object$strata)){
        pp <- Fill(sfit$surv,ntimes)
        p <- matrix(rep(pp,npat),
                    ncol=ifelse(bytimes,ntimes,npat),
                    nrow=ifelse(bytimes,npat,ntimes),
                    byrow=bytimes)
    }
    else{
        covars <- attr(terms(eval.parent(object$call$formula)),"term.labels")
        if (!all(match(covars,names(newdata),nomatch=FALSE)))
            stop("Not all strata defining variables occur in newdata.")
        ## FIXME there are different ways to build strata levels
        ## how can we test which one was used???
        stratdat <- newdata[,covars,drop=FALSE]
        names(stratdat) <- covars
        NewStratVerb <- survival::strata(stratdat)
        NewStrat <- interaction(stratdat,sep=" ")
        levs <- levels(sfit$strata)
        #    print(levs)
        #    print(levels(NewStrat))
        #    print(levels(NewStratVerb))
        if (!all(choose <- match(NewStratVerb,levs,nomatch=F))
            &&
            !all(choose <- match(NewStrat,levs,nomatch=F)))
            stop("Not all strata levels in newdata occur in fit.")
        survlist <- split(sfit$surv,sfit$strata)
        pp <- lapply(survlist[choose],Fill,ntimes)
        p <- matrix(unlist(pp,use.names=FALSE),
                    ncol=ifelse(bytimes,ntimes,npat),
                    nrow=ifelse(bytimes,npat,ntimes),
                    byrow=bytimes)
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictRestrictedMeanTime.survfit <- function(object,newdata,times,...){
    p <- predict.survfit(object,newdata=newdata,times=times,bytimes=TRUE,fill="last")
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}


## library randomSurvivalForest
## predictRestrictedMeanTime.rsf <- function(object,newdata,times,...){
## p <- predict.rsf(object,newdata=newdata,times=times,bytimes=TRUE,fill="last")
## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
## stop("Prediction failed")
## p
## }


##' @export
predictRestrictedMeanTime.psm <- function(object,newdata,times,...){
    if (length(times)==1){
        p <- rms::survest(object,times=c(0,times),newdata=newdata,what="survival",conf.int=FALSE)[,2]
    }else{
        p <- rms::survest(object,times=times,newdata=newdata,what="survival",conf.int=FALSE)
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}


##' @export
predictRestrictedMeanTime.riskRegression <- function(object,newdata,times,...){
    if (missing(times))stop("Argument times is missing")
    temp <- predict(object,newdata=newdata)
    pos <- prodlim::sindex(jump.times=temp$time,eval.times=times)
    p <- cbind(1,1-temp$cuminc)[,pos+1,drop=FALSE]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictRestrictedMeanTime.rfsrc <- function(object, newdata, times, ...){
    ptemp <- predict(object,newdata=newdata,importance="none",...)$survival
    pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
    p <- cbind(1,ptemp)[,pos+1,drop=FALSE]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}


# methods for uncensored regression
# --------------------------------------------------------------------

predictProb <- function(object,newdata,times,...){
  UseMethod("predictProb",object)
}

##' @export
predictProb.glm <- function(object,newdata,times,...){
  ## no censoring -- only normal family with mu=0 and sd=sd(y)
  N <- NROW(newdata)
  NT <- length(times)
  if (!(unclass(family(object))$family=="gaussian"))
    stop("Currently only gaussian family implemented for glm.")
  betax <- predict(object,newdata=newdata,se.fit=FALSE)
  ##   print(betax[1:10])
  pred.matrix <- matrix(rep(times,N),byrow=TRUE,ncol=NT,nrow=N)
  p <- 1-pnorm(pred.matrix - betax,mean=0,sd=sqrt(var(object$y)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}


##' @export
predictProb.ols <- function(object,newdata,times,...){
    ## no censoring -- only normal family with mu=0 and sd=sd(y)
    N <- NROW(newdata)
    NT <- length(times)
    if (!(unclass(family(object))$family=="gaussian"))
        stop("Currently only gaussian family implemented.")
    betax <- predict(object,newdata=newdata,type="lp",se.fit=FALSE)
    ##   print(betax[1:10])
    pred.matrix <- matrix(rep(times,N),byrow=TRUE,ncol=NT,nrow=N)
    p <- 1-pnorm(pred.matrix - betax,mean=0,sd=sqrt(var(object$y)))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictProb.randomForest <- function(object,newdata,times,...){
  ## no censoring -- only normal family with mu=0 and sd=sd(y)
  N <- NROW(newdata)
  NT <- length(times)
  predMean <- predict(object,newdata=newdata,se.fit=FALSE)
  pred.matrix <- matrix(rep(times,N),byrow=TRUE,ncol=NT,nrow=N)
  p <- 1-pnorm(pred.matrix - predMean,mean=0,sd=sqrt(var(object$y)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}



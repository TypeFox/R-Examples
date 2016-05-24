#' Predicting survival probabilities
#' 
#' Function to extract survival probability predictions from various modeling
#' approaches. The most prominent one is the Cox regression model which can be
#' fitted for example with `coxph' and with `cph'.
#' 
#' The function predictSurvProb is a generic function that means it invokes
#' specifically designed functions depending on the 'class' of the first
#' argument.
#' 
#' The function \code{pec} requires survival probabilities for each row in
#' newdata at requested times. These probabilities are extracted from a fitted
#' model of class \code{CLASS} with the function \code{predictSurvProb.CLASS}.
#' 
#' Currently there are \code{predictSurvProb} methods for objects of class cph
#' (library rms), coxph (library survival), aalen (library timereg), cox.aalen
#' (library timereg), 
#' rpart (library rpart), product.limit (library prodlim),
#' survfit (library survival), psm (library rms)
#' 
#' @aliases predictSurvProb predictSurvProb.aalen
#' predictSurvProb.riskRegression predictSurvProb.cox.aalen
#' predictSurvProb.coxph predictSurvProb.cph predictSurvProb.default
#' predictSurvProb.rfsrc predictSurvProb.matrix predictSurvProb.pecCtree
#' predictSurvProb.pecCforest predictSurvProb.prodlim predictSurvProb.psm
#' predictSurvProb.selectCox predictSurvProb.survfit 
#' predictSurvProb.pecRpart
#' @usage
#' \method{predictSurvProb}{aalen}(object,newdata,times,...)
#' \method{predictSurvProb}{riskRegression}(object,newdata,times,...)
#' \method{predictSurvProb}{cox.aalen}(object,newdata,times,...)
#' \method{predictSurvProb}{cph}(object,newdata,times,...)
#' \method{predictSurvProb}{coxph}(object,newdata,times,...)
#' \method{predictSurvProb}{matrix}(object,newdata,times,...)
#' \method{predictSurvProb}{selectCox}(object,newdata,times,...)
#' \method{predictSurvProb}{pecCforest}(object,newdata,times,...)
#' \method{predictSurvProb}{prodlim}(object,newdata,times,...)
#' \method{predictSurvProb}{psm}(object,newdata,times,...)
#' \method{predictSurvProb}{survfit}(object,newdata,times,...)
#' \method{predictSurvProb}{pecRpart}(object,newdata,times,...)
#' #' \method{predictSurvProb}{pecCtree}(object,newdata,times,...)
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
#' a specific \code{predictSurvProb} S3 method has to be written. For examples,
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
##' # Extract predicted survival probabilities 
##' # at selected time-points:
##' ttt <- quantile(d$time)
##' # for selected predictor values:
##' ndat <- data.frame(X1=c(0.25,0.25,-0.05,0.05),X2=c(0,1,0,1))
##' # as follows
##' predictSurvProb(coxmodel,newdata=ndat,times=ttt)
##' 
##' # stratified cox model
##' sfit <- coxph(Surv(time,status)~strata(X1)+X2,data=d,y=TRUE)
##' predictSurvProb(sfit,newdata=d[1:3,],times=c(1,3,5,10))
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
##' psurv <- predictSurvProb(fitCox,newdata=valdat,times=seq(0,60,12))
##' ## This is a matrix with survival probabilities
##' ## one column for each of the 5 time points
##' ## one row for each validation set individual
##' 
##' # Do the same for a randomSurvivalForest model
##' library(randomForestSRC)
##' rsfmodel <- rfsrc(Surv(time,status)~X1+X2,data=d)
##' predictSurvProb(rsfmodel,newdata=ndat,times=ttt)
##' 
##' ## Cox with ridge option
##' f1 <- coxph(Surv(time,status)~X1+X2,data=learndat)
##' f2 <- coxph(Surv(time,status)~ridge(X1)+ridge(X2),data=learndat)
##' plot(predictSurvProb(f1,newdata=valdat,times=10),
##'      pec:::predictSurvProb.coxph(f2,newdata=valdat,times=10),
##'      xlim=c(0,1),
##'      ylim=c(0,1),
##'      xlab="Unpenalized predicted survival chance at 10",
##'      ylab="Ridge predicted survival chance at 10")
##'
##' 
#' @export 
predictSurvProb <- function(object,newdata,times,...){
    UseMethod("predictSurvProb",object)
}

##' @export 
predictSurvProb.default <- function(object,newdata,times,...){
  stop("No method for evaluating predicted probabilities from objects in class: ",class(object),call.=FALSE)
}


##' @export 
predictSurvProb.numeric <- function(object,newdata,times,...){
    if (NROW(object) != NROW(newdata))
        ## || NCOL(object) != length(times))
      stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(object)," x ",NCOL(object),"\n\n",sep=""))
  object
}

##' @export 
predictSurvProb.matrix <- function(object,newdata,times,...){
    if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(object)," x ",NCOL(object),"\n\n",sep=""))
        ## stop(paste("Prediction matrix has wrong dimensionss: ",NROW(object)," rows and ",NCOL(object)," columns.\n But requested are predicted probabilities for ",NROW(newdata), " subjects (rows) in newdata and ",NCOL(newdata)," time points (columns)",sep=""))
    }
    object
}

##' @export 
predictSurvProb.aalen <- function(object,newdata,times,...){
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
predictSurvProb.cox.aalen <- function(object,newdata,times,...){
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

#' Combines the rpart result with a stratified Kaplan-Meier (prodlim) to predict survival
#'
#' 
#' @title Predict survival based on rpart tree object
#' @param formula passed to rpart
#' @param data passed to rpart
#' @param ... passed to rpart
#' @return list with three elements: ctree and call
#' @examples
#' library(prodlim)
#' library(rpart)
#' library(survival)
#' set.seed(50)
#' d <- SimSurv(50)
#' nd <- data.frame(X1=c(0,1,0),X2=c(-1,0,1))
#' f <- pecRpart(Surv(time,status)~X1+X2,data=d)
#' predictSurvProb(f,newdata=nd,times=c(3,8))
#' @export 
pecRpart <- function(formula,data,...){
    robj <- rpart::rpart(formula=formula,data=data,...)
    nclass <- length(unique(robj$where))
    data$rpartFactor <- factor(predict(robj,newdata=data,...))
    form <- update(formula,paste(".~","rpartFactor",sep=""))
    survfit <- prodlim::prodlim(form,data=data)
    out <- list(rpart=robj,survfit=survfit,levels=levels(data$rpartFactor))
    class(out) <- "pecRpart"
    out
}

##' @export
predictSurvProb.pecRpart <- function(object,newdata,times,...){
    newdata$rpartFactor <- factor(predict(object$rpart,newdata=newdata),
                                  levels=object$levels)
    p <- predictSurvProb(object$survfit,newdata=newdata,times=times)
    p
}
    
##' @export 
predictSurvProb.coxph <- function(object,newdata,times,...){
    ## baselineHazard.coxph(object,times)
    ## require(survival)
    ## new feature of the survival package requires that the
    ## original data are included
    ## survival.survfit.coxph <- getFromNamespace("survfit.coxph",ns="survival")
    ## survival.summary.survfit <- getFromNamespace("summary.survfit",ns="survival")
     ## b <- function(x){browser()}
     ## b()
    survfit.object <- survival::survfit(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
    if (is.null(attr(object$terms,"specials")$strata)){
        ## case without strata 
        inflated.pred <- summary(survfit.object,times=times)$surv
        p <- t(inflated.pred)        
    } else{
        ## case with strata 
        inflated.pred <- summary(survfit.object,times=times)
        plist <- split(inflated.pred$surv,inflated.pred$strata)
        p <- do.call("rbind",lapply(plist,function(x){
            beyond <- length(times)-length(x)
            c(x,rep(NA,beyond))
        }))
        ## p <- matrix(inflated.pred,ncol=length(times),byrow=TRUE)
    }
    if ((miss.time <- (length(times) - NCOL(p)))>0)
        p <- cbind(p,matrix(rep(NA,miss.time*NROW(p)),nrow=NROW(p)))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictSurvProb.coxph.penal <- function(object,newdata,times,...){
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
predictSurvProb.cph <- function(object,newdata,times,...){
    if (!match("surv",names(object),nomatch=0)) stop("Argument missing: set surv=TRUE in the call to cph!")
    p <- rms::survest(object,times=times,newdata=newdata,se.fit=FALSE,what="survival")$surv
    if (is.null(dim(p))) p <- matrix(p,nrow=NROW(newdata))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictSurvProb.selectCox <- function(object,newdata,times,...){
    predictSurvProb(object[[1]],newdata=newdata,times=times,...)
}

##' @export
predictSurvProb.prodlim <- function(object,newdata,times,...){
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
predictSurvProb.survfit <- function(object,newdata,times,...){
    p <- predict.survfit(object,newdata=newdata,times=times,bytimes=TRUE,fill="last")
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}


## library randomSurvivalForest
## predictSurvProb.rsf <- function(object,newdata,times,...){
## p <- predict.rsf(object,newdata=newdata,times=times,bytimes=TRUE,fill="last")
## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
## stop("Prediction failed")
## p
## }


##' @export
predictSurvProb.psm <- function(object,newdata,times,...){
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
predictSurvProb.riskRegression <- function(object,newdata,times,...){
    if (missing(times))stop("Argument times is missing")
    temp <- predict(object,newdata=newdata)
    pos <- prodlim::sindex(jump.times=temp$time,eval.times=times)
    p <- cbind(1,1-temp$cuminc)[,pos+1,drop=FALSE]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictSurvProb.rfsrc <- function(object, newdata, times, ...){
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


## update.cox <- function(object,tstar,data){
## object$call$data <- data[data$time>tstar,]
## update <- eval(object$call)
## class(update) <- "dynamicCox"
## update
## }
##' @export
## predictProb.dynamicCox <- function(object,newdata,cutpoints,learn.data,...){
## p <- matrix(1,nrow=NROW(newdata),ncol=length(cutpoints))
## p
## }

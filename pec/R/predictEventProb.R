# methods for competing risk regression
# --------------------------------------------------------------------
#' Predicting event probabilities (cumulative incidences) in competing risk
#' models.
#' 
#' Function to extract event probability predictions from various modeling
#' approaches. The most prominent one is the combination of cause-specific Cox
#' regression models which can be fitted with the function \code{cumincCox}
#' from the package \code{compRisk}.
#' 
#' The function predictEventProb is a generic function that means it invokes
#' specifically designed functions depending on the 'class' of the first
#' argument.
#' 
#' See \code{\link{predictSurvProb}}.
#' 
#' @aliases predictEventProb predictEventProb.CauseSpecificCox
#' predictEventProb.riskRegression predictEventProb.FGR
#' predictEventProb.prodlim predictEventProb.rfsrc
#' @param object A fitted model from which to extract predicted event
#' probabilities
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted event probabilities.
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
#' @param \dots Additional arguments that are passed on to the current method.
#' @return A matrix with as many rows as \code{NROW(newdata)} and as many
#' columns as \code{length(times)}. Each entry should be a probability and in
#' rows the values should be increasing.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso See \code{\link{predictSurvProb}}.
#' @keywords survival
#' @examples
#' 
#' library(pec)
#' library(CoxBoost)
#' library(survival)
#' library(riskRegression)
#' library(prodlim)
#' train <- SimCompRisk(100)
#' test <- SimCompRisk(10)
#' cox.fit  <- CSC(Hist(time,cause)~X1+X2,data=train)
#' predictEventProb(cox.fit,newdata=test,times=seq(1:10),cause=1)
#' ## cb.fit <- coxboost(Hist(time,cause)~X1+X2,cause=1,data=train,stepno=10)
#' ## predictEventProb(cb.fit,newdata=test,times=seq(1:10),cause=1)
#'
#' ## with strata
#' cox.fit2  <- CSC(list(Hist(time,cause)~strata(X1)+X2,Hist(time,cause)~X1+X2),data=train)
#' predictEventProb(cox.fit2,newdata=test,times=seq(1:10),cause=1)
#' 
#' @export 
predictEventProb <- function(object,newdata,times,cause,...){
  UseMethod("predictEventProb",object)
}

##' @export
predictEventProb.matrix <- function(object,newdata,times,...){
  if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
    stop(paste("Prediction matrix has wrong dimensions: ",
               NROW(object),
               " rows and ",
               NCOL(object),
               " columns.\n But requested are predicted probabilities for ",
               NROW(newdata),
               " subjects (rows) in newdata and ",
               length(times),
               " time points (columns)",
               sep=""))
  }
  object
}

##' @export 
predictEventProb.prodlim <- function(object,newdata,times,cause,...){
    ## require(prodlim)
    p <- predict(object=object,cause=cause,type="cuminc",newdata=newdata,times=times,mode="matrix",level.chaos=1)
    ## if the model has no covariates
    ## then all cases get the same prediction
    ## in this exceptional case we proceed a vector
    if (NROW(p)==1 && NROW(newdata)>=1)
        p <- as.vector(p)
    ## p[is.na(p)] <- 0
    if (is.null(dim(p)))
        {if (length(p)!=length(times))
             stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
     }
    else{
        if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    }
    p
}

##' @export 
predictEventProb.FGR <- function(object,newdata,times,cause,...){
    ## require(cmprsk)
    ## predict.crr <- cmprsk:::predict.crr
    p <- predict(object=object,newdata=newdata,times=times)
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
                stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export
predictEventProb.riskRegression <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata,times=times)
  pos <- prodlim::sindex(jump.times=temp$time,eval.times=times)
  p <- cbind(0,temp$risk)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}

##' @export 
predictEventProb.ARR <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata,times=times)
  pos <- prodlim::sindex(jump.times=temp$time,eval.times=times)
  p <- cbind(0,temp$P1)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}


##' @export 
predictEventProb.CauseSpecificCox <- function (object, newdata, times, cause, ...) {
    survtype <- object$survtype
    N <- NROW(newdata)
    ## suppressMessages(browser())
    NC <- length(object$model)
    ## browser(skipCalls=4)
    if (length(cause)>1)
        stop(paste0("Can only predict one cause. Provided are: ",
                    paste(cause,collapse=", "),
                    sep=""))
    eTimes <- object$eventTimes
    if (missing(cause))
        cause <- object$theCause
    causes <- object$causes
    stopifnot(match(as.character(cause),causes,nomatch=0)!=0)
    if (survtype=="survival"){
        if (object$theCause!=cause)
            stop("Object can be used to predict cause ",object$theCause," but not ",cause,".\nNote: the cause can be specified in CSC(...,cause=).")
    }
    # predict cumulative cause specific hazards
    trycumhaz1 <- try(cumHaz1 <- -log(predictSurvProb(object$models[[paste("Cause",cause)]],times=eTimes,newdata=newdata)),silent=FALSE)
    ## trycumhaz1[is.infinite(trycumhaz1)] <- NA
    if (inherits(trycumhaz1,"try-error")==TRUE)
        stop("Prediction of cause-specific Cox model failed")
    if (length(eTimes)==1)
        Haz1 <- cumHaz1
    else
        Haz1 <- t(apply(cbind(0,cumHaz1),1,diff))
    ## it may happen that surv = 0 which implies cumhaz = inf 
    ## Haz1[is.infinite(Haz1) | is.na(Haz1)] <- 0
    if (survtype=="hazard"){
        ## browser()
        cumHazOther <- lapply(causes[-match(cause,causes)],function(c){
                                        trycumhaz <- try(cumHaz.c <- -log(predictSurvProb(object$models[[paste("Cause",c)]],times=eTimes,newdata=newdata)),silent=FALSE)
                                        if (inherits(trycumhaz,"try-error")==TRUE)
                                            stop("Prediction of cause-specific Cox model failed")
                                        cumHaz.c
                                    })
        lagsurv <- exp(-cumHaz1 - Reduce("+",cumHazOther))
        cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
    }
    else{
        tdiff <- min(diff(eTimes))/2
        trylagsurv <- try(lagsurv <- predictSurvProb(object$models[["OverallSurvival"]],times=eTimes-tdiff,newdata=newdata),silent=FALSE)
        if (inherits(trylagsurv,"try-error")==TRUE)
            stop("Prediction of overall curvival Cox model failed")
        cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
    }
    pos <- prodlim::sindex(jump.times=eTimes, eval.times=times)
    p <- cbind(0,cuminc1)[,pos+1,drop=FALSE]
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}

##' @export 
predictEventProb.rfsrc <- function(object, newdata, times, cause, ...){
  if (missing(cause)) stop("missing cause")
  if (!is.numeric(cause)) stop("cause is not numeric")
  cif <- predict(object,newdata=newdata,importance="none",...)$cif[,,cause,drop=TRUE]
  pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
  p <- cbind(0,cif)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
            stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  p
}


## predictUpdateProb.CSC <- function (object, newdata,times,horizon, cause, ...) {
  ## survtype <- object$survtype
  ## N <- NROW(newdata)
  ## NC <- length(object$model)
  ## eTimes <- object$eventTimes
  ## if (missing(cause))
    ## cause <- object$theCause
  ## causes <- object$causes
  ## stopifnot(match(as.character(cause),causes,nomatch=0)!=0)
  ## # predict cumulative cause specific hazards
  ## cumHaz1 <- -log(predictSurvProb(object$models[[paste("Cause",cause)]],times=eTimes,newdata=newdata))
  ## Haz1 <- t(apply(cbind(0,cumHaz1),1,diff))
  ## if (survtype=="hazard"){
    ## cumHazOther <- lapply(causes[-match(cause,causes)],function(c){
      ## -log(predictSurvProb(object$models[[paste("Cause",c)]],times=eTimes,newdata=newdata))
    ## })
    ## lagsurv <- exp(-cumHaz1- do.call("+",cumHazOther))
    ## cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
  ## }
  ## else{
    ## tdiff <- min(diff(eTimes))/2
    ## lagsurv <- predictSurvProb(object$models[["OverallSurvival"]],times=eTimes-tdiff,newdata=newdata)
    ## cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
  ## }
  ## pos <- prodlim::sindex(jump.times=eTimes, eval.times=times)
  ## cbind(0,cuminc1)[,pos+1,drop=FALSE]
## }


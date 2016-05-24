# methods for competing risk regression
# --------------------------------------------------------------------
#' Predicting life years lost (cumulative cumulative incidences) in competing
#' risk models.
#' 
#' Function to extract predicted life years lost from various modeling
#' approaches. The most prominent one is the combination of cause-specific Cox
#' regression models which can be fitted with the function \code{cumincCox}
#' from the package \code{compRisk}.
#' 
#' The function predictLifeYearsLost is a generic function that means it
#' invokes specifically designed functions depending on the 'class' of the
#' first argument.
#' 
#' See \code{\link{predictSurvProb}}.
#' 
#' @aliases predictLifeYearsLost predictLifeYearsLost.CauseSpecificCox
#' predictLifeYearsLost.riskRegression predictLifeYearsLost.FGR
#' predictLifeYearsLost.prodlim predictLifeYearsLost.rfsrc
#' @param object A fitted model from which to extract predicted event
#' probabilities
#' @param newdata A data frame containing predictor variable combinations for
#' which to compute predicted event probabilities.
#' @param times A vector of times in the range of the response variable, for
#' which the cumulative incidences event probabilities are computed.
#' @param cause Identifies the cause of interest among the competing events.
#' @param \dots Additional arguments that are passed on to the current method.
#' @return A matrix with as many rows as \code{NROW(newdata)} and as many
#' columns as \code{length(times)}. Each entry should be a positive value and
#' in rows the values should be increasing.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{predictSurvProb}}, \code{\link{predictEventProb}}.
#' @keywords survival
#' @examples
#' 
#' library(pec)
#' library(riskRegression)
#' library(survival)
#' library(prodlim)
#' train <- SimCompRisk(100)
#' test <- SimCompRisk(10)
#' fit <- CSC(Hist(time,cause)~X1+X2,data=train,cause=1)
#' predictLifeYearsLost(fit,newdata=test,times=seq(1:10),cv=FALSE,cause=1)
#' 
#' @export
predictLifeYearsLost <- function(object,newdata,times,cause,...){
  UseMethod("predictLifeYearsLost",object)
}

##' @export
predictLifeYearsLost.matrix <- function(object,newdata,times,...){
  if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
    stop(paste("Life-years-lost matrix has wrong dimensions: ",
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
predictLifeYearsLost.prodlim <- function(object,newdata,times,cause,...){
  ## require(prodlim)
  time.interest <- object$time
  cif <- predict(object=object,cause=cause,type="cuminc",newdata=newdata,times=time.interest,mode="matrix",level.chaos=1)
  ## if the model has no covariates
  ## then all cases get the same cif
  ## in this exceptional case we proceed a vector
  if (NROW(cif)==1 && NROW(newdata)>1)
    cif <- as.vector(cif)
  pos <- prodlim::sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
            stop(paste("\nLYL matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(lyl)," x ",NCOL(lyl),"\n\n",sep=""))
  lyl
}

##' @export
predictLifeYearsLost.FGR <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  time.interest <- sort(unique(object$crrFit$uftime))
  cif <- predict(object,newdata=newdata,times=time.interest)
  pos <- prodlim::sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
            stop(paste("\nLYL matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(lyl)," x ",NCOL(lyl),"\n\n",sep=""))
  lyl
}

##' @export
predictLifeYearsLost.riskRegression <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  time.interest <- object$time
  cif <- predict(object,newdata=newdata,times=time.interest)
  pos <- prodlim::sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
            stop(paste("\nLYL matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(lyl)," x ",NCOL(lyl),"\n\n",sep=""))
  lyl
}

##' @export
predictLifeYearsLost.ARR <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  time.interest <- object$time
  cif <- predict(object,newdata=newdata,times=time.interest)
  pos <- prodlim::sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
            stop(paste("\nLYL matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(lyl)," x ",NCOL(lyl),"\n\n",sep=""))
  lyl
}


##' @export
predictLifeYearsLost.CauseSpecificCox <- function (object, newdata, times, cause, ...) {
  survtype <- object$survtype
  N <- NROW(newdata)
  NC <- length(object$model)
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
  cumHaz1 <- -log(predictSurvProb(object$models[[paste("Cause",cause)]],times=eTimes,newdata=newdata))
  if (length(eTimes)==1)
    Haz1 <- cumHaz1
  else
    Haz1 <- t(apply(cbind(0,cumHaz1),1,diff))
  if (survtype=="hazard"){
    cumHazOther <- lapply(causes[-match(cause,causes)],function(c){
      -log(predictSurvProb(object$models[[paste("Cause",c)]],times=eTimes,newdata=newdata))
    })
    lagsurv <- exp(-cumHaz1 - Reduce("+",cumHazOther))
    cif <- t(apply(lagsurv*Haz1,1,cumsum))
  }
  else{
    tdiff <- min(diff(eTimes))/2
    lagsurv <- predictSurvProb(object$models[["OverallSurvival"]],times=eTimes-tdiff,newdata=newdata)
    cif <- t(apply(lagsurv*Haz1,1,cumsum))
  }
  pos <- prodlim::sindex(jump.times=eTimes,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, eTimes)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
            stop(paste("\nLYL matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(lyl)," x ",NCOL(lyl),"\n\n",sep=""))
  lyl
}


##' @export
predictLifeYearsLost.rfsrc <- function(object, newdata, times, cause, ...){
    if (missing(cause)) stop("missing cause")
    cif <- predict(object,newdata=newdata,importance="none",...)$cif[,,cause,drop=TRUE]
    pos <- prodlim::sindex(jump.times=object$time.interest,eval.times=times)
    lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
        pos.j <- 1:(pos[j]+1)
        p <- cbind(0,cif)[,pos.j,drop=FALSE]
        time.diff <- diff(c(0, object$time.interest)[pos.j])
        apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
    })), ncol = length(pos))
    if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
        stop(paste("\nLYL matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(lyl)," x ",NCOL(lyl),"\n\n",sep=""))
    lyl
}


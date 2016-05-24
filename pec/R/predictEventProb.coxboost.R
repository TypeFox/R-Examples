#' Formula interface for function \code{CoxBoost} of package \code{CoxBoost}.
#' 
#' Formula interface for function \code{CoxBoost} of package \code{CoxBoost}.
#' 
#' See \code{CoxBoost}.
#' @aliases coxboost
#' @param formula An event-history formula for competing risks of the
#' form \code{Hist(time,status)~sex+age} where \code{status} defines
#' competing events and right censored data. The code for right
#' censored can be controlled with argument \code{cens.code}, see man
#' page the function \code{\link{Hist}}.
#' @param data A data.frame in which the variables of formula are
#' defined.
#' @param cv If \code{TRUE} perform cross-validation to optimize the
#' parameter \code{stepno}. This calls the function \code{cv.CoxBoost}
#' whose arguments are prefix controlled, that is \code{cv.K=7} sets
#' the argument \code{K} of \code{cv.CoxBoost} to \code{7}.  If
#' \code{FALSE} use \code{stepno}.
#' @param cause The cause of interest in competing risk models.
#' @param penalty See \code{CoxBoost}.
#' @param ... Arguments passed to either \code{CoxBoost} via
#' \code{CoxBoost.arg} or to \code{cv.CoxBoost} via
#' \code{cv.CoxBoost.arg}.
#' @return See \code{CoxBoost}.
#' @author Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @seealso See \code{CoxBoost}.
#' @references See \code{CoxBoost}.
#' @keywords survival
#' @export
coxboost <- function(formula,data,cv=TRUE,cause=1,penalty,...){
  call <- match.call(expand.dots=TRUE)
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[2]=="Hist")) stop("The left hand side of formula look like this: Hist(time,event).")
  actual.terms <- terms(formula,data=data)
  formula <- eval(call$formula)
  response <- model.response(model.frame(formula,data))
  Time <- as.numeric(response[,"time"])
  if (attr(response,"model")=="competing.risks"){
    ## adapt the event variable
    Event <- rep(2,NROW(response))
    thisCause <- as.numeric(response[,"event"]==cause)
    Event[thisCause==1] <- 1
    Status <- as.numeric(response[,"status"])
    Event[Status==0] <- 0
  }
  else{
    ## survival
    Event <- as.numeric(response[,"status"])
  }
  X <- model.matrix(actual.terms,data=data)[,-c(1),drop=FALSE]## remove intercept
  if (NCOL(X)<=1) stop("CoxBoost needs at least two covariates.")
  if (missing(penalty)) penalty <- sum(Event==1)*(9)
  cv.defaults=list(maxstepno=200,K=10,penalty=penalty)
  CoxBoost.defaults=list(stepno=100,penalty=penalty)
  args <- prodlim::SmartControl(call= list(...),
                       keys=c("cv","CoxBoost"),
                       ignore=c("formula","data","cv","cause"),
                       forced=list("cv"=list(time=Time,status=Event,x=X),"CoxBoost"=list(time=Time,status=Event,x=X)),
                       defaults=list("cv"=cv.defaults,"CoxBoost"=CoxBoost.defaults),
                       ignore.case=FALSE,
                       replaceDefaults=FALSE,
                       verbose=TRUE)
  if (cv==TRUE){
    cv.step <- do.call("cv.CoxBoost",args$cv)
    args$CoxBoost$stepno <- cv.step$optimal.step
  }
  cb <- do.call("CoxBoost",args$CoxBoost)
  out <- list(coxboost=cb,
              stepno=args$CoxBoost$stepno,
              call=call,
              formula=formula,
              response=response)
  class(out) <- "coxboost"
  out
}

##' @export
predictSurvProb.coxboost <- function(object,newdata,times,...) {
  newcova <- model.matrix(terms(object$formula,data=newdata),
                          data=model.frame(object$formula,data=newdata,na.action=na.fail))[,-c(1)]
  newcova <- newcova[,object$coxboost$xnames]
  p <- predict(object$coxboost,newcova,type="risk",times=times)
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

##' @export
predictEventProb.coxboost <- function(object,newdata,times,cause,...){
  if (missing(cause)) stop("missing cause")
  if (attr(object$response,"model")!="competing.risks") stop("Not a competing risk object")
  newcova <- model.matrix(terms(object$formula,data=newdata),
                          data=model.frame(object$formula,data=newdata,na.action=na.fail))[,-c(1)]
  newcova <- newcova[,object$coxboost$xnames]
  p <- predict(object$coxboost,newdata=newcova,type="CIF",times=times)
  if (is.null(dim(p))) {
    if (length(p)!=length(times))
      stop("Prediction failed (wrong number of times)")
  }
  else{
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
  }
  p
}

##' @export
predictLifeYearsLost.coxboost <- function(object,newdata,times,cause,...){
    if (missing(cause)) stop("missing cause")
    ## if (cause!=1) stop("CoxBoost can only predict cause 1")
    if (attr(object$response,"model")!="competing.risks") stop("Not a competing risk object")
    newcova <- model.matrix(terms(object$formula,data=newdata),
                            data=model.frame(object$formula,data=newdata))[,-c(1)]
    time.interest <- sort(unique(object$coxboost$time))
    cif <- predict(object$coxboost,newdata=newcova,type="CIF",times=time.interest)
    pos <- prodlim::sindex(jump.times=time.interest,eval.times=times)
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


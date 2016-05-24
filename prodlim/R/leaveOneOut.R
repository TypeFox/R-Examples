#' Compute leave-one-out estimates
#'
#' This function is the work-horse for \code{jackknife}
#' @title Compute jackknife pseudo values.
#' @aliases leaveOneOut leaveOneOut.survival leaveOneOut.competing.risks 
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{jackknife}}
#' 
#' @param object Object of class \code{"prodlim"}.
#' @param times time points at which to compute leave-one-out
#' event/survival probabilities.
#' @param cause For competing risks the cause of interest. 
#' @param lag For survival models only. If \code{TRUE} lag the result, i.e. compute
#' S(t-) instead of S(t).
#' @param ... not used
#' @export
leaveOneOut <- function(object,times,cause,lag=FALSE,...){
    if (object$model=="survival")
        leaveOneOut.survival(object=object,times=times,lag=lag,...)
    else if (object$model=="competing.risks")
        leaveOneOut.competing.risks(object=object,times=times,cause=cause,...)
    else stop("No method for jackknifing this object.")
}

#' @export 
leaveOneOut.survival <- function(object,times,lag=FALSE,...){
  stopifnot(object$covariate.type==1)
  mr <- object$model.response
  time <- object$time
  Y <- object$n.risk
  D <- object$n.event
  Y <- Y[D>0]
  time <- time[D>0]
  D <- D[D>0]
  NU <- length(time)
  obstimes <- mr[,"time"]
  status <- mr[,"status"]
  N <- length(obstimes)
  ##
  S <- predict(object,times=time,newdata=mr)
  ## idea: find the at-risk set for pseudo-value k by
  ##       substracting 1 in the period where subj k is
  ##       at risk. need the position of obstime.k in time ...
  ## pos <- match(obstimes,time)
  ## if (useC==TRUE){
  loo <- .C("loo_surv",
            Y = as.double(Y),
            D=as.double(D),
            time=as.double(time),
            obsT=as.double(obstimes),
            status=as.double(status),
            S=double(NU*N),
            N=as.integer(N),
            NT=as.integer(NU),
            PACKAGE="prodlim")$S
  out <- matrix(loo,nrow=N,ncol=NU,byrow=FALSE)
  ## }
  ## else{
  pos <- sindex(jump.times=time,eval.times=obstimes)
  ## loo2 <- do.call("rbind",lapply(1:N,function(k){
  ## Dk <- D
  ## if (status[k]==1) Dk[pos[k]] <- Dk[pos[k]]-1
  ## Yk <- Y-c(rep(1,pos[k]),rep(0,NU-pos[k]))
  ## cumprod(1-Dk/Yk)}))
  ## }
  ## out <- loo
  if (!missing(times)){
    found <- sindex(jump.times=time,eval.times=times)+1
    if (lag==FALSE)
      out <- cbind(1,out)[,found,drop=TRUE]
    else
      out <- cbind(1,cbind(1,out))[,found,drop=TRUE]
  }
  out
}
#' @export 
leaveOneOut.competing.risks <- function(object,times,cause,...){
  stopifnot(object$covariate.type==1)
  mr <- object$model.response
  states <- attr(mr,"states")
  if (missing(cause)) {
    C <- 1
    cause <- states[1]
  }
  else{
    C <- match(cause,states,nomatch=0)
    if (length(C)>1 || C==0) stop("Cause must match exactly one of the names of object$n.event.")
  }
  D <- object$n.event[[C]]
  #  it is sufficient to consider time points where events occur
  time <- object$time[D>0]
  Y <- object$n.risk[D>0]
  sFit <- prodlim(Hist(time,status)~1,data=data.frame(unclass(mr)))
  S <- sFit$surv[D>0]
  D <- D[D>0]
  lagSk <- leaveOneOut.survival(sFit,times=time,lag=1)
  NU <- length(time)
  obstimes <- mr[,"time"]
  status <- mr[,"status"]
  E <- getEvent(mr)
  N <- length(obstimes)
  ## idea: see leaveOneOut.survival
  ## browser()
  ## if (useC==TRUE){
  ## print(cbind(time=time,Y=Y,D=D))
  loo <- .C("loo_comprisk",
            Y = as.double(Y),
            D=as.double(D),
            time=as.double(time),
            obsT=as.double(obstimes),
            status=as.double(status*(E==cause)),
            lagSurv=as.double(lagSk),
            F=double(NU*N),
            N=as.integer(N),
            NT=as.integer(NU),
            PACKAGE="prodlim")$F
  out <- matrix(loo,nrow=N,ncol=NU,byrow=FALSE)
  ## browser()
  ## }
  ## else{
  ## pos <- sindex(jump.times=time,eval.times=obstimes)
  ## loo <- do.call("rbind",lapply(1:N,function(k){
  ## Dk <- D
  ## if (status[k]==1 && E[k]==cause) Dk[pos[k]] <- Dk[pos[k]]-1
  ## Yk <- Y-c(rep(1,pos[k]),rep(0,NU-pos[k]))
  ## Sk <- as.numeric(lagSk[k,,drop=TRUE])
  ## Hk <- Dk/Yk
  ## Fk <- cumsum(Sk*Hk)
  ## Fk
  ## }))
  ## out <- loo
  ## }
  if (!missing(times)){
    found <- sindex(jump.times=time,eval.times=times)+1
    out <- cbind(0,out)[,found,drop=TRUE]
  }
  out
}

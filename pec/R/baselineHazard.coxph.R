baselineHazard.coxph <- function(object,x,y,times=NULL){
  stopifnot(class(object)=="coxph")
  if (is.null(object$x)) stop("You have to say `x=TRUE' in the call to coxph")
  if (!is.null(object$strata)){
    yList <- split(object$y,object$strata)
    xList <- split(object$x,object$strata)
    return(lapply(1:length(yList),function(s){
      object$x <- cbind(xList[[s]])
      object$y <- yList[[s]]
      object$strata <- NULL
      baselineHazard.coxph(object=object,times=times)
    }))}
  ## browser()
  beta <- coef(object)
  if (missing(x)) x <- object$x
  if (missing(y)) y <- object$y
  ## elp <- exp(apply(x,1,function(y)sum(y*beta)))
  elp <- exp(x%*%beta)
  response <- unclass(object$y)
  time <- response[,1,drop=TRUE]
  status <- response[,2,drop=TRUE]
  jumptimes <- sort(time[status!=0])
  dNn <- table(jumptimes)
  S0 <- .C("SNull",
           time=as.double(time),
           jumptimes=as.double(jumptimes),
           elp=as.double(elp),
           S=double(length(jumptimes)),
           N=as.integer(length(time)),
           NJ=as.integer(length(jumptimes)),
           ## DUP=FALSE,
           package="pec")$S
  ## S01 <- sapply(jumptimes,function(s){
  ## sum(elp * (time>=s))
  ## })
  Lambda <-   cumsum((1/S0) * dNn )
  if (is.null(times)){
    data.frame(time=as.vector(unique(jumptimes)),cumhazard=as.vector(Lambda))
  }
  else{
    data.frame(time=times,cumhazard=c(0,Lambda)[1+prodlim::sindex(jump.times=jumptimes,eval.times=times)])
  }
}


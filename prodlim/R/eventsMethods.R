events <- function(object,...){
  UseMethod("events",object)
}


events.prodlim <- function(object){
  events.Hist(object$model.response)
}


events.Hist <- function(object,...){
  model <- attr(object,"model")
  cens.code <- attr(object,"cens.code")
  states <- attr(object,"states")
  if (model=="survival"){
    factor(object[,"status",drop=TRUE],levels=c(cens.code,states),labels=c("unknown",states))
  }
  else{
    if (model=="competing.risks"){
      D <- object[,"status",drop=TRUE]
      cens.type <- attr(object,"cens.type")
      E <- object[,"event",drop=TRUE]
      names(E) <- NULL
      stupid.stupid.factor.levels <- as.integer(factor(c(states,cens.code),levels=c(states,cens.code)))
      ## sorted.stupid.stupid.factor.levels <- c(stupid.stupid.factor.levels[-1],stupid.stupid.factor.levels[1])
      ## events <- factor(E,levels=sorted.stupid.stupid.factor.levels,labels=c(states,"unknown"))
      events <- factor(E,levels=stupid.stupid.factor.levels,labels=c(states,"unknown"))
      events
    }
    else stop("No event.Hist function for multi.state models")
  }
}




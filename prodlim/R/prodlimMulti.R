prodlimMulti <- function(response,size.strata,N,NU,cotype,force.multistate){
  ##   original function by Matthias `Wang' Wangler 
  is.event <- response[,"status"]!=0
  if (force.multistate==TRUE){
    to    <- response[,"status"]
    from <- rep(0,length(to))
  }
  else{
    to    <- response[,"event"]
    from  <- response[,"from"]
  }
  state.names <- unique(c(from, to[response[,"status"]!=0]))
  ns <- length(state.names)
  cens <- FALSE
  if(length(to[is.event])>0) cens <- TRUE
  from <- as.integer(factor(from,levels=state.names)) - 1
  from <- as.numeric(from)
  to[is.event] <- as.integer(factor(to[is.event], levels=state.names)) - 1
  to[!is.event] <- ns
  to <- as.numeric(to)
  states <- sort(unique(c(from, to[is.event])))   
  ## possible transitions
  tra <- unique(cbind(from[is.event], to[is.event]))
  sorted <- order(tra[,1],tra[,2])    
  tra <- matrix(tra[sorted,], ncol=2)
  tra <- cbind(0:(length(tra[,1])-1),tra)
  colnames(tra) <- c("row","from", "to")
  ntra <- nrow(tra)
  trow <- match(paste(from,to), paste(tra[,"from"],tra[,"to"]), nomatch=0) - 1
  cens.in <- sort(unique(from[!is.event]))
  nci <- length(cens.in)
  cpos <- match(paste(from,to), paste(cens.in, ns), nomatch = 0) - 1
  ## start distribution (all are starting in state 0 !!!)
  if( cotype > 1 ) {
    #    nr.start <- table(from,co$covariates$strata[,1])[1,]
    nr.start <- size.strata  ## WANG???
  }
  else{nr.start <- length(from[from==0])}
  fit <- .C("prodlim_multistates",
            as.integer(N),
            as.integer(ns),
            as.integer(length(is.event)),
            as.integer(size.strata),
            as.integer(ntra),
            as.integer(tra[,"from"]),
            as.integer(tra[,"to"]),
            as.integer(trow),
            as.integer(nci),
            as.integer(cens.in),
            as.integer(cpos),
            as.double(response[,"time"]),
            as.integer(response[,"status"]),
            as.integer(nr.start),
            time=double(N),
            hazard=double(N*ns*ns),
            prob=double(N*ns*ns),
            nevent=integer(N*ns*ns),
            ncens=integer(N*ns),
            nrisk=integer(N*ns),
            first.strata=integer(NU),
            ntimes.strata=integer(NU),
            PACKAGE="prodlim")
  tra[,"from"] <- state.names[tra[,"from"]+1]
  tra[,"to"] <- state.names[tra[,"to"]+1]
  cens.in <- state.names[cens.in+1]
  NT <- sum(fit$ntimes.strata)
  res <- list("time"=fit$time[1:NT],"hazard"=fit$hazard[1:(NT*ns*ns)],"prob"=fit$prob[1:(NT*ns*ns)],"nevent"=fit$nevent[1:(NT*ns*ns)],"ncens"=fit$ncens[1:(NT*ns)],"nrisk"=nrisk <- fit$nrisk[1:(NT*ns)],"first.strata"=fit$first.strata,"size.strata"=fit$ntimes.strata,"uniquetrans"=tra,"cens.in"=cens.in,"states"=states,"state.names"=state.names,"model"="multi.states")
  res
}

#' Summary of data
#'
#' Compute statistical summaries (mean, quantile, variance, survival rate,...)
#' 
#' See http://simulx.webpopix.org/stamlx for more details.      
#' @param r a data frame
#' @param FUN a string, or a vector of strings, with the name of the functions to apply to the result of the simulation
#' @param probs a vector of quantiles  between 0 and 1. Only used when "quantile" has been defined in \code{FUN} 
#' @param surv.time a scalar or a vector of times. Only used when "event" has been defined in \code{type} 
#' 
#' @return A data frame. 
#' @examples
#' \dontrun{
#' modelPK <- inlineModel("
#' [LONGITUDINAL] 
#' input={V,Cl,alpha, beta,b}
#' 
#' EQUATION:
#' C = pkmodel(V, Cl)
#' h = alpha*exp(beta*C)
#' g = b*C
#' 
#' DEFINITION:
#' y = {distribution=normal, prediction=C, sd=g}
#' e = {type=event, maxEventNumber=1, rightCensoringTime=30, hazard=h}
#' 
#' [INDIVIDUAL]
#' input={V_pop,Cl_pop,omega_V,omega_Cl}
#' 
#' DEFINITION:
#' V     = {distribution=lognormal,   prediction=V_pop,    sd=omega_V}
#' Cl    = {distribution=lognormal,   prediction=Cl_pop,   sd=omega_Cl}
#' ")
#' 
#' adm  <- list(amount=100, time=0)
#' p <- c(V_pop=10, Cl_pop=1, omega_V=0.2, omega_Cl=0.2, alpha=0.02, beta=0.1, b=0.1)
#' out.y <- list(name=c('y'), time=seq(0,to=25,by=5))
#' out.e <- list(name=c('e'), time=0)
#' out.p <- c("V", "Cl")
#' out   <- list(out.y, out.e, out.p)
#' g <- list(size=100, level='individual')
#' res1 <- simulx(model=modelPK, treatment=adm, parameter=p, output=out, group=g)
#' 
#' statmlx(res1$parameter, FUN = "mean", probs = c(0.05, 0.5, 0.95))
#' statmlx(res1$parameter, FUN = "quantile", probs = c(0.05, 0.5, 0.95))
#' statmlx(res1$parameter, FUN = c("sd", "quantile"), probs = c(0.05, 0.95))
#' statmlx(res1$y, FUN = c("mean", "sd", "quantile"), probs = c(0.05, 0.95))
#' statmlx(res1$e, surv.time=c(10,20))
#' 
#' res2 <- simulx(model=modelPK, treatment=adm, parameter=p, output=out, group=g, nrep=3)
#' statmlx(res2$parameter, FUN = c("sd", "quantile"), probs = c(0.05, 0.95))
#' statmlx(res2$y, FUN = c("mean", "sd", "quantile"), probs = c(0.05, 0.95))
#' statmlx(res2$e, surv.time=c(10,20,30))
#' }
#' @importFrom stats aggregate
#' @export


statmlx <- function(r, FUN="mean", probs=c(0.05, 0.5, 0.95), surv.time=NULL)
{
  if (!is.null(surv.time))
    type="event"
  else
    type="continuous"
  
  if (any(!(FUN %in% c("mean","sd","median","var","quantile"))))
    stop("\n\n possible values for 'FUN' are {'mean','sd','median','var','quantile'} ")
  # if (any(!(type %in% c('continuous','event'))))
  # stop("\n\n possible values for 'type' are {'continuous','event'} ")
  # if (any(type=="survival") && is.null(time)) 
  #   stop("\n\n a time value, or a vector of times, should be provided when 'survival' is defined as 'FUN'")
  
  r <- r[!is.na(names(r))]
  
  if (type=="event"){
    r <- r[r$time>0,]
    it <- which(names(r) == "time")
    ie <- which(!(names(r) %in% c("pop","rep","id","time","group")))
    ig <- r[-c(it,ie)]
    uid <-  cumsum(!duplicated(ig)) 
    y <- cbind(uid,r)
    
    rrr <- list(uid)
    ev.nb <- aggregate(r[ie], by=rrr, "sum")
    names(ev.nb)[2] <- "nbEv"
    rev <- cbind(unique(ig), ev.nb[2]) 
    nev <- names(rev)
    ev.min <- aggregate(r[it], by=rrr, "min")
    if (!is.null(surv.time)){
      for (k in (1:length(surv.time))){
        survk <- (ev.min$time >= surv.time[k])*1
        rev   <- cbind(rev, survk)
      }
      nev <- c(nev,paste0("S",surv.time))
    }
    names(rev) <- nev
    # ev.t1 <- ev.min
    # ev.t1[ev.nb==0] <- NA
    
    r <- rev
  }
  
  list.FUN <- c("mean","sd","median","var") 
  l.FUN <- length(FUN)
  n.FUN <- NULL
  for (j in (1:l.FUN)){
    if (FUN[j] %in% list.FUN)
      n.FUN <- c(n.FUN,FUN[j])
    else 
      n.FUN <- c(n.FUN,paste0("p",probs*100))
  }
  
  rrr <- list()
  nr <- NULL
  if (!is.null(r$pop)){
    rrr[[length(rrr)+1]] <- r$pop 
    nr <- c(nr, "pop")
  }
  if (!is.null(r$rep)){
    rrr[[length(rrr)+1]] <- r$rep 
    nr <- c(nr, "rep")
  }
  if (!is.null(r$group)){
    rrr[[length(rrr)+1]] <- r$group 
    nr <- c(nr, "group")
  }  
  if (!is.null(r$time)){
    rrr[[length(rrr)+1]] <- r$time 
    nr <- c(nr, "time")
  }
  r$pop <- r$rep <- r$group <- r$time <- r$id <- NULL
  
  sr <- NULL
  for (j in (1:l.FUN)){
    statj <- FUN[j]
    if (statj %in% list.FUN) {
      if (length(rrr)>0)
        X <- aggregate(r, by=rrr, statj)
      else
        X <- as.data.frame(t(sapply(r,statj)))
      n.srj <- paste0(names(r),".",statj)
    } else {
      if (length(rrr)>0){
        aX <- aggregate(r, by=rrr, statj, probs=probs)
        nX <- length(aX)
        X <- NULL
        for (k in (1:nX))
          X <- cbind(X, aX[[k]])
        X <- as.data.frame(X)
      } else {
        aX <- sapply(r,statj, probs=probs)
        X <- as.data.frame(t(as.vector(aX)))
      }
      aa <- apply(expand.grid(names(r),paste0("p",probs*100)), 1, paste, collapse=".")  
      n.srj <- t(matrix(aa,nrow=ncol(r)))
    }
    names(X) <- c(nr,n.srj)
    if (j==1)
      sr <- X
    else 
      sr <- merge(sr, X)
  }
  lo <- NULL
  if (!is.null(sr$pop)) 
  {
    lo <- c(lo, "pop")
    sr$pop <- as.factor(sr$pop)
  }
  if (!is.null(sr$rep))
  {
    lo <- c(lo, "rep")
    sr$rep <- as.factor(sr$rep)
  }
  if (!is.null(sr$group)) 
  {
    lo <- c(lo, "group")
    sr$group <- as.factor(sr$group)
  }
  if (!is.null(sr$time)) 
    lo <- c(lo, "time")
  if (!is.null(lo))
  {
    lo <- paste(lo,collapse=",")
    eval(parse(text=paste0("sr <- sr[with(sr, order(",lo,")), ]")))
  }
  return(sr)
}




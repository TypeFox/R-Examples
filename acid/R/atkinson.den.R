atkinson.den <-
function (incs,dens,epsilon = 1,pm0=NA,lower=NULL,upper=NULL,zero.approx=NULL){ 
  if(length(incs)!=length(dens)) print("incs and dens must have the same length!")
  incs       <- as.numeric(incs)
  if(is.null(lower)) lower<-min(incs)
  if(is.null(upper)) upper<-max(incs)
  if(min(incs)<lower) print("WARNING: min(incs)<lower int limit!")
  if(!is.na(pm0)&incs[1]==0) incs<-incs[-1]
  n.inc      <- length(incs)
  incs.diff  <- diff(incs)
  av.incs    <- (incs[-1]+incs[-n.inc])/2
  dens       <- as.numeric(dens)
  if(!is.na(pm0)&dens[1]==pm0) dens<-dens[-1]
  av.dens    <- (dens[-1]+dens[-n.inc])/2
  if(is.null(zero.approx)){
    fi.inner <- av.dens*incs.diff
    fi.outer <- c(dens[1]/2*(incs[1]-lower),dens[n.inc]/2*(upper-incs[n.inc]))
    m.inner.int<- t(av.incs)%*%(fi.inner)
    m.outer.int<- incs[1]*fi.outer[1]+incs[n.inc]*fi.outer[2]
    m          <- m.inner.int+m.outer.int 
    if (is.null(epsilon)) epsilon <- 1
    if (epsilon == 1) {
      if(!is.na(pm0)){ A <- 1 - (exp(weighted.mean(log(c(0,incs[1],av.incs,incs[n.inc])),c(pm0,fi.outer[1],fi.inner,fi.outer[2])))/m)
      }else A <- 1 - (exp(weighted.mean(log(c(incs[1],av.incs,incs[n.inc])),c(fi.outer[1],fi.inner,fi.outer[2])))/m)
    }else {
      av.A <- (av.incs/m)^(1 - epsilon)
      A.innter.int <- t(av.A)%*%fi.inner
      A.outer.int  <- (incs[1]/m)^(1 - epsilon)*fi.outer[1]+(incs[n.inc]/m)^(1 - epsilon)*fi.outer[2]
      A<- 1- (A.innter.int+A.outer.int)^(1/(1-epsilon))
    }
  }else{
    fi.inner <- av.dens*incs.diff
    fi.outer <- c(dens[1]/2*(incs[1]-lower),dens[n.inc]/2*(upper-incs[n.inc]),pm0)
    m.inner.int<- t(av.incs)%*%(fi.inner)
    m.outer.int<- incs[1]*fi.outer[1]+incs[n.inc]*fi.outer[2]+pm0*zero.approx
    m          <- m.inner.int+m.outer.int 
    if (is.null(epsilon)) epsilon <- 1
    if (epsilon == 1) {
      A <- 1 - (exp(weighted.mean(log(c(zero.approx,incs[1],av.incs,incs[n.inc])),c(pm0,fi.outer[1],fi.inner,fi.outer[2])))/m)
    }else {
      av.A <- (av.incs/m)^(1 - epsilon)
      A.innter.int <- t(av.A)%*%fi.inner
      A.outer.int  <- (incs[1]/m)^(1 - epsilon)*fi.outer[1]+(incs[n.inc]/m)^(1 - epsilon)*fi.outer[2]+(zero.approx/m)^(1 - epsilon)*fi.outer[3]
      A<- 1- (A.innter.int+A.outer.int)^(1/(1-epsilon))
    }
  }
  list(AIM=A,epsilon=epsilon,mean=m,pm0=pm0,lower=lower,upper=upper,zero.approx=zero.approx)
}

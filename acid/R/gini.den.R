gini.den <-
function (incs,dens,pm0=NA,lower=NULL,upper=NULL){ 
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
  fi.inner <- av.dens*incs.diff
  fi.outer <- c(dens[1]/2*(incs[1]-lower),dens[n.inc]/2*(upper-incs[n.inc]),pm0)
  
  if(!is.na(pm0)&incs[1]==0) Gini<-weighted.gini(x=c(0,lower,av.incs,upper),w=c(fi.outer[3],fi.outer[1],fi.inner,fi.outer[2]))$Gini
  else Gini<-weighted.gini(x=c(lower,av.incs,upper),w=c(fi.outer[1],fi.inner,fi.outer[2]))$Gini
  list(Gini=Gini,pm0=pm0,lower=lower,upper=upper)
}

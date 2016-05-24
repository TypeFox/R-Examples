################### scalecheck #######################
scalecheck<-function(par, lower=lower, upper=upper, bdmsk=NULL, dowarn=TRUE){
   # a function to check the initial parameters and bounds for inputs to optimization codes
   # Arguments:
   #   par -- starting parameters supplied 
   #    lower, upper -- lower and upper bounds supplied
   #    bdmsk -- bounds and masks indicator
   #
   # Returns:
   #   list(lpratio, lbratio) -- the log of the ratio of largest to smallest parameters
   #      and bounds intervals (upper-lower) in absolute value (ignoring Inf, NULL, NA)
   ######################################
   if (is.null(par)) { stop("Null parameter vector") }
   npar<-length(par)
   if (is.null(bdmsk)) bdmsk <- rep(1,npar) # ensure bdmsk defined
   if (is.null(lower)) { 
      if (dowarn) warning("Null lower bounds vector")
      lower<-rep(-Inf,npar)
   }
   if (is.null(upper)) { 
      if (dowarn) warning("Null upper bounds vector")
      upper<-rep(Inf,npar)
   }
   # Ignore masks
   lower<-lower[which(bdmsk==1)]
   upper<-upper[which(bdmsk==1)]
   newpar<-par[which(bdmsk==1)]
   # Now check the parameters
   newpar<-abs(par[which(is.finite(par))])
   logpar<-log10(newpar[which(newpar>0)]) # Change 20100711
   newlower<-abs(lower[which(is.finite(lower))])
   loglower<-log10(newlower[which(newlower>0)]) # Change 20100711
   newupper<-abs(upper[which(is.finite(upper))])
   logupper<-log10(newupper[which(newupper>0)]) # Change 20100711
   bddiff<-upper-lower
   bddiff<-bddiff[which(is.finite(bddiff))]
   lbd<-log10(bddiff[which(bddiff>0)]) # Change 20100711
   lpratio<-max(logpar) - min(logpar)
   if (length(lbd) > 0) {
      lbratio<-max(lbd)-min(lbd)
   } else { 
      lbratio<-NA
   }
   ratios<-list(lpratio=lpratio,lbratio=lbratio) # return(ratios)
}
################### end scalecheck #######################

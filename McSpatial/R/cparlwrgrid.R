cparlwrgrid <- function(form,nonpar,window=0,bandwidth=0,kern="tcub",method="gcv",print=TRUE,distance="Mahal",targetobs=NULL,data=NULL) {

  nw = length(window)
  nb = length(bandwidth)
  if ((nw==1)&(nb==1)) { cat("Must provide a vector of window or bandwidth values","\n") }
  if ((nw>1)&(nb>1))   { cat("Both window and bandwidth vectors specified; will use window vector","\n") }

  minval = 9999999
  minterm = 0
  icross = ifelse(method=="gcv",FALSE,TRUE)

  if (nw>1) {
    for (iw in window) {
      fit1 <- cparlwr(form,nonpar,window=iw,bandwidth=0,kern=kern,distance=distance,targetobs=targetobs,data=data)
      hval = ifelse(icross==TRUE,fit1$cv,fit1$gcv)
      if (print==TRUE) {print(c(iw,hval))}
      if (hval<minval) {
        minh = iw
        minval = hval
        outfit <- fit1 
      }
     }
   }

  if (nb>1) {
    for (ib in bandwidth) {
      fit1 <- cparlwr(form,nonpar,window=0,bandwidth=ib,kern=kern,distance=distance,targetobs=targetobs,data=data)
      hval = ifelse(icross==TRUE,fit1$cv,fit1$gcv)
      if (print==TRUE) {print(c(ib,hval))}
      if (hval<minval) {
        minh = ib
        minval = hval
        outfit <- fit1 
      }
    }
   }

  out <- list(outfit$target,outfit$ytarget,outfit$xcoef.target,outfit$xcoef.target.se,
      outfit$yhat,outfit$xcoef,outfit$xcoef.se,
      outfit$df1,outfit$df2,outfit$sig2,outfit$cv,outfit$gcv,outfit$infl,minh)
  names(out) <- c("target","ytarget","xcoef.target","xcoef.target.se",
    "yhat","xcoef","xcoef.se","df1","df2","sig2","cv","gcv","infl","minh")


  cat("\n","h = ",minh,"Function Value = ",minval,"\n") 
  return(out)
}


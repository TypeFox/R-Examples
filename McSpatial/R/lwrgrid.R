lwrgrid <- function(form,window=0,bandwidth=0,kern="tcub",method="gcv",print=TRUE,distance="Mahal",target=NULL,data=NULL) {

  nw = length(window)
  nb = length(bandwidth)
  if ((nw==1)&(nb==1)) { cat("Must provide a vector of window or bandwidth values","\n") }
  if ((nw>1)&(nb>1))   { cat("Both window and bandwidth vectors specified; will use window vector","\n") }

  minval = 9999999
  minterm = 0
  icross = ifelse(method=="gcv",FALSE,TRUE)

  if (nw>1) {
    for (iw in window) {
      fit1 <- lwr(form,window=iw,bandwidth=0,kern=kern,distance=distance,target=target,data=data)
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
      fit1 <- lwr(form,window=0,bandwidth=ib,kern=kern,distance=distance,target=target,data=data)
      hval = ifelse(icross==TRUE,fit1$cv,fit1$gcv)
      if (print==TRUE) {print(c(ib,hval))}
      if (hval<minval) {
        minh = ib
        minval = hval
        outfit <- fit1 
      }
    }
   }
  
  out <- list(outfit$target,outfit$ytarget,outfit$dtarget1,outfit$dtarget2,outfit$ytarget.se,outfit$dtarget1.se,outfit$dtarget2.se,
      outfit$yhat,outfit$dhat1,outfit$dhat2,outfit$yhat.se,outfit$dhat1.se,outfit$dhat2.se,outfit$df1,outfit$df2,
      outfit$sig2,outfit$cv,outfit$gcv,outfit$infl,minh)
  names(out) <- c("target","ytarget","dtarget1","dtarget2","ytarget.se","dtarget1.se","dtarget2.se",
                  "yhat","dhat1","dhat2","yhat.se","dhat1.se","dhat2.se","df1","df2","sig2","cv","gcv","infl","minh")


  cat("\n","h = ",minh,"Function Value = ",minval,"\n") 
  return(out)
}


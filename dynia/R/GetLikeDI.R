GetLikeDI <-
function(delta,z,T,xint=NA,itype=c("step","pulse"),...){

  
  n<-length(z)
  n2 <- n - T+1
  itype=match.arg(itype)
  ###Check the type of intervention##
  if (itype=="pulse")
  {xi <- c(rep(0, n - n2), delta^(0:(n2 - 1)))}
    else if (itype=="step") 
  {xi <- c(rep(0, n - n2), cumsum(delta^(0:(n2 - 1))))}

  ###Check the existance of covariate###
  if (exists("xint",mode="numeric")==F) xreg2<-xi
     else 
     { 
       xreg2 <- matrix(c(xi, xint), ncol = 2)
       dimnames(xreg2)[[2]] <- c("omega", "xint")
     }
   out<-arima(z,...,xreg=xreg2)
   as.numeric(logLik(out))
}

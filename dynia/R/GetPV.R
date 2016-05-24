GetPV <-
function(delta0,z,T,xint=NA,itype=c("step","pulse"),...){
  
  n<-length(z)
  n2 <- n - T
  itype=match.arg(itype)
  if (exists("xint",mode="numeric")==F) 
  {
    LLRestricted <- GetLikeDI(delta=delta0,z=z,T=T,itype=itype,...)
    ans <- optimize(f = GetLikeDI, interval = c(0.1, 2), maximum = TRUE,z=z,T=T
                    ,itype=itype,...) 
    LLFull <- ans$objective

   }
  else 
  { 
    LLRestricted <- GetLikeDI(delta=delta0,z=z,T=T,xint=xint,itype=itype,...)
    ans <- optimize(f = GetLikeDI, interval = c(0.1, 2), maximum = TRUE,z=z,T=T
                    ,xint=xint,itype=itype,...)
    LLFull <- ans$objective

  }


  X <- 2 * (LLFull - LLRestricted)
  pval <- 1 - pchisq(X, 1)
 return(pval)
  #list(delta = ans$maximum, pvalue = pval, Int.Model=out$Intervention)
}

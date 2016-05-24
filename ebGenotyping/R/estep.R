estep <-
function(mu,delta,pm1,p0,dat,cvg){
  pp1 <- 1-pm1-p0
  a <- outer(delta,mu,"+")
  zm1 <- exp(dat*a-cvg*log((1+exp(a))/2))*pm1
  zp1 <- exp((cvg-dat)*a-cvg*log((1+exp(a))/2))*pp1
  z0 <- p0
  inf.index <- which(is.infinite(zm1))
    if (length(inf.index)>0){
      zm1[inf.index] <- 1e+100*pm1
    }
  inf.index <- which(is.infinite(zp1))
    if (length(inf.index)>0){
      zp1[inf.index] <- 1e+100*pp1
    }  
  zsum <- zm1+z0+zp1
  zm1 <- zm1/zsum
  z0 <- z0/zsum
  zp1 <- zp1/zsum
  return(list(zRR=zm1,zRV=z0,zVV=zp1))
}

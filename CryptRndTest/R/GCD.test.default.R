GCD.test.default=function(x,B=32,KS=TRUE,CSQ=TRUE,AD=TRUE,JB=TRUE,test.k=TRUE,test.g=TRUE,mu,sd,alpha=0.05){
  options(warn=-1)
  check(test=4,x=x,B=B,alpha=alpha,KS=KS,CSQ=CSQ,AD=AD,JB=JB,test.k=test.k,test.g=test.g,mu=mu,sd=sd)
  
  res.tbl=GCD.test.main(x=x,B=B,KS=KS,CSQ=CSQ,AD=AD,JB=JB,test.k=test.k,test.g=test.g,mu=mu,sd=sd,alpha=alpha)
  res.tbl$call = match.call()
  class(res.tbl) = c("GCD.test","CryptRndTest")
  res.tbl
}
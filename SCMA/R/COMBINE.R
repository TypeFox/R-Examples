combine<-function(method, pvalues=read.table(file.choose(new=FALSE))){

  n<-nrow(pvalues)
  
  if(method=="x"){
    Q<-prod(pvalues)
    mult<-numeric()
    for(k in 0:(n-1)){
      mult<-c(mult,((-log(Q))^k)/factorial(k))
    }
    lambda<-Q*sum(mult)
  }
  
  if(method=="+"){
    S<-sum(pvalues)
    add<-numeric()
    for(k in 0:floor(S)){
      add<-c(add,((-1)^k)*(choose(n,k))*((S-k)^n))
    }
    lambda<-sum(add)/factorial(n)
  }
  
  return(lambda)

} 
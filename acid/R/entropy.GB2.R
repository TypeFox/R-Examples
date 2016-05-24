entropy.GB2 <-
function(b,a,p,q,alpha=NULL,ylim=c(0,1000000),zeroapprox=0.01){ ## - see Jenkins (2009): Distributionally-Sensitive Inequality Indices ... and Cowell (2000), pp.109-110 for details.
  if (is.null(alpha)) alpha <- 1
  entropy<-rep(NA,length(b))
  for(j in 1:length(b)){ 
    if (alpha == 0){
      entropyj  <- lgamma(p[j]+1/a[j])+lgamma(q[j]-1/a[j])-lgamma(p[j])-lgamma(q[j])-digamma(p[j])/a+digamma(q[j])/a
    }else if (alpha == 1){ 
      entropyj  <- digamma(p+1/a)/a-digamma(q-1/a)/a-lgamma(p+1/a)-lgamma(q-1/a)+lgamma(p)+lgamma(q)
    }else if (alpha == -1){
      entropyj  <- -1/2+ ( gamma(p[j]-1/a[j])*gamma(q[j]+1/a[j])*gamma(p[j]+1/a[j])*gamma(q[j]-1/a[j]) )/(2*(gamma(p[j]))^2*(gamma(q[j]))^2)
    }else if (alpha == 2){
      entropyj  <- -1/2+ (gamma(p[j])*gamma(q[j])*gamma(p[j]+2/a[j])*gamma(q[j]-2/a[j]))/(2*(gamma(p[j]+1/a))^2*(gamma(q[j]-1/a))^2)
    }else{
      den       <- dGB2(y,b[j],a[j],p[j],q[j])
      meanGB2   <- b[j]*beta(p[j]+1/a[j],q[j]-1/a[j])/beta(p[j],q[j])
      yi        <- y[-length(y)]+(max(y)-min(y))/(length(y)-1)/2
      deni      <- dGB2(yi,b[j],a[j],p[j],q[j])
      y         <- seq(ylim[1],ylim[2],length.out=10001)+zeroapprox
      int_alpha <- t(yi^alpha)%*%deni*(max(y)-min(y))/(length(y)-1)
      entropyj  <- (int_alpha/meanGB2^alpha-1)/(alpha*(alpha-1))
    }
    entropy[j]<-entropyj
  }  
  return(entropy)
}

seqadd.confint=function(L, alpha, presample, addsample){
  n0=length(presample)
  s0=sd(presample)
  z=(L/2)^2/qt(1-alpha/2,n0-1)^2
  n=max(n0+1,ceiling(s0^2/z))
  a=z*(s0^2)
  if (n<1/a) n=ceiling(1/a)
  g=function(x){((n/n0)/(n-n0))*x^2-(2/(n-n0))*x+1/(n-n0)-z*s0^2}
  b=seq(-1,1,0.05)
  lb=length(b)
  k=max(which(g(b)*g(1)<0))
  a0=uniroot(g,c(b[k],1))$root
  a1=(1-a0)/(n-n0)
  ln=a0*mean(presample)+a1*sum(addsample)
  list("Significance level"=alpha,
       "Length of confidence interval"=L,
       "Length of presample"=n0,
       "Number of additional observations"=n-n0,
       "Total number of observations"=n,
       "confidence interval"=c(ln-L/2,ln+L/2))
}


seqadd.normal=function(L, alpha, presample){
  n0=length(presample)
  s0=sd(presample)
  z=(L/2)^2/qt(1-alpha/2,n0-1)^2
  n=max(n0+1,ceiling(s0^2/z))
  list("Significance level"=alpha,
       "Length of confidence interval"=L,
       "Length of presample"=n0,
       "Standard deviation of presample"=s0,
       "Additional number of observations"=n-n0)
}

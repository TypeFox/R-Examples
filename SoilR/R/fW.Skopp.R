#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.Skopp<- structure(
  function #Effects of moisture on decomposition rates according to the function proposed by Skopp et al. 1990
  ### Calculates the effects of relative soil water content on decomposition rates.
  ##references<< J. Skopp, M. D. Jawson, and J. W. Doran. 1990.
  ##Steady-state aerobic microbial activity as a function of soil water content. Soil Sci. Soc. Am. J., 54(6):1619-1625

(  rwc, ##<< relative water content
   alpha=2, ##<< Empirical parameter
   beta=2, ##<< Empirical parameter
   f=1.3, ##<< Empirical parameter
   g=0.8 ##<< Empirical parameter
   
)
  {
  p=pmin(alpha*rwc^f, beta*(1-rwc)^g)
  return(p)
}
,
ex=function(){
  th=seq(0,1,0.01)
  xi=fW.Skopp(rwc=th)
  plot(th,xi,type="l",main="Effects of soil water content on decomposition rates",
       xlab="Relative water content",ylab=expression(xi))
}
)

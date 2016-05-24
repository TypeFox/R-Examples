#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.Moyano<- structure(
  function #Effects of moisture on decomposition rates according to the function proposed by Moyano et al. (2013)
  ### Calculates the effects of water content on decomposition rates.
  ##references<< F. E. Moyano, S. Manzoni, C. Chenu. 2013
  ##Responses of soil heterotrophic respiration to moisture availability: An exploration of processes and models. 
  ##Soil Biology and Biochemistry, Volume 59, April 2013, Pages 72-85
  
  (theta,     ##<< A scalar or vector containing values of volumetric soil water content.
   a=3.11, ##<< Empirical parameter
   b=2.42 ##<< Empirical parameter
  )
  {
    a*theta-b*theta^2
  }
  ,
  ex=function(){
    th=seq(0,1,0.01)
    xi=fW.Moyano(theta=th)
    plot(th,xi,type="l",main="Effects of soil water content on decomposition rates",
         xlab="Volumetric soil water content (cm3 cm-3)",ylab=expression(xi))
  }
)

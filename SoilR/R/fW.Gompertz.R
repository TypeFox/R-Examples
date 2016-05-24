#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.Gompertz<- structure(
  function #Effects of moisture on decomposition rates according to the Gompertz function
  ### Calculates the effects of water content on decomposition rates.
  ##references<< I. Janssens, S. Dore, D. Epron, H. Lankreijer, N. Buchmann, B. Longdoz, J. Brossaud, L. Montagnani. 2003.
  ##Climatic Influences on Seasonal and Spatial Differences in Soil CO2 Efflux. In Valentini, R. (Ed.)
  ##Fluxes of Carbon, Water and Energy of European Forests. pp 235-253. Springer. 
  
  (theta,     ##<< A scalar or vector containing values of volumetric soil water content.
   a=0.824, ##<< Empirical parameter
   b=0.308 ##<< Empirical parameter
  )
  {
  exp(-exp(a-b*theta*100))
  }
  ,
  ex=function(){
    th=seq(0,1,0.01)
    xi=fW.Gompertz(theta=th)
    plot(th,xi,type="l",main="Effects of soil water content on decomposition rates",
         xlab="Volumetric soil water content (cm3 cm-3)",ylab=expression(xi))
  }
)

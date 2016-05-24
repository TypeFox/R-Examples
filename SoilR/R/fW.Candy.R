#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.Candy<- structure(
  function #Effects of moisture on decomposition rates according to the Candy model
  ### Calculates the effects of water content and pore volume on decomposition rates.
  ##references<< J. Bauer, M. Herbst, J.A. Huisman, L. Weiherm\"uller, H. Vereecken. 2008.
  ##Sensitivity of simulated soil heterotrophic respiration to temperature and moisture reduction functions.
  ##Geoderma, Volume 145, Issues 1-2, 15 May 2008, Pages 17-27.
  
  (theta,     ##<< A scalar or vector containing values of volumetric soil water content.
   PV   ##<< A scalar or vector containing values of pore volume.
  )
  {
    Mi=theta/PV # Moisture index
    
    fw=ifelse(Mi<=0.5, 4*Mi*(1-Mi),1)
    
    return(fw)
  }
,
  ex=function(){
    th=seq(0,1,0.01)
    xi1=fW.Candy(theta=th,PV=0.4)
    xi2=fW.Candy(theta=th,PV=0.6)
    xi3=fW.Candy(theta=th,PV=0.8)
    plot(th,xi1,type="l",main="Effects of soil water content and pore volume on decomposition rates",
         xlab="Volumetric soil water content (cm3 cm-3)",ylab=expression(xi))
    lines(th,xi2,col=2)
    lines(th,xi3,col=3)
    legend("bottomright",c("Pore volume = 0.4","Pore volume = 0.6", "Pore volume = 0.8"),lty=1,col=1:3)
  }
)

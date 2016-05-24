#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.LandT<-structure(
  function #Effects of temperature on decomposition rates according to a function proposed by Lloyd and Taylor (1994)
    ### Calculates the effects of temperature on decomposition rates according to a function proposed by Lloyd and Taylor (1994).
    ##references<<  Lloyd, J., and J. A. Taylor (1994), On the Temperature Dependence of Soil Respiration, 
    ##Functional Ecology, 8(3), 315-323.
    (Temp     ##<< A scalar or vector containing values of soil temperature for which the effects on decomposition rates are calculated
     )
   {
      exp(308.56*((1/56.02)-(1/((Temp+273)-227.13))))
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Temperature=0:50
      plot(Temperature,fT.LandT(Temperature),type="l",
           ylab="f(T) (unitless)", 
           main="Effects of temperature on decomposition 
           rates according to the Lloyd and Taylor function")
    }
)

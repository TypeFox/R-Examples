#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.Century2<-structure(
  function #Effects of temperature on decomposition rates according the the CENTURY model
    ### Calculates the effects of temperature on decomposition rates according to the CENTURY model.
    ##references<<  Adair, E. C., W. J. Parton, S. J. D. Grosso, W. L. Silver, M. E. Harmon, S. A. Hall, I. C. Burke, and S. C. Hart. 2008. 
    ##Simple three-pool model accurately describes patterns of long-term litter decomposition in diverse climates. Global Change Biology 14:2636-2660.
    (Temp,     ##<< A scalar or vector containing values of temperature for which the effects on decomposition rates are calculated.
     Tmax=45,   ##<< A scalar defining the maximum temperature in degrees C.
     Topt=35   ##<< A scalar defining the optimum temperature for the decomposition process in degrees C.
     )
   {
      3.439423*exp((0.2/2.63)*(1-(((Tmax-Temp)/(Tmax-Topt))^2.63))*((Tmax-Temp)/(Tmax-Topt))^0.2)
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Temperature=0:50
      plot(Temperature,fT.Century2(Temperature),type="l",
           ylab="f(T) (unitless)",
           main="Effects of temperature on decomposition rates according to the Century model")
    }
)

#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.Standcarb<-structure(
  function #Effects of temperature on decomposition rates according to the StandCarb model
    ### Calculates the effects of temperature on decomposition rates according to the StandCarb model.
    ##references<< Harmon, M. E., and J. B. Domingo (2001), A users guide to STANDCARB version 2.0: 
    ##A model to simulate carbon stores in forest stands. Oregon State University, Corvallis.

    (Temp,      ##<< A scalar or vector containing values of temperature for which the effects on decomposition rates are calculated.
     Topt=45,   ##<< A scalar representing the optimum temperature for decomposition.
     Tlag=4,    ##<< A scalar that determines the lag of the response curve.
     Tshape=15, ##<< A scalar that determines the shape of the response curve.
     Q10=2      ##<< A scalar. Temperature coefficient Q10.
     )
   {
     exp(-1*(Temp/(Topt+Tlag))^Tshape)*Q10^((Temp-10)/10)
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Temperature=0:50
      plot(Temperature,fT.Standcarb(Temperature),type="l",ylab="f(T) (unitless)", 
           main="Effects of temperature on decomposition rates according to the StandCarb model")
    }
)

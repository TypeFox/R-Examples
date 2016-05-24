#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.Demeter<-structure(
  function #Effects of moisture on decomposition rates according to the DEMETER model
    ### Calculates the effects of soil moisture on decomposition rates according to the DEMETER model.
    ##references<< Foley, J. A. (1995), An equilibrium model of the terrestrial carbon budget, Tellus B, 47(3), 310-319.

    (M,         ##<< A scalar or vector containing values of soil moisture for which the effects on decomposition rates are calculated.
     Msat=100   ##<< A scalar representing saturated soil moisture. 
     )
   {
     0.25+0.75*(M/Msat)
      ### A scalar or a vector containing the effects of moisture on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Moisture=0:100
      plot(Moisture,fW.Demeter(Moisture),type="l",ylab="f(W) (unitless)", 
           main="Effects of soil moisture on decomposition rates according to the DEMETER model")
     }
)

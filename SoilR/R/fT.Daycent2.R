#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.Daycent2<-structure(
  function #Effects of temperature on decomposition rates according to the DAYCENT model
    ### Calculates the effects of temperature on decomposition rates according to the Daycent/Century models.
    ##references<<  Del Grosso, S. J., W. J. Parton, A. R. Mosier, E. A. Holland, E. Pendall, D. S. Schimel, and D. S. Ojima (2005), 
    ##Modeling soil CO2 emissions from ecosystems, Biogeochemistry, 73(1), 71-91.
    (Temp     ##<< A scalar or vector containing values of soil temperature for which the effects on decomposition rates are calculated.
     )
   {
      0.56+(1.46*atan(pi*0.0309*(Temp-15.7)))/pi
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Temperature=0:50
      plot(Temperature,fT.Daycent2(Temperature),type="l",
           ylab="f(T) (unitless)",
           main="Effects of temperature on decomposition rates according to the DAYCENT model")
    }
)

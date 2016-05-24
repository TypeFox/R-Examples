#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.KB<-structure(
  function #Effects of temperature on decomposition rates according to a model proposed by M. Kirschbaum (1995)
    ### Calculates the effects of temperature on decomposition rates according to a model proposed by Kirschbaum (1995).
    ##references<<  Kirschbaum, M. U. F. (1995), The temperature dependence of soil organic matter decomposition,
    ##and the effect of global warming on soil organic C storage, Soil Biology and Biochemistry, 27(6), 753-760.
    (Temp     ##<< a scalar or vector containing values of soil temperature for which the effects on decomposition rates are calculated
     )
   {
      exp(-3.764+0.204*Temp*(1-0.5*Temp/36.9))
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Temperature=0:50
      plot(Temperature,fT.KB(Temperature),type="l",ylab="f(T) (unitless)", 
           main="Effects of temperature on decomposition rates according to the Kirschbaum function")
    }
)

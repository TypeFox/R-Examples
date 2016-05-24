#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.RothC<-structure(
  function #Effects of temperature on decomposition rates according to the functions included in the RothC model
    ### Calculates the effects of temperature on decomposition rates according to the functions included in the RothC model.
    ##references<<  Jenkinson, D. S., S. P. S. Andrew, J. M. Lynch, M. J. Goss, and P. B. Tinker (1990), 
    ##The Turnover of Organic Carbon and Nitrogen in Soil, Philosophical Transactions: Biological Sciences, 329(1255), 361-368.
    (Temp     ##<< A scalar or vector containing values of temperature for which the effects on decomposition rates are calculated.
     )
   {
      47.9/(1+exp(106/(ifelse(Temp >= -18.3, Temp, NA)+18.3)))
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    ##note<< This function returns NA for Temp <= -18.3
    }
    ,
    ex=function(){
      Temperature=0:50
      plot(Temperature,fT.RothC(Temperature),type="l",ylab="f(T) (unitless)", 
           main="Effects of temperature on decomposition rates according to the RothC model")
    }
)

#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.Arrhenius<-structure(
  function #Effects of temperature on decomposition rates according the Arrhenius equation
    ### Calculates the effects of temperature on decomposition rates according to the Arrhenius equation.
    (Temp,     ##<< A scalar or vector containing values of temperature (in degrees Kelvin) for which the effects on decomposition rates are calculated.
     A=1000,  ##<< A scalar defining the pre-exponential factor.
     Ea=75000,   ##<< A scalar defining the activation energy in units of J mol^-1.
     Re=8.3144621   ##<< A scalar defining the universal gas contant in units of J K^-1 mol^-1.
     )
   {
      A*exp(-Ea/(Re*Temp))
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Temperature=273:300
      plot(Temperature,fT.Arrhenius(Temperature),type="l",ylab="f(T) (unitless)", xlab="Temperature (K)",
           main="Effects of temperature on decomposition rates according to the Arrhenius equation")
    }
)

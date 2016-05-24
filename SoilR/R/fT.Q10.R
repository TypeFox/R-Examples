#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.Q10<-structure(
  function #Effects of temperature on decomposition rates according to a Q10 function
    ### Calculates the effects of temperature on decomposition rates according to the modified Van't Hoff function (Q10 function).
    (Temp,     ##<< A scalar or vector containing values of temperature for which the effects on decomposition rates are calculated.
     k_ref=1,      ##<< A scalar representing the value of the decomposition rate at a reference temperature vaule.
     T_ref=10,  ##<< A scalar representing the reference temperature.
     Q10=2     ##<< A scalar. Temperature coefficient Q10.
     )
   {
     k_ref*Q10^((Temp-T_ref)/10)
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Temperature=0:50
      plot(Temperature,fT.Q10(Temperature),type="l",ylab="f(T) (unitless)", 
           main="Effects of temperature on decomposition rates according to a Q10 function")
      lines(Temperature, fT.Q10(Temperature,Q10=2.2),col=2)
      lines(Temperature, fT.Q10(Temperature,Q10=1.4),col=4)
      legend("topleft",c("Q10 = 2", "Q10 = 2.2", "Q10 = 1.4"),lty=c(1,1,1),col=c(1,2,4),bty="n")
    }
)

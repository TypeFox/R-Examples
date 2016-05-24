#
# vim:set ff=unix expandtab ts=2 sw=2:
fT.Century1<-structure(
  function #Effects of temperature on decomposition rates according the the CENTURY model
    ### Calculates the effects of temperature on decomposition rates according to the CENTURY model.
    ##references<<  Burke, I. C., J. P. Kaye, S. P. Bird, S. A. Hall, R. L. McCulley, and G. L. Sommerville. 2003. 
    ##Evaluating and testing models of terrestrial biogeochemistry: the role of temperature in controlling decomposition. 
    ##Pages 235-253 in C. D. Canham, J. J. Cole, and W. K. Lauenroth, editors. Models in ecosystem science. Princeton University Press, Princeton.
    (Temp,     ##<< A scalar or vector containing values of temperature for which the effects on decomposition rates are calculated.
     Tmax=45,   ##<< A scalar defining the maximum temperature in degrees C.
     Topt=35   ##<< A scalar defining the optimum temperature for the decomposition process in degrees C.
     )
   {
      (((Tmax-Temp)/(Tmax-Topt))^0.2)*exp((0.2/2.63)*(1-((Tmax-Temp)/(Tmax-Topt))^2.63))
      ### A scalar or a vector containing the effects of temperature on decomposition rates (unitless).
    }
    ,
    ex=function(){
      Temperature=0:50
      plot(Temperature,fT.Century1(Temperature),type="l",ylab="f(T) (unitless)",
           main="Effects of temperature on decomposition rates according to the Century model")
    }
)

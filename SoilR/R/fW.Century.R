#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.Century<- structure(
  function #Effects of moisture on decomposition rates according to the CENTURY model
    ### Calculates the effects of precipitation and potential evapotranspiration on decomposition rates.
    ##references<<  Adair, E. C., W. J. Parton, S. J. D. Grosso, W. L. Silver, M. E. Harmon, S. A. Hall, I. C. Burke, and S. C. Hart (2008), 
    ##Simple three-pool model accurately describes patterns of long-term litter decomposition in diverse climates, Global Change Biology, 14(11), 2636-2660.
    ## \code{\\} Parton, W. J., J. A. Morgan, R. H. Kelly, and D. S. Ojima (2001), Modeling soil C responses to environmental change in grassland systems, 
    ##in The potential of U.S. grazing lands to sequester carbon and mitigate the greenhouse effect, edited by R. F. Follett, J. M. Kimble and R. Lal, pp. 371-398, Lewis Publishers, Boca Raton.

     (PPT,     ##<< A scalar or vector containing values of monthly precipitation.
     PET   ##<< A scalar or vector containing values of potential evapotranspiration.
     )
   {
      1/(1+30*exp(-8.5*(PPT/PET)))
      ### A scalar or a vector containing the effects of precipitation and potential evapotranspiration on decomposition rates (unitless).
   }
    ,
    ex=function(){
      PPT=seq(0,1500,by=10)
      PET=rep(1500,length(PPT))
      PPT.PET=fW.Century(PPT,PET)
      plot(PPT/PET,PPT.PET, 
           ylab="f(PPT, PET) (unitless)", 
           main="Effects of precipitation and potential evapotranspiration on decomposition rates")
    }
)

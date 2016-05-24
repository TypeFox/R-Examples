#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.Daycent2<-structure(
  function # Effects of moisture on decomposition rates according to the DAYCENT model
    ### Calculates the effects of volumetric water content on decomposition rates according to the Daycent/Century models.
    ##references<<  Del Grosso, S. J., W. J. Parton, A. R. Mosier, E. A. Holland, E. Pendall, D. S. Schimel, and D. S. Ojima (2005), 
    ##Modeling soil CO2 emissions from ecosystems, Biogeochemistry, 73(1), 71-91.
    (W,      ##<< A scalar or vector of volumetric water content in percentage. 
     WP=0,   ##<< A scalar representing the wilting point in percentage.
     FC=100  ##<< A scalar representing the field capacity in percentage. 
     )
   {
     RWC=(W-WP)*100/(FC-WP)
     fRWC=5*(0.287+(atan(pi*0.009*(RWC-17.47)))/pi)
     return(data.frame(RWC,fRWC))
     ### A data frame with values of relative water content (RWC) and the effects of RWC on decomposition rates (fRWC).
    }
   ,
  ex=function(){
     W=10:90
     fW=fW.Daycent2(W,WP=10,FC=90)
     plot(fW,type="l",ylim=c(0,6)) #This plot reproduces Figure 1b, in del Grosso et al. (2005)
     }
)

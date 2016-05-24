#
# vim:set ff=unix expandtab ts=2 sw=2:
fW.Standcarb<-structure(
  function #Effects of moisture on decomposition rates according to the StandCarb model
    ### Calculates the effects of moisture on decomposition rates according to the StandCarb model.
    ##references<< Harmon, M. E., and J. B. Domingo (2001), A users guide to STANDCARB version 2.0: 
    ##A model to simulate carbon stores in forest stands. Oregon State University, Corvallis.

    (Moist,      ##<< A scalar or vector containing values of moisture content of a litter or soil pool (%).
     MatricShape=5,  ##<< A scalar that determines when matric limit is reduced to the point that decay can begin to occur. 
     MatricLag=0,    ##<< A scalar used to offset the curve to the left or right.
     MoistMin=30,  ##<< A scalar determining the minimum moisture content.
     MoistMax=350 ,      ##<< A scalar determining the maximum moisture content without diffusion limitations.
     DiffuseShape=15,   ##<< A scalar that determines the range of moisture contents where diffusion is not limiting.
     DiffuseLag=4  ##<< A scalar used to shift the point when moisture begins to limit diffusion.
     )
   {
      IncreaseRate=3/MoistMin  ##<< a parameter determining the point at which the matric limitation ends
      MatricLimit=(1-exp(-IncreaseRate*(Moist+MatricLag)))^MatricShape
      
      DiffuseLimit=exp(-1*(Moist/(MoistMax+DiffuseLag))^DiffuseShape)
      
      MoistDecayIndex=MatricLimit*DiffuseLimit
      
      return(data.frame(MatricLimit,DiffuseLimit,MoistDecayIndex))
     
      ### A data frame with limitation due to water potential (MatricLimit), limitation due to oxygen diffusion (DiffuseLimit), and the overall limitation of moisture on decomposition rates (MoistDecayIndex).
    }
    ,
    ex=function(){
      MC=0:500
      DeadFoliage=fW.Standcarb(MC)
      DeadBranch=fW.Standcarb(MC,MoistMax=200)
      DeadWood=fW.Standcarb(MC,MoistMax=150)
      StableSoil=fW.Standcarb(MC,MoistMin=15,MoistMax=100)
      plot(MC,DeadFoliage$MoistDecayIndex,type="l",xlab="Moisture Content (%)",
           ylab="f(W) (unitless)",
           main="Effects of moisture on decomposition rates according to the StandCarb model")
      lines(MC,DeadBranch$MoistDecayIndex,col=4)
      lines(MC,DeadWood$MoistDecayIndex,col=3)
      lines(MC,StableSoil$MoistDecayIndex,col=2)
      legend("topright",c("Dead Foliage","Dead Branch","Dead Wood","Stable Soil"),
             lty=c(1,1,1),col=c(1,4,3,2),bty="n")
    }
)

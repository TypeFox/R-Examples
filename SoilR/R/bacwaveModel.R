#
# vim:set ff=unix expandtab ts=2 sw=2:
bacwaveModel<-structure(
  function # Implementation of the microbial model Bacwave (bacterial waves)
  ### This function implements the microbial model Bacwave (bacterial waves), a two-pool model with a bacterial and a substrate pool. It is a special case of the general nonlinear model.
  ##details<< This implementation containts default parameters presented in Zelenev et al. (2000). It produces nonlinear damped oscillations in the form of a stable focus.
  ##references<< Zelenev, V.V., A.H.C. van Bruggen, A.M. Semenov. 2000. ``BACWAVE,'' a spatial-temporal model for traveling waves of bacterial populations
  ## in response to a moving carbon source in soil. Microbail Ecology 40: 260-272.
  (t, ##<< vector of times (in hours) to calculate a solution.
   umax=0.063, ##<< a scalar representing the maximal relative growth rate of bacteria (hr-1)
   ks=3.0, ##<< a scalar representing the substrate constant for growth (ug C /ml soil solution)
   theta=0.23, ##<< a scalar representing soil water content (ml solution/cm3 soil)
   Dmax=0.26, ##<< a scalar representing the maximal relative death rate of bacteria (hr-1)
   kd=14.5, ##<< a scalar representing the substrate constant for death of bacteria (ug C/ml soil solution)
   kr=0.4, ##<< a scalar representing the fraction of death biomass recycling to substrate (unitless)
   Y=0.44, ##<< a scalar representing the yield coefficient for bacteria (ug C/ugC)
   ival=c(S0=0.5,X0=1.5), ##<< a vector of length 2 with the initial values for the substrate and the bacterial pools (ug C/cm3)
   BGF=0.15, ##<< a scalar representing the constant background flux of substrate (ug C/cm3 soil/hr)
   ExuM=8, ##<< a scalar representing the maximal exudation rate (ug C/(hr cm3 soil))
   ExuT=0.8 ##<< a scalar representing the  time constant for exudation, responsible for duration of exudation (1/hr).
    )
    {
    t_start=min(t)
    t_end=max(t)
    nr=2
    
    #function with system of equations
    f=function(C,t){
      S=C[[1]] #Substrate pool
      X=C[[2]] #Bacterial pool
      O=matrix(byrow=TRUE,nrow=2,c((X/Y)*(umax*S)/((ks*theta)+S),
                                   (Dmax*kd*X)/(kd+(S/theta))))
      return(O)
    }
    
    #List of transfer coefficients for T matrix
    alpha=list()
    alpha[["1_to_2"]]=function(C,t){
      Y
    }
    alpha[["2_to_1"]]=function(C,t){
      kr
    }
    
    Anl=new("TransportDecompositionOperator",t_start,Inf,nr,alpha,f)
    
    inputrates=BoundInFlux(
      function(t){
        matrix(
          nrow=nr,
          ncol=1,
          c(BGF+(ExuM*exp(-ExuT*t)),  0)
        )
      },
      t_start,
      t_end
    )
    
    modnl=GeneralNlModel(t, Anl, ival, inputrates, deSolve.lsoda.wrapper)
    
    return(modnl)
    ### An object of class NlModel that can be further queried.
  }
  ,
  ex=function(){
  
    hours=seq(0,800,0.1)
    
    #Run the model with default parameter values
    bcmodel=bacwaveModel(t=hours)
    Cpools=getC(bcmodel)
    
    #Time solution
    matplot(hours,Cpools,type="l",ylab="Concentrations",xlab="Hours",lty=1,ylim=c(0,max(Cpools)*1.2))
    legend("topleft",c("Substrate", "Microbial biomass"),lty=1,col=c(1,2),bty="n")
    
    #State-space diagram
    plot(Cpools[,2],Cpools[,1],type="l",ylab="Substrate",xlab="Microbial biomass")
    
    #Microbial biomass over time
    plot(hours,Cpools[,2],type="l",col=2,xlab="Hours",ylab="Microbial biomass")
    
  }
)

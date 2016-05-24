#
# vim:set ff=unix expandtab ts=2 sw=2:
TwopMMmodel<-structure(
  function # Implementation of a two-pool Michaelis-Menten model
  ### This function implements a two-pool Michaelis-Meneten model with a microbial biomass and a substrate pool.
  ##details<< This implementation is similar to the model described in Manzoni and Porporato (2007). 
  ##references<< Manzoni, S, A. Porporato (2007). A theoretical analysis of nonlinearities and feedbacks in soil carbon and nitrogen cycles.
  ## Soil Biology and Biochemistry 39: 1542-1556.
  (t, ##<< vector of times (in days) to calculate a solution.
   ks=0.000018, ##<< a scalar representing SOM decomposition rate (m3 d-1 (gCB)-1)
   kb=0.007, ##<< a scalar representing microbial decay rate (d-1)
   Km=900, ##<< a scalar representing the Michaelis constant (g m-3)
   r=0.6, ##<< a scalar representing the respired carbon fraction (unitless)
   Af=1, ##<< a scalar representing the Activity factor; i.e. a temperature and moisture modifier (unitless)
   ADD=3.2, ##<< a scalar representing the annual C input to the soil (g m-3 d-1)
   ival # a vector of length 2 with the initial values of the SOM pool and the microbial biomass pool (g m-3)
  )
{
    t_start=min(t)
    t_end=max(t)
    nr=2
    if(length(ival)!=2) stop("The vector of initial values ival must be of length 2")
    
    
    #function with system of equations
    f=function(C,t){
      S=C[1] #SOM pool
      B=C[2] #Microbial biomass pool
      O=matrix(byrow=TRUE,nrow=2,ncol=1,c((Af*ks)*B*(S/(Km+S)),
                                   kb*B))
      return(O)
    }
    
    #List of transfer coefficients for T matrix
    alpha=list()
    alpha[["1_to_2"]]=function(C,t){
      1-r
    }
    alpha[["2_to_1"]]=function(C,t){
      1
    }
    
    Anl=new("TransportDecompositionOperator",t_start,Inf,nr,alpha,f)
    
    inputrates=BoundInFlux(
      function(t){
        matrix(
          nrow=nr,
          ncol=1,
          c(ADD,  0)
        )
      },
      t_start,
      t_end
    )

    modnl=GeneralNlModel( t, Anl, ival, inputrates, deSolve.lsoda.wrapper)
    
    return(modnl)
    ### An object of class NlModel that can be further queried.
  }
  ,
  ex=function(){
    
    days=seq(0,1000,0.5)
    
    #Run the model with default parameter values
    MMmodel=TwopMMmodel(t=days,ival=c(100,10))
    Cpools=getC(MMmodel)
    
    #Time solution
    matplot(days,Cpools,type="l",ylab="Concentrations",xlab="Days",lty=1,ylim=c(0,max(Cpools)*1.2))
    legend("topleft",c("SOM-C", "Microbial biomass"),lty=1,col=c(1,2),bty="n")
    
    ks=0.000018
    kb=0.007
    r=0.6
    ADD=3.2
    
    #Analytical solution of fixed points
    Cs=kb/((1-r)*ks)
    abline(h=Cs,lty=2)
    
    Cb=(ADD*(1-r))/(r*kb)
    abline(h=Cb,lty=2,col=2)
    
    #State-space diagram
    plot(Cpools[,2],Cpools[,1],type="l",ylab="SOM-C",xlab="Microbial biomass")
    
    #Microbial biomass over time
    plot(days,Cpools[,2],type="l",col=2,xlab="Days",ylab="Microbial biomass")
    
    #The default parameterization exhaust the microbial biomass.
    #A different behavior is obtained by increasing ks and decreasing kb
    MMmodel=TwopMMmodel(t=days,ival=c(972,304) ,Af=3,kb=0.0000001)
    Cpools=getC(MMmodel)
    
    matplot(days,Cpools,type="l",ylab="Concentrations",xlab="Days",lty=1,ylim=c(0,max(Cpools)*1.2))
    legend("topleft",c("SOM-C", "Microbial biomass"),lty=1,col=c(1,2),bty="n")
    
    plot(Cpools[,2],Cpools[,1],type="l",ylab="SOM-C",xlab="Microbial biomass")
    
    plot(days,Cpools[,2],type="l",col=2,xlab="Days",ylab="Microbial biomass")
    
  }
)

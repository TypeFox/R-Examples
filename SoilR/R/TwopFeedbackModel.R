#
# vim:set ff=unix expandtab ts=2 sw=2:
TwopFeedbackModel<-structure(
    function #Implementation of a two pool model with feedback structure
    ### This function creates a model for two pools connected with feedback. It is a wrapper for the more general function \code{\link{GeneralModel}}.
    ##references<< Sierra, C.A., M. Mueller, S.E. Trumbore. 2012. Models of soil organic matter decomposition: the SoilR package version 1.0. Geoscientific Model Development 5, 1045-1060.
     (t,    	##<< A vector containing the points in time where the solution is sought.
      ks,	##<< A vector of length 2 with the values of the decomposition rate for pools 1 and 2. 
      a21, ##<< A scalar with the value of the transfer rate from pool 1 to pool 2.
      a12, ##<< A scalar with the value of the transfer rate from pool 2 to pool 1.
      C0,	##<< A vector of length 2 containing the initial amount of carbon for the 2 pools.
      In,     ##<< A data.frame object specifying the amount of litter inputs by time. 
      xi=1,  ##<< A scalar or data.frame object specifying the external (environmental and/or edaphic) effects on decomposition rates.
      solver=deSolve.lsoda.wrapper,  ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.
      pass=FALSE  ##<< if TRUE forces the constructor to create the model even if it is invalid 
    )	
    { 
      t_start=min(t)
      t_end=max(t)
      if(length(ks)!=2) stop("ks must be of length = 2")
      if(length(C0)!=2) stop("the vector with initial conditions must be of length = 2")
      
      if(length(In)==1){
          inputFluxes=BoundInFlux(
            function(t){matrix(nrow=2,ncol=1,c(In,0))},
            t_start,
            t_end
          )
      }
      if(class(In)=="data.frame"){
         x=In[,1]  
         y=In[,2]  
         inputFlux=splinefun(x,y)
         inputFluxes=BoundInFlux(
            function(t){matrix(nrow=2,ncol=1,c(inputFlux(t),0))},
            min(x),
            max(x)
         )
        }
      A=-1*abs(diag(ks))
      A[2,1]=a21
      A[1,2]=a12
      
      if(length(xi)==1) fX=function(t){xi}
      if(class(xi)=="data.frame"){
        X=xi[,1]
        Y=xi[,2]
        fX=function(t){as.numeric(spline(X,Y,xout=t)[2])}
       }
      Af=BoundLinDecompOp(
        function(t){fX(t)*A},
        t_start,
        t_end
      )
      Mod=GeneralModel(t=t,A=Af,ivList=C0,inputFluxes=inputFluxes,solver,pass)
     return(Mod)
### A Model Object that can be further queried 
      ##seealso<< \code{\link{TwopParallelModel}}, \code{\link{TwopSeriesModel}} 
    }
    ,
    ex=function(){
    
    #This example show the difference between the three types of two-pool models  
    times=seq(0,20,by=0.1)
    ks=c(k1=0.8,k2=0.00605)
    C0=c(C10=5,C20=5)
      
    Temp=rnorm(times,15,2)
    WC=runif(times,10,20)
    TempEffect=data.frame(times,fT=fT.Daycent1(Temp))
    MoistEffect=data.frame(times, fW=fW.Daycent2(WC)[2])
    
    Inmean=1
    InRand=data.frame(times,Random.inputs=rnorm(length(times),Inmean,0.2))
    InSin=data.frame(times,Inmean+0.5*sin(times*pi*2))
          
    Parallel=TwopParallelModel(t=times,ks=ks,C0=C0,In=Inmean,gam=0.9,
                               xi=(fT.Daycent1(15)*fW.Demeter(15)))
    Series=TwopSeriesModel(t=times,ks=ks,a21=0.2*ks[1],C0=C0,In=InSin,
                           xi=(fT.Daycent1(15)*fW.Demeter(15)))
    Feedback=TwopFeedbackModel(t=times,ks=ks,a21=0.2*ks[1],a12=0.5*ks[2],C0=C0,
                               In=InRand,xi=MoistEffect)
    
    CtP=getC(Parallel)
    CtS=getC(Series)
    CtF=getC(Feedback)
    
    RtP=getReleaseFlux(Parallel)
    RtS=getReleaseFlux(Series)
    RtF=getReleaseFlux(Feedback)
    
    par(mfrow=c(2,1),mar=c(4,4,1,1))
    plot(times,rowSums(CtP),type="l",ylim=c(0,20),ylab="Carbon stocks (arbitrary units)",xlab=" ")
    lines(times,rowSums(CtS),col=2)
    lines(times,rowSums(CtF),col=3)
    legend("topleft",c("Two-pool Parallel","Two-pool Series","Two-pool Feedback"),
           lty=c(1,1,1),col=c(1,2,3),bty="n")
          
    plot(times,rowSums(RtP),type="l",ylim=c(0,3),ylab="Carbon release (arbitrary units)", xlab="Time")
    lines(times,rowSums(RtS),col=2)
    lines(times,rowSums(RtF),col=3)
    par(mfrow=c(1,1))

}
)

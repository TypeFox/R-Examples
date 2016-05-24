#
# vim:set ff=unix expandtab ts=2 sw=2:
TwopSeriesModel<-structure(
    function #Implementation of a two pool model with series structure
    ### This function creates a model for two pools connected in series. It is a wrapper for the more general function \code{\link{GeneralModel}}.
    ##references<< Sierra, C.A., M. Mueller, S.E. Trumbore. 2012. Models of soil organic matter decomposition: the SoilR package version 1.0. Geoscientific Model Development 5, 1045-1060.
     (t,  		##<< A vector containing the points in time where the solution is sought.
      ks,	##<< A vector of length 2 with the values of the decomposition rate for pools 1 and 2. 
      a21, ##<< A scalar with the value of the transfer rate from pool 1 to pool 2.
      C0,	##<< A vector of length 2 containing the initial amount of carbon for the 2 pools.
      In,     ##<< A scalar or a data.frame object specifying the amount of litter inputs by time. 
      xi=1,   ##<< A scalar or a data.frame specifying the external (environmental and/or edaphic) effects on decomposition rates. 
      solver=deSolve.lsoda.wrapper,  ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.
 pass=FALSE  ##<< if TRUE Forces the constructor to create the model even if it is invalid 
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
      
      if(length(xi)==1) fX=function(t){xi}
      if(class(xi)=="data.frame"){
      X=xi[,1]
      Y=xi[,2]
      fX=splinefun(X,Y)
      }
      Af=BoundLinDecompOp(
        function(t){fX(t)*A},
        t_start,
        t_end
      )
      Mod=GeneralModel(t=t,A=Af,ivList=C0,inputFluxes=inputFluxes,pass=pass)
     return(Mod)
### A Model Object that can be further queried 
      ##seealso<< \code{\link{TwopParallelModel}},\code{\link{TwopFeedbackModel}} 
    }
    ,
    ex=function(){
      t_start=0 
      t_end=10 
      tn=50
      timestep=(t_end-t_start)/tn 
      t=seq(t_start,t_end,timestep) 
      ks=c(k1=0.8,k2=0.4)
      a21=0.5
      C0=c(C10=100,C20=150)
      In = 30
      
      Temp=rnorm(t,15,1)
      TempEffect=data.frame(t,fT.Daycent1(Temp))
      
      Ex1=TwopSeriesModel(t,ks,a21,C0,In,xi=TempEffect)
      Ct=getC(Ex1)
      Rt=getReleaseFlux(Ex1)
      
      plot(t,rowSums(Ct),type="l",ylab="Carbon stocks (arbitrary units)",
           xlab="Time (arbitrary units)",lwd=2,ylim=c(0,sum(Ct[1,]))) 
      lines(t,Ct[,1],col=2)
      lines(t,Ct[,2],col=4) 
      legend("bottomright",c("Total C","C in pool 1", "C in pool 2"),
             lty=c(1,1,1),col=c(1,2,4),lwd=c(2,1,1),bty="n")

      plot(t,rowSums(Rt),type="l",ylab="Carbon released (arbitrary units)",
           xlab="Time (arbitrary units)",lwd=2,ylim=c(0,sum(Rt[1,]))) 
      lines(t,Rt[,1],col=2)
      lines(t,Rt[,2],col=4) 
      legend("topright",c("Total C release","C release from pool 1", "C release from pool 2"),
             lty=c(1,1,1),col=c(1,2,4),lwd=c(2,1,1),bty="n")
      
      Inr=data.frame(t,Random.inputs=rnorm(length(t),30,5))
      plot(Inr)
      
      Ex2=TwopSeriesModel(t,ks,a21,C0,In=Inr,xi=fT.Q10(15))
      Ctr=getC(Ex2)
      Rtr=getReleaseFlux(Ex2)
      
      plot(t,rowSums(Ctr),type="l",ylab="Carbon stocks (arbitrary units)",
           xlab="Time (arbitrary units)",lwd=2,ylim=c(0,sum(Ctr[1,]))) 
      lines(t,Ctr[,1],col=2)
      lines(t,Ctr[,2],col=4) 
      legend("topright",c("Total C","C in pool 1", "C in pool 2"),
             lty=c(1,1,1),col=c(1,2,4),lwd=c(2,1,1),bty="n")

      plot(t,rowSums(Rtr),type="l",ylab="Carbon released (arbitrary units)",
           xlab="Time (arbitrary units)",lwd=2,ylim=c(0,sum(Rtr[1,]))) 
      lines(t,Rtr[,1],col=2)
      lines(t,Rtr[,2],col=4) 
      legend("topright",c("Total C release","C release from pool 1", "C release from pool 2"),
             lty=c(1,1,1),col=c(1,2,4),lwd=c(2,1,1),bty="n")
}
)

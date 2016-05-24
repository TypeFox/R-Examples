#
# vim:set ff=unix expandtab ts=2 sw=2:
TwopParallelModel<-structure(
    function #Implementation of a two pool model with parallel structure
    ### This function creates a model for two independent (parallel) pools. 
    ### It is a wrapper for the more general function \code{\link{ParallelModel}} that can handle an arbitrary number of pools.
    ##references<< Sierra, C.A., M. Mueller, S.E. Trumbore. 2012. Models of soil organic matter decomposition: the SoilR package version 1.0. Geoscientific Model Development 5, 1045-1060.
     (t,			##<< A vector containing the points in time where the solution is sought.
     ks,	##<< A vector of length 2 containing the decomposition rates for the 2 pools. 
     C0,	##<< A vector of length 2 containing the initial amount of carbon for the 2 pools.
     In,     ##<< A scalar or a data.frame object specifying the amount of litter inputs by time.
     gam,  ##<< A scalar representing the partitioning coefficient, i.e. the proportion from the total amount of inputs that goes to pool 1.
     xi=1,   ##<< A scalar or a data.frame specifying the external (environmental and/or edaphic) effects on decomposition rates. 
     solver=deSolve.lsoda.wrapper, ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.
     pass=FALSE  ##<< Forces the constructor to create the model even if it is invalid 
    )	
    { 
      t_start=min(t)
      t_stop=max(t)
      if(length(ks)!=2) stop("ks must be of length = 2")
      if(length(C0)!=2) stop("the vector with initial conditions must be of length = 2")
      if(gam > 1 | gam < 0) stop("The the partitioning coefficient gam is outside the interval [0,1]")
      
      if(length(In)==1) inputrates=BoundInFlux(
        function(t){matrix(nrow=2,ncol=1,c(gam*In,(1-gam)*In))},
        t_start,
        t_stop
      )
      if(class(In)=="data.frame"){
         x=In[,1]  
         y=In[,2]  
         inputrate=function(t0){as.numeric(spline(x,y,xout=t0)[2])}
         inputrates=BoundInFlux(
            function(t){
                matrix(nrow=2,ncol=1,
                    c(
                        gam*inputrate(t),
                        (1-gam)*inputrate(t)
                    )
                )
            },
            min(x),
            max(x)
         )   
        }
      
      if(length(xi)==1) fX=function(t){xi}
      if(class(xi)=="data.frame"){
      X=xi[,1]
      Y=xi[,2]
      fX=function(t){as.numeric(spline(X,Y,xout=t)[2])}
      }

      coeffs=TimeMap.new(
        t_start,
        t_stop,
        function(t){fX(t)*(-1*abs(ks))}
      )
      res=ParallelModel(t,coeffs,startvalues=C0,inputrates,solver,pass=pass)
      ### A Model Object that can be further queried 
      ##seealso<< \code{\link{ThreepParallelModel}} 
    }
    ,
    ex=function(){
      t_start=0 
      t_end=10 
      tn=50
      timestep=(t_end-t_start)/tn 
      t=seq(t_start,t_end,timestep) 

      Ex=TwopParallelModel(t,ks=c(k1=0.5,k2=0.2),C0=c(c10=100, c20=150),In=10,gam=0.7,xi=0.5)
      Ct=getC(Ex)

      plot(t,rowSums(Ct),type="l",lwd=2,
           ylab="Carbon stocks (arbitrary units)",xlab="Time",ylim=c(0,sum(Ct[1,]))) 
      lines(t,Ct[,1],col=2)
      lines(t,Ct[,2],col=4)
      legend("topright",c("Total C","C in pool 1", "C in pool 2"),
             lty=c(1,1,1),col=c(1,2,4),lwd=c(2,1,1),bty="n")
      
      Rt=getReleaseFlux(Ex)
      plot(t,rowSums(Rt),type="l",ylab="Carbon released (arbitrary units)",
           xlab="Time",lwd=2,ylim=c(0,sum(Rt[1,]))) 
      lines(t,Rt[,1],col=2)
      lines(t,Rt[,2],col=4) 
      legend("topleft",c("Total C release","C release from pool 1", "C release from pool 2"),
             lty=c(1,1,1),col=c(1,2,4),lwd=c(2,1,1),bty="n")

}
)

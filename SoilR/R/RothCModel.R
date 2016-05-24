#
# vim:set ff=unix expandtab ts=2 sw=2:
RothCModel<-structure(
    function #Implementation of the RothCModel
    ### This function implements the RothC model of Jenkinson et al. It is a wrapper for the more general function \code{\link{GeneralModel}}.
    ##references<< Jenkinson, D. S., S. P. S. Andrew, J. M. Lynch, M. J. Goss, and P. B. Tinker. 1990. The Turnover of Organic Carbon and Nitrogen in Soil. 
    ##Philosophical Transactions: Biological Sciences 329:361-368.
    ##Sierra, C.A., M. Mueller, S.E. Trumbore. 2012. Models of soil organic matter decomposition: the SoilR package version 1.0. Geoscientific Model Development 5, 1045-1060.
    (t,      ##<< A vector containing the points in time where the solution is sought.
      ks=c(k.DPM=10,k.RPM=0.3,k.BIO=0.66,k.HUM=0.02,k.IOM=0),	##<< A vector of lenght 5 containing the values of the decomposition rates for the different pools
      C0=c(0,0,0,0,2.7),	##<< A vector of length 5 containing the initial amount of carbon for the 5 pools.
      In=1.7,    ##<< A scalar or data.frame object specifying the amount of litter inputs by time. 
      DR=1.44, ##<< A scalar representing the ratio of decomposable plant material to resistant plant material (DPM/RPM).
      clay=23.4, ##<< Percent clay in mineral soil. 
      xi=1,  ##<< A scalar or data.frame object specifying the external (environmental and/or edaphic) effects on decomposition rates.
      solver=deSolve.lsoda.wrapper,  ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.
      pass=FALSE  ##<< if TRUE forces the constructor to create the model even if it is invalid 
    )	
    { 
      t_start=min(t)
      t_end=max(t)
      if(length(ks)!=5) stop("ks must be of length = 5")
      if(length(C0)!=5) stop("the vector with initial conditions must be of length = 5")
      
      if(length(In)==1){
          inputFluxes=BoundInFlux(
            function(t){matrix(nrow=5,ncol=1,c(In*(DR/(DR+1)),In*(1/(DR+1)),0,0,0))},
            t_start,
            t_end
        )
      }
      if(class(In)=="data.frame"){
         x=In[,1]  
         y=In[,2]  
         inputFlux=splinefun(x,y)
          inputFluxes=BoundInFlux(
            function(t){matrix(nrow=5,ncol=1,c(inputFlux(t)*(DR/(DR+1)),inputFlux(t)*(1/(DR+1)),0,0,0))},
            min(x),
            max(x)
          )
        }

      x=1.67*(1.85+1.60*exp(-0.0786*clay))
      B=0.46/(x+1) # Proportion that goes to the BIO pool
      H=0.54/(x+1) # Proportion that goes to the HUM pool

      ai3=B*ks
      ai4=H*ks

      A=diag(-ks)
      A[3,]=A[3,]+ai3
      A[4,]=A[4,]+ai4

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
      Mod=GeneralModel(t=t,A=Af,ivList=C0,inputFluxes=inputFluxes,solverfunc=solver,pass=pass)
     return(Mod)
### A Model Object that can be further queried 
      ##seealso<< \code{\link{ICBMModel}} 
    }
    ,
    ex=function(){
      t=0:500 
      Ex=RothCModel(t)
      Ct=getC(Ex)
      Rt=getReleaseFlux(Ex)
      
     plot(
       t,
       Ct[,1],
       type="l",
       col=1,
        ylim=c(0,25),
       ylab=expression(paste("Carbon stores (Mg C", ha^-1,")")),
       xlab="Time (years)",
       lwd=2
     ) 
     lines(t,Ct[,2],col=2,lwd=2,lty=2) 
     lines(t,Ct[,3],col=3,lwd=2,lty=3)
     lines(t,Ct[,4],col=4,lwd=2,lty=4)
     lines(t,Ct[,5],col=5,lwd=2,lty=5)
     lines(t,rowSums(Ct),lwd=2)
     legend(
        "topright",
        c(
          "Pool 1, DPM",
          "Pool 2, RPM",
          "Pool 3, BIO",
          "Pool 4, HUM",
          "Pool 5, IOM",
          "Total Carbon"
        ),
        lty=c(1:5,1),
        lwd=rep(2,5),
        col=c(1,2,3,4,5,"black")
        ,bty="n"
    )

     plot(t,Rt[,1],type="l",ylim=c(0,2),ylab="Respiration (Mg C ha-1 yr-1)",xlab="Time") 
     lines(t,Rt[,2],col=2) 
     lines(t,Rt[,3],col=3) 
     lines(t,Rt[,4],col=4) 
     lines(t,Rt[,5],col=5) 
     lines(t,rowSums(Rt),lwd=2) 
     legend(
        "topright",
        c("Pool 1, DPM", 
          "Pool 2, RPM",
          "Pool 3, BIO",
          "Pool 4, HUM",
          "Pool 5, IOM",
          "Total Respiration"
        ), 
        lty=c(1,1,1,1,1,1),
        lwd=c(1, 1,1,1,1,2),
        col=c(1,2,3,4,5,1),
        bty="n"
      )

}
)

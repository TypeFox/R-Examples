#
# vim:set ff=unix expandtab ts=2 sw=2:
Yasso07Model<-structure(
    function #Implementation of the Yasso07 model
    ### This function creates a model for five pools as described in Tuomi et al. (2008)
    ##references<< Tuomi, M., Thum, T., Jarvinen, H., Fronzek, S., Berg, B., Harmon, M., Trofymow, J., Sevanto, S., and Liski, J. (2009).
    ##Leaf litter decomposition-estimates of global variability based on Yasso07 model. Ecological Modelling, 220:3362 - 3371.
     (t,      ##<< A vector containing the points in time where the solution is sought.
      ks=c(kA=0.66, kW=4.3, kE=0.35, kN=0.22, kH=0.0033),	##<< A vector of lenght 5 containing the values of the decomposition rates for each pool.
      p=c(p1=0.32,p2=0.01,p3=0.93,p4=0.34,p5=0,p6=0,p7=0,p8=0,p9=0.01,p10=0,p11=0,p12=0.92,pH=0.04), ##<< A vector of length 13 containing transfer coefficients among different pools.
      C0,	##<< A vector containing the initial amount of carbon for the 5 pools. The length of this vector must be 5.
      In,     ##<< A single scalar or data.frame object specifying the amount of litter inputs by time
      xi=1,  ##<< A scalar or data.frame object specifying the external (environmental and/or edaphic) effects on decomposition rates.
      solver=deSolve.lsoda.wrapper,  ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.
      pass=FALSE  ##<< if TRUE forces the constructor to create the model even if it is invalid 
    )	
    { 
      t_start=min(t)
      t_end=max(t)
      if(length(ks)!=5) stop("ks must be of length = 5")
      if(length(C0)!=5) stop("the vector with initial conditions must be of length = 5")
      if(length(p)!=13) stop("The vector of transfer coefficients p must be of length = 13")
      
      if(length(In)==1){
          inputFluxes=BoundInFlux(
            function(t){matrix(nrow=5,ncol=1,c(In,0,0,0,0))},
            t_start,
            t_end
        )
      }
      if(class(In)=="data.frame"){
       inputFluxes=BoundInFlux(In)
      }
      A1=abs(diag(ks))
      Ap=diag(-1,5,5)
      Ap[1,2]=p[1]
      Ap[1,3]=p[2]
      Ap[1,4]=p[3]
      Ap[2,1]=p[4]
      Ap[2,3]=p[5]
      Ap[2,4]=p[6]
      Ap[3,1]=p[7]
      Ap[3,2]=p[8]
      Ap[3,4]=p[9]
      Ap[4,1]=p[10]
      Ap[4,2]=p[11]
      Ap[4,3]=p[12]
      Ap[5,1:4]=p[13]
      
      A=Ap%*%A1
      
      if(length(xi)==1){
	fX=function(t){xi}
	Af=BoundLinDecompOp(function(t) fX(t)*A,t_start,t_end)
	}
      if(class(xi)=="data.frame"){
	      X=xi[,1]
      	Y=xi[,2]
        fX=splinefun(X,Y)
	Af=BoundLinDecompOp(function(t){fX(t)*A},min(X),max(X))
       }
      Mod=GeneralModel(t=t,A=Af,ivList=C0,inputFluxes=inputFluxes,solver,pass)
     return(Mod)
### A Model Object that can be further queried 
      ##seealso<< \code{\link{ICBMModel}}, \code{\link{RothCModel}}
    }
    ,
    ex=function(){
      years=seq(0,50,0.1) 
      C0=rep(100,5)
      In=0
      
      Ex1=Yasso07Model(t=years,C0=C0,In=In)
      Ct=getC(Ex1)
      Rt=getReleaseFlux(Ex1)
      
      plotCPool(years,Ct,col=1:5,xlab="years",ylab="C pool",
                ylim=c(0,max(Ct)))
      legend("topright",c("xA","xW","xE","xN","xH"),lty=1,col=1:5,bty="n")
      
      plotCPool(years,Rt,col=1:5,xlab="years",ylab="Respiration",ylim=c(0,50))
      legend("topright",c("xA","xW","xE","xN","xH"),lty=1,col=1:5,bty="n")

}
)

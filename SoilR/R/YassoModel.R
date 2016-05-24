#
# vim:set ff=unix expandtab ts=2 sw=2:
YassoModel<-structure(
  function #Implementation of the Yasso model.
  ### This function creates a model for seven pools as described in Liski et al. (2005). Model not yet implemented due to lack of data in original publication: values of vector p not completely described in paper. 0.1 was assumed. 
  ##references<< Liski, J., Palosuo, T., Peltoniemi, M., and Sievanen, R. (2005).
  ##Carbon and decomposition model Yasso for forest soils. Ecological Modelling, 189:168-182.
  (t,      ##<< A vector containing the points in time where the solution is sought.
   ks=c(a_fwl=0.54, a_cwl=0.03, k_ext=0.48, k_cel=0.3, k_lig=0.22, k_hum1=0.012,k_hum2=0.0012),  ##<< A vector of lenght 7 containing the values of the exposure and decomposition rates for each pool.
   p=c(fwl_ext=0.1, cwl_ext=0.1, fwl_cel=0.1, cwl_cel=0.1, fwl_lig=0.1, cwl_lig=0.1, pext=0.05, pcel=0.24, plig=0.77, phum1=0.51), ##<< A vector of containing transfer coefficients among different pools.
   C0,	##<< A vector containing the initial amount of carbon for the 7 pools. The length of this vector must be 7.
   In=c(u_fwl=0.0758, u_cwl=0.0866, u_nwl_cnwl_ext=0.251*0.3,u_nwl_cnwl_cel=0.251*0.3,u_nwl_cnwl_lig=0.251*0.3,0,0),     ##<< A vector of constatn litter inputs. 
   xi=1,  ##<< A scalar or data.frame object specifying the external (environmental and/or edaphic) effects on decomposition rates.
   solver=deSolve.lsoda.wrapper,  ##<< A function that solves the system of ODEs. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.
   pass=FALSE  ##<< if TRUE forces the constructor to create the model even if it is invalid 
  )	
  { 
    t_start=min(t)
    t_end=max(t)
    if(length(ks)!=7) stop("ks must be of length = 7")
    if(length(C0)!=7) stop("the vector with initial conditions must be of length = 7")
    if(length(p)!=10) stop("The vector of transfer coefficients p must be of length = 10")
    
    if(length(In)==7){
      inputFluxes=BoundInFlux(
        function(t){matrix(nrow=7,ncol=1,In)},
        t_start,
        t_end
      )
    }
    if(class(In)=="data.frame") stop("Inputs must be a vector of length 7")
      
    A1=abs(diag(ks))
    Ap=diag(-1,7,7)
    Ap[3,1]=p[1]
    Ap[3,2]=p[2]
    Ap[4,1]=p[3]
    Ap[4,2]=p[4]
    Ap[5,1]=p[5]
    Ap[5,2]=p[6]
    Ap[5,3]=p[7]
    Ap[5,4]=p[8]
    Ap[6,5]=p[9]
    Ap[7,6]=p[10]
    
    A=Ap%*%A1
    
    if(length(xi)==1){
      fX=function(t){xi}
      Af=BoundLinDecompOp(function(t){fX(t)*A},t_start,t_end)
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
    ##seealso<< \code{\link{ThreepParallelModel}}, \code{\link{ThreepSeriesModel}}
  }
  ,
  ex=function(){
    years=seq(0,500,0.5) 
    C0=rep(100,7)
    
    Ex1=YassoModel(t=years,C0=C0)
    Ct=getC(Ex1)
    Rt=getReleaseFlux(Ex1)
    
    plotCPool(years,Ct,col=1:7,xlab="years",ylab="C pool",ylim=c(0,200))
    legend("topright",c("fwl","cwl","ext","cel","lig","hum1","hum2"),lty=1,col=1:7,bty="n")
    
    plotCPool(years,Rt,col=1:7,xlab="years",ylab="Respiration",ylim=c(0,50))
    legend("topright",c("fwl","cwl","ext","cel","lig","hum1","hum2"),lty=1,col=1:7,bty="n")
    
  }
)

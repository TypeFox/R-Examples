#
# vim:set ff=unix expandtab ts=2 sw=2:
GeneralNlModel=structure(function #Use this function to create objects of class NlModel.
### The function creates a numerical model  for n arbitrarily connected pools.
### It is one of the constructors of class NlModel.                       
### It is used by some more specialized wrapper functions, but can also be used directly. 

(t,			##<< A vector containing the points in time where the solution is sought.
 TO,			##<< A object describing the model decay rates for the n pools, connection and feedback coefficients. The number of pools n must be consistent with the number of initial values and input fluxes. 
 ivList,		##<< A numeric vector containing the initial amount of carbon for the n pools. The length of this vector is equal to the number of pools. 
 inputFluxes, ##<< A TimeMap object consisting of a vector valued function describing the inputs to the pools as funtions of time \code{\link{TimeMap.new}}.
 solverfunc=deSolve.lsoda.wrapper,		##<< The function used by to actually solve the ODE system.  
 pass=FALSE  ##<< Forces the constructor to create the model even if it is invalid. If set to TRUE, does not enforce the requirements for a biologically meaningful model, e.g. does not check if negative values of respiration are calculated.
 )
{
   # first we test that the input values are compatible with each other
   obj=new(Class="NlModel",t,TO,ivList,inputFluxes,solverfunc,pass)
   return(obj)
   ### The function returns an object of class NlModel that can be further queried by other functions such as \code{\link{getC}}.
   ##seealso<< \code{\link{GeneralModel}}.
}
,ex=function(){
 t_start=0
  t_end=20
  tn=100
  timestep=(t_end-t_start)/tn
  t=seq(t_start,t_end,timestep)
  k1=1/2
  k2=1/3
  Km=0.5
  nr=2
  
  alpha=list()
  alpha[["1_to_2"]]=function(C,t){
    1/5
  }
  alpha[["2_to_1"]]=function(C,t){
    1/6
  }

  f=function(C,t){
    # The only thing to take care of is that we release a vector of the same
    # size as C
   S=C[[1]]
   M=C[[2]]
   O=matrix(byrow=TRUE,nrow=2,c(k1*M*(S/(Km+S)),
                                k2*M))
    return(O)
  }
  Anl=new("TransportDecompositionOperator",t_start,Inf,nr,alpha,f)
 
  
  c01=3
  c02=2
  iv=c(c01,c02)
  inputrates=new("TimeMap",t_start,t_end,function(t){return(matrix(
    nrow=nr,
    ncol=1,
    c( 2,  2)
  ))})
  #################################################################################
  # we check if we can reproduce the linear decomposition operator from the
  # nonlinear one
#  Tr=getTransferMatrix(Anl) #this is a function of C and t
#  
#          
#  #################################################################################
#  # build the two models (linear and nonlinear)
#  mod=GeneralModel( t, A,iv, inputrates, deSolve.lsoda.wrapper) 
#  modnl=GeneralNlModel( t, Anl, iv, inputrates, deSolve.lsoda.wrapper)
#
#   Ynonlin=getC(modnl) 
#  lt1=2
#  lt2=4
#    plot(t,Ynonlin[,1],type="l",lty=lt1,col=1,
#      ylab="Concentrations",xlab="Time",ylim=c(min(Ynonlin),max(Ynonlin)))
#    lines(t,Ynonlin[,2],type="l",lty=lt2,col=2)
#    legend("topleft",c("Pool 1", "Pool 2"),lty=c(lt1,lt2),col=c(1,2))

}       
)

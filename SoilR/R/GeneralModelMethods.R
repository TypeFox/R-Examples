#
# vim:set ff=unix expandtab ts=2 sw=2:
setMethod(f="GeneralModel",
  signature=c(
    "numeric",
    "ANY",
    "numeric"#,
    #"ANY"
  ),
  definition=function # a constructor for class Model
  ### This method tries to create a Model object from any combination of arguments 
  ### that can be converted into  the required set of building blocks for a model
  ### for n arbitrarily connected pools.
  (t,			##<< A vector containing the points in time where the solution is sought.
   A,			##<< something that can be converted to any of the available DecompositionOperator classes
   ivList,		##<< A vector containing the initial amount of carbon for the n pools. The length of this vector is equal to the number of pools and thus equal to the length of k. This is checked by an internal  function. 
   inputFluxes, ##<<  something that can be converted to any of the available InFlux classes
   solverfunc=deSolve.lsoda.wrapper,		##<< The function used by to actually solve the ODE system. This can be \code{\link{deSolve.lsoda.wrapper}} or any other user provided function with the same interface. 
   pass=FALSE  ##<< Forces the constructor to create the model even if it is invalid 
   )
  {
     obj=Model(t,A,ivList,inputFluxes,solverfunc,pass)
     #obj=new(Class="Model",t,DecompositionOperator(A),ivList,InFlux(inputFluxes),solverfunc,pass)
     return(obj)
     ### A model object that can be further queried. 
     ##seealso<< \code{\link{TwopParallelModel}}, \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}} 
  }
)


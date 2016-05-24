#
# vim:set ff=unix expandtab ts=2 sw=2:
setMethod(f="GeneralModel_14",
  signature=c(
    t="numeric",
    A="ANY",
    ivList="numeric",
    initialValF="ANY",
    inputFluxes="ANY",
    inputFc="ANY",
    Fc="missing",
    di="numeric",
    lambda="missing",
    solverfunc="function",
    pass="logical"
  ),
  definition=function # a constructor for class Model_14
  ### This method tries to create a Model object from any combination of arguments 
  ### that can be converted into  the required set of building blocks for a model
  ### for n arbitrarily connected pools.
(t,	##<< A vector containing the points in time where the solution is sought.
 A,	
 ivList,
 initialValF, 
 inputFluxes, 
 inputFc,
 di=-0.0001209681, 
 solverfunc=deSolve.lsoda.wrapper,
 pass=FALSE  
 )
  {
  obj=Model_14(
    t=t,
    A=A,
    ivList=ivList,
    initialValF=initialValF,
    inputFluxes=inputFluxes,
    inputFc=inputFc,
    c14DecayRate=di,
    solverfunc=solverfunc,
    pass=pass
  )
  return(obj)
  ### A model object that can be further queried. 
  ##seealso<< \code{\link{TwopParallelModel}}, \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}} 
  }
)
setMethod(f="GeneralModel_14",
  signature=c(
    t="numeric",
    A="ANY",
    ivList="numeric",
    initialValF="ANY",
    inputFluxes="ANY",
    inputFc="ANY",
    Fc="missing",
    di="numeric",
    lambda="missing",
    solverfunc="function",
    pass="missing"
  ),
  definition=function # a constructor for class Model <- 14
  ### A wrapper for \code{\link{Model_14}} with pass=FALSE
(t,	
 A,	
 ivList,
 initialValF, 
 inputFluxes, 
 inputFc,
 di=-0.0001209681, 
 solverfunc=deSolve.lsoda.wrapper,
 pass=FALSE  
 )
  {
    obj=Model_14(
      t=t,
      A=A,
      ivList=ivList,
      initialValF=initialValF,
      inputFluxes=inputFluxes,
      inputFc=inputFc,
      c14DecayRate=di,
      solverfunc=solverfunc,
      pass=FALSE
    )
     return(obj)
  }
)
setMethod(f="GeneralModel_14",
  signature=c(
    t="numeric",
    A="ANY",
    ivList="numeric",
    initialValF="ANY",
    inputFluxes="ANY",
    inputFc="ANY",
    Fc="missing",
    di="numeric",
    lambda="missing",
    solverfunc="missing",
    pass="logical"
  ),
  definition=function # a constructor for class Model_ 14
  ### A wrapper for \code{\link{Model_14}} with pass=FALSE
  (
    t,	
    A,	
    ivList,
    initialValF, 
    inputFluxes, 
    inputFc,
    di=-0.0001209681, 
    pass
   )
  {
    obj=Model_14(
      t=t,
      A=A,
      ivList=ivList,
      initialValF=initialValF,
      inputFluxes=inputFluxes,
      inputFc=inputFc,
      c14DecayRate=di,
      solverfunc=deSolve.lsoda.wrapper,
      pass=pass
    )
     return(obj)
     ### A model object that can be further queried. 
     ##seealso<< \code{\link{TwopParallelModel}}, \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}} 
  }
)

setMethod(f="GeneralModel_14",
  signature=c(
    t="numeric",
    A="ANY",
    ivList="numeric",
    initialValF="ANY",
    inputFluxes="ANY",
    inputFc="ANY",
    Fc="missing",
    di="numeric",
    lambda="missing",
    solverfunc="missing",
    pass="missing"
  ),
  definition=function # a constructor for class Model
  ### A wrapper for \code{\link{Model_14}} with pass=FALSE
  (
    t,	##<< A vector containing the points in time where the solution is sought.
    A,	
    ivList,
    initialValF, 
    inputFluxes, 
    inputFc,
    di=-0.0001209681, 
    pass=FALSE  
  )
  {
    obj=Model_14(
      t=t,
      A=A,
      ivList=ivList,
      initialValF=initialValF,
      inputFluxes=inputFluxes,
      inputFc=inputFc,
      c14DecayRate=di,
      solverfunc=deSolve.lsoda.wrapper,
      pass=FALSE
    )
    return(obj)
    ### A model object that can be further queried. 
    ##seealso<< \code{\link{TwopParallelModel}}, \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}} 
  }
)

setMethod(f="GeneralModel_14",
  signature=c(
    t="numeric",
    A="ANY",
    ivList="numeric",
    initialValF="ANY",
    inputFluxes="ANY",
    inputFc="missing",
    Fc="ANY",
    di="numeric",
    lambda="missing",
    solverfunc="missing",
    pass="missing"
  ),
  definition=function # a constructor for class Model
  ### A wrapper for \code{\link{Model_14}} with pass=FALSE
  (
    t,	##<< A vector containing the points in time where the solution is sought.
    A,	
    ivList,
    initialValF, 
    inputFluxes, 
    Fc,
    di=-0.0001209681, 
    pass=FALSE  
  )
  {
    warning("The argument Fc has been renamed to inputFc. The use of Fc is deprecated and will be removed in the next version. 
    Please change your code to use the new keyword argument inputFc.")
    obj=Model_14(
      t=t,
      A=A,
      ivList=ivList,
      initialValF=initialValF,
      inputFluxes=inputFluxes,
      inputFc=Fc,
      c14DecayRate=di,
      solverfunc=deSolve.lsoda.wrapper,
      pass=FALSE
    )
    return(obj)
    ### A model object that can be further queried. 
    ##seealso<< \code{\link{TwopParallelModel}}, \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}} 
  }
)
setMethod(f="GeneralModel_14",
  signature=c(
    t="numeric",
    A="ANY",
    ivList="numeric",
    initialValF="ANY",
    inputFluxes="ANY",
    inputFc="ANY",
    Fc="missing",
    di="missing",
    lambda="missing",
    solverfunc="missing",
    pass="logical"
  ),
  definition=function # a constructor for class Model
  ### A wrapper for \code{\link{Model_14}} with pass=FALSE
  (
    t,	##<< A vector containing the points in time where the solution is sought.
    A,	
    ivList,
    initialValF, 
    inputFluxes, 
    inputFc,
    pass
  )
  {
    warning("The argument Fc has been renamed to inputFc. The use of Fc is deprecated and will be removed in the next version. 
    Please change your code to use the new keyword argument inputFc.")
    obj=Model_14(
      t=t,
      A=A,
      ivList=ivList,
      initialValF=initialValF,
      inputFluxes=inputFluxes,
      inputFc=inputFc,
      c14DecayRate=-0.0001209681,
      solverfunc=deSolve.lsoda.wrapper,
      pass=pass
    )
    return(obj)
    ### A model object that can be further queried. 
    ##seealso<< \code{\link{TwopParallelModel}}, \code{\link{TwopSeriesModel}}, \code{\link{TwopFeedbackModel}} 
  }
)

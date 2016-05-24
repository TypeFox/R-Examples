#' Structural model: one-compartment, oral absorption, multiple bolus dose, parameterized using KE.
#' 
#' This is a structural model function that encodes a  model that is 
#' one-compartment, oral absorption, multiple bolus dose, parameterized using KE.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#' 
#' @param model_switch a vector of values, the same size as \code{xt}, identifying which model 
#' response should be computed for the 
#' corresponding xt value.  Used for multiple response models.
#' @param xt a vector of independent variable values (often time).
#' @param parameters A named list of parameter values.
#' @param poped.db a poped database.  This can be used to extract information that may be needed in the model file.
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.md.KE.R
#' 
#' @export

ff.PK.1.comp.oral.md.KE <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  ## -- Analytic solution for both mutiple and single dosing
  with(as.list(parameters),{
    y=xt
    N = floor(xt/TAU)+1
    y=(DOSE*Favail/V)*(KA/(KA - KE)) * 
      (exp(-KE * (xt - (N - 1) * TAU)) * (1 - exp(-N * KE * TAU))/(1 - exp(-KE * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
    return(list( y=y,poped.db=poped.db))
  })
}

#' Structural model: one-compartment, oral absorption, multiple bolus dose, parameterized using CL.
#' 
#' This is a structural model function that encodes a  model that is 
#' one-compartment, oral absorption, multiple bolus dose, parameterized using CL.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.md.CL.R
#' 
#' @export
ff.PK.1.comp.oral.md.CL <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  ## -- Analytic solution for both mutiple and single dosing
  with(as.list(parameters),{
    y=xt
    N = floor(xt/TAU)+1
    y=(DOSE*Favail/V)*(KA/(KA - CL/V)) * 
      (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
    return(list( y=y,poped.db=poped.db))
  })
}

#' Structural model: one-compartment, oral absorption, single bolus dose, parameterized using KE.
#' 
#' This is a structural model function that encodes a  model that is 
#' one-compartment, oral absorption, single bolus dose, parameterized using KE.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.sd.KE.R
#' 
#' @export
ff.PK.1.comp.oral.sd.KE <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-KE)))*(exp(-KE*xt)-exp(-KA*xt))
    return(list( y= y,poped.db=poped.db))
  })
}

#' Structural model: one-compartment, oral absorption, single bolus dose, parameterized using CL.
#' 
#' This is a structural model function that encodes a  model that is 
#' one-compartment, oral absorption, single bolus dose, parameterized using CL.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.sd.CL.R
#' 
#' @export
ff.PK.1.comp.oral.sd.CL <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption
  with(as.list(parameters),{
    y=xt
    y=(DOSE*Favail*KA/(V*(KA-CL/V)))*(exp(-CL/V*xt)-exp(-KA*xt))
    return(list( y= y,poped.db=poped.db))
  })
}

#' Structural model: one-compartment, single bolus IV dose, parameterized using CL driving an EMAX model with a direct efect.
#' 
#' This is a structural model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PKPD.1.comp.sd.CL.emax.R
#' 
#' @export
ff.PKPD.1.comp.sd.CL.emax <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    y=xt
    MS <- model_switch
    
    # PK model
    CONC = DOSE/V*exp(-CL/V*xt) 
    
    # PD model
    EFF = E0 + CONC*EMAX/(EC50 + CONC)
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return(list( y= y,poped.db=poped.db))
  })
}

#' Structural model: one-compartment, oral absoprtion, multiple bolus dose, 
#' parameterized using CL driving an inhibitory IMAX model with a direct efect.
#' 
#' This is a structural model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{ff_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family structural_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PKPD.1.comp.oral.md.CL.imax.R
#' 
#' @export
ff.PKPD.1.comp.oral.md.CL.imax <- function(model_switch,xt,parameters,poped.db){
  ##-- Model: One comp first order absorption + inhibitory imax
  ## -- works for both mutiple and single dosing  
  with(as.list(parameters),{
    
    y=xt
    MS <- model_switch
    
    # PK model
    returnArgs=ff.PK.1.comp.oral.md.CL(model_switch,xt,parameters,poped.db)
    CONC=returnArgs$y
    
    # PD model
    EFF = E0*(1 - CONC*IMAX/(IC50 + CONC))
    
    y[MS==1] = CONC[MS==1]
    y[MS==2] = EFF[MS==2]
    
    return(list( y= y,poped.db=poped.db))
  })
}

#' RUV model:  
#' Additive and Proportional.
#' 
#' This is a residual unexplained variability (RUV) model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{fError_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' @param epsi A matrix with the same number of rows as the \code{xt} vector, colums match the numbers defined in this 
#' function.
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family RUV_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.md.CL.R
#' @export
feps.add.prop <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Additive + Proportional 
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  y = y*(1+epsi[,1])+epsi[,2]
  
  return(list( y= y,poped.db =poped.db )) 
}

#' RUV model:  
#' Additive .
#' 
#' This is a residual unexplained variability (RUV) model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{fError_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' @param epsi A matrix with the same number of rows as the \code{xt} vector, colums match the numbers defined in this 
#' function.
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family RUV_models
#' 
#' @example tests/testthat/examples_fcn_doc/examples_feps.add.R
#' @export
feps.add <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Additive 
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  y = y+epsi[,1]
  
  return(list( y= y,poped.db =poped.db )) 
}

#' RUV model:  
#' Proportional.
#' 
#' This is a residual unexplained variability (RUV) model function that encodes the model described above.
#' The function is suitable for input to the \code{\link{create.poped.database}} function using the
#'  \code{fError_file} argument.
#'
#' @inheritParams ff.PK.1.comp.oral.md.KE
#' @param epsi A matrix with the same number of rows as the \code{xt} vector, colums match the numbers defined in this 
#' function.
#' 
#' @return A list consisting of:
#' \enumerate{
#' \item y the values of the model at the specified points.
#' \item poped.db A (potentially modified) poped database.
#' }
#' 
#' @family models 
#' @family RUV_models
#' 
#' @example tests/testthat/examples_fcn_doc/warfarin_basic.R
#' @example tests/testthat/examples_fcn_doc/examples_ff.PK.1.comp.oral.sd.CL.R
#' 
#' @export
feps.prop <- function(model_switch,xt,parameters,epsi,poped.db){
  ## -- Residual Error function
  ## -- Proportional 
  returnArgs <- do.call(poped.db$model$ff_pointer,list(model_switch,xt,parameters,poped.db)) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  y = y*(1+epsi[,1])
  
  return(list( y= y,poped.db =poped.db )) 
}


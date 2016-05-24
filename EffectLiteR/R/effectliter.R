


#' Estimate average and conditional effects
#' 
#' This function is the main funtion of the package and can be used to estimate
#' average and conditional effects of a treatment variable on an outcome variable,
#' taking into account any number of continuous and categorical covariates.
#' It automatically generates lavaan syntax for a multi-group structural equation
#' model, runs the model using lavaan, and extracts various average and conditional
#' effects of interest.
#' 
#' @param y Dependent variable (character string). Can be the name of a manifest variable or of a latent variable.
#' @param x Treatment variable (character string) treated as categorical variable.
#' @param k Vector of manifest variables treated as categorical covariates (character vector).
#' @param z Vector of continuous covariates (character vector). Names of both manifest and latent variables are allowed.
#' @param control Value of \code{x} that is used as control group.
#' @param measurement Measurement model. The measurement model is lavaan syntax (character string), that will be appended before the automatically generated lavaan input. It can be used to specify a measurement for a latent outcome variable and/or latent covariates. See also the example and \code{\link[EffectLiteR]{generateMeasurementModel}}.
#' @param data A data frame. 
#' @param fixed.cell logical. If \code{FALSE} (default), the group sizes are treated as stochastic rather than fixed.
#' @param fixed.z logical. If \code{FALSE} (default), the continuous covariates are treated as stochastic rather than fixed. fixed.z 
#' @param missing Missing data handling. Will be passed on to \code{\link[lavaan]{sem}}.
#' @param se Type of standard errors. Will be 
#' passed on to \code{\link[lavaan]{sem}}.
#' @param bootstrap Number of bootstrap draws, if bootstrapping is used. Will be 
#' passed on to \code{\link[lavaan]{sem}}.
#' @param mimic Will be passed on to \code{\link[lavaan]{sem}}.
#' @param syntax.only logical. If \code{TRUE}, only syntax is returned and the model 
#' will not be estimated.
#' @param interactions character. Can be one of \code{c("all","none","2-way","X:K","X:Z")} and indicates the type of interaction used in the parameterization of the regression.
#' @param propscore Vector of covariates (character vector) that will be used to compute (multiple) propensity scores based on a multinomial regression without interactions. Alternatively, the user can specify a formula with the treatment variable as dependent variable for more control over the propensity score model.
#' @param ids Formula specifying cluster ID variables. Will be passed on to \code{\link[lavaan.survey]{lavaan.survey}}. See \code{\link[survey]{svydesign}} for details.
#' @param weights Formula to specify sampling weights. Currently only one weight variable is supported. Will be passed on to \code{\link[lavaan.survey]{lavaan.survey}}. See \code{\link[survey]{svydesign}} for details. Note: Only use weights if you know what you are doing. For example, some conditional treatment effects may require different weights than average effects.
#' @param homoscedasticity logical. If \code{TRUE}, residual variances of the dependent variable are assumed to be homogeneous across cells.
#' @param add Character string that will be pasted at the end of the generated lavaan syntax. Can for example be used to add additional (in-) equality constraints or to compute user-defined conditional effects.
#' @param ... Further arguments passed to \code{\link[lavaan]{sem}}. Currently not used.
#' @return Object of class effectlite.
#' @examples
#' ## Example with one categorical covariate
#' m1 <- effectLite(y="y", x="x", k="z", control="0", data=nonortho)
#' print(m1) 
#' 
#' ## Example with one categorical and one continuous covariate
#' m1 <- effectLite(y="dv", x="x", k=c("k1"), z=c("z1"), control="control", data=example01)
#' print(m1)
#' 
#' ## Example with latent outcome and latent covariate
#' measurement <- '
#' eta2 =~ 1*CPM12 + 1*CPM22
#' eta1 =~ 1*CPM11 + 1*CPM21
#' CPM11 + CPM12 ~ 0*1
#' CPM21 ~ c(m,m)*1
#' CPM22 ~ c(p,p)*1'
#'
#' m1 <- effectLite(y="eta2", x="x", z=c("eta1"), control="0", 
#'                  measurement=measurement, data=example02lv)
#' print(m1)
#' 
#'\dontrun{
#' ## Example with cluster variable and sampling weights
#' m1 <- effectLite(y="y", x="x", z="z", fixed.cell=TRUE, control="0", 
#'                     syntax.only=F, data=example_multilevel, 
#'                     ids=~cid, weights=~weights)
#' print(m1)
#' }
#' @export
#' @import lavaan
effectLite <- function(y, x, k=NULL, z=NULL, control="0", 
                       measurement=character(), data, fixed.cell=FALSE, 
                       fixed.z=FALSE, missing="listwise", se="standard", 
                       bootstrap=1000L, mimic="lavaan", syntax.only=FALSE, interactions="all", 
                       propscore=NULL, ids=~0, weights=NULL, 
                       homoscedasticity=FALSE, add=character(),...){
  
  obj <- new("effectlite")
  ##TODO make use of ldots argument: sem_args <- list(...) and then pass it to sem() 
  ##TODO change such that first class input is generated, then class parnames...
  obj@call <- match.call()
  obj@input <- createInput(y,x,k,z,propscore,control,measurement,data, 
                           fixed.cell, fixed.z, missing, se, bootstrap,
                           mimic, interactions, ids, weights, homoscedasticity,
                           add)
  obj@input <- computePropensityScore(obj@input)
  obj@parnames <- createParNames(obj)  
  obj@lavaansyntax <- createLavaanSyntax(obj)
  
  if(syntax.only){
    res <- obj@lavaansyntax@model    
  }else{
    obj@results <- computeResults(obj)
    res <- obj
  }
  
  return(res)  
}



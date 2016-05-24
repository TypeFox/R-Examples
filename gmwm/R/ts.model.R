# Copyright (C) 2014 - 2015  James Balamuta, Stephane Guerrier, Roberto Molinari
#
# This file is part of GMWM R Methods Package
#
# The `gmwm` R package is free software: you can redistribute it and/or modify it
# under the terms of the Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
# included within the packages source as the LICENSE file.
#
# The `gmwm` R package is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the Attribution-NonCommercial-ShareAlike 4.0 International 
# (CC BY-NC-SA 4.0) along with `gmwm`.  If not, see <http://www.smac-group.com/licensing/>.

#' @title Create an Autoregressive 1 [AR(1)] Process
#' @description Setups the necessary backend for the AR1 process.
#' @param phi A \code{double} value for the \eqn{\phi}{phi} of an AR1 process.
#' @param sigma2 A \code{double} value for the variance, \eqn{\sigma ^2}{sigma^2}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{Used in summary: "AR1","SIGMA2"}
#'  \item{theta}{\eqn{\phi}{phi}, \eqn{\sigma^2}{sigma^2}}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{"AR1"}
#'  \item{obj.desc}{Depth of Parameters e.g. list(1,1)}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' AR1()
#' AR1(phi=.32, sigma2=1.3)
AR1 = function(phi = NULL, sigma2 = 1) {
  starting = FALSE;
  if(is.null(phi)){
    phi = 0;
    sigma2 = 1;
    starting = TRUE;
  }
  if(length(phi) != 1 & length(sigma2) != 1){
    stop("Bad AR1 model submitted. Must be double values for two parameters.")
  }
  out = structure(list(process.desc = c("AR1","SIGMA2"),
                       theta = c(phi,sigma2),
                       plength = 2,
                       desc = "AR1",
                       obj.desc = list(c(1,1)),
                       starting = starting), class = "ts.model")
  invisible(out)
}

#' @title Create a Gauss-Markov (GM) Process
#' @description Setups the necessary backend for the GM process.
#' @param beta A \code{double} value for the \eqn{\beta}{beta} of an GM process.
#' @param sigma2_gm A \code{double} value for the variance, \eqn{\sigma ^2_{gm}}{sigma^2[gm]}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{Used in summary: "BETA","SIGMA2"}
#'  \item{theta}{\eqn{\beta}{beta}, \eqn{\sigma ^2_{gm}}{sigma^2[gm]}}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{"GM"}
#'  \item{obj.desc}{Depth of Parameters e.g. list(1,1)}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @details 
#' When supplying values for \eqn{\beta}{beta} and \eqn{\sigma ^2_{gm}}{sigma^2[gm]},
#' these parameters should be of a GM process and NOT of an AR1. That is,
#' do not supply AR1 parameters such as \eqn{\phi}{phi}, \eqn{\sigma^2}{sigma^2}.
#' 
#' Internally, GM parameters are converted to AR1 using the `freq` 
#' supplied when creating data objects (\link[gmwm]{imu}, \link[gmwm]{gts})
#' or specifying a `freq` parameter in \link[gmwm]{gmwm} or \link[gmwm]{gmwm.imu}.
#' 
#' The `freq` of a data object takes precedence over the `freq` set when modeling.
#' @author JJB
#' @examples
#' GM()
#' GM(beta=.32, sigma2_gm=1.3)
GM = function(beta = NULL, sigma2_gm = 1) {
  starting = FALSE;
  if(is.null(beta)){
    beta = 0;
    sigma2_gm = 1;
    starting = TRUE;
  }
  if(length(beta) != 1 & length(sigma2_gm) != 1){
    stop("Bad GM model submitted. Must be double values for two parameters.")
  }
  out = structure(list(process.desc = c("BETA","SIGMA2_GM"),
                       theta = c(beta,sigma2_gm),
                       plength = 2,
                       desc = "GM",
                       obj.desc = list(c(1,1)),
                       starting = starting), class = "ts.model")
  invisible(out)
}


#' @title Create an Autoregressive P [AR(P)] Process
#' @description Setups the necessary backend for the AR(P) process.
#' @param phi A \code{vector} with double values for the \eqn{\phi}{phi} of an AR(P) process.
#' @param sigma2 A \code{double} value for the variance, \eqn{\sigma ^2}{sigma^2}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{Used in summary: "AR-1","AR-2", ..., "AR-P", "SIGMA2"}
#'  \item{theta}{\eqn{\phi_1}{phi[[1]]}, \eqn{\phi_2}{phi[[2]]}, ..., \eqn{\phi_p}{phi[[p]]}, \eqn{\sigma^2}{sigma^2}}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{"AR"}
#'  \item{obj.desc}{Depth of Parameters e.g. list(p,1)}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' AR(1) # Slower version of AR1()
#' AR(phi=.32, sigma=1.3) # Slower version of AR1()
#' AR(2) # Equivalent to ARMA(2,0).
AR = function(phi = NULL, sigma2 = 1) {
  starting = FALSE;
  
  if(is.null(phi)){
    stop("Must supply either a whole number of a string of numbers for phi parameter in AR().")
  }
  
  p = length(phi)
  if(p == 1 & is.whole(phi)){
   p = phi
   phi = rep(0,p)
   starting = TRUE;
  }

  out = structure(list(process.desc = c(paste0("AR-",1:p), "SIGMA2"),
                      theta = c(phi,sigma2),
                      plength = 2,
                      desc = "ARMA", #update to AR when backend supports it!
                      obj.desc = list(c(p,0,1)), # Remove the 0, when the backend supports it!
                      starting = starting), class = "ts.model")
  invisible(out)
}


#' @title Create an Moving Average Q [MA(Q)] Process
#' @description Setups the necessary backend for the MA(Q) process.
#' @param theta A \code{vector} with double values for the \eqn{\theta}{theta} of an MA(Q) process.
#' @param sigma2 A \code{double} value for the variance, \eqn{\sigma ^2}{sigma^2}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{Used in summary: "MA-1","MA-2", ..., "MA-Q", "SIGMA2"}
#'  \item{theta}{\eqn{\theta_1}{theta[[1]]}, \eqn{\theta_2}{theta[[2]]}, ..., \eqn{\theta_q}{theta[[q]]}, \eqn{\sigma^2}{sigma^2}}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{"MA"}
#'  \item{obj.desc}{Depth of Parameters e.g. list(q,1)}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' MA(1) # One theta
#' MA(2) # Two thetas!
#' 
#' MA(theta=.32, sigma=1.3) # 1 theta with a specific value.
#' MA(theta=c(.3,.5), sigma=.3) # 2 thetas with specific values.
MA = function(theta = NULL, sigma2 = 1) {
  starting = FALSE;
  if(is.null(theta)){
    stop("Must supply either a whole number of a string of doubles for theta parameter in MA().")
  }
  
  q = length(theta)
  if(q == 1 & is.whole(theta)){
    q = theta
    theta = rep(1,q)
    starting = TRUE;
  }
  
  out = structure(list(process.desc = c(paste0("MA-",1:q), "SIGMA2"),
                      theta = c(theta,sigma2),
                      plength = 2,
                      desc = "ARMA", # Move to MA when backend supports it!
                      obj.desc = list(c(0,q,1)), # Remove the 0 when backend supports it!
                      starting = starting), class = "ts.model")
  invisible(out)
}


#' @title Create an Quantisation Noise (QN) Process
#' @description Sets up the necessary backend for the QN process.
#' @param q2 A \code{double} value for the \eqn{Q^2}{Q^2} of a QN process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{Used in summary: "QN"}
#'  \item{theta}{\eqn{Q^2}{Q^2}}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{y desc replicated x times}
#'  \item{obj.desc}{Depth of Parameters e.g. list(1)}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' QN()
#' QN(q2=3.4)
QN = function(q2 = NULL) {
  starting = FALSE
  if(is.null(q2)){
    q2 = 2
    starting = TRUE
  }
  if(length(q2) != 1){
    stop("Bad QN model submitted. Must be a double that indicates the Q2 value.")
  }
  out = structure(list(process.desc = "QN",
                       theta = q2,
                       plength = 1,
                       desc = "QN",
                       obj.desc = list(1),
                       starting = starting), class = "ts.model")
  invisible(out)
}

#' @title Create an White Noise (WN) Process
#' @description Sets up the necessary backend for the WN process.
#' @param sigma2 A \code{double} value for the variance, \eqn{\sigma ^2}{sigma^2}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{Used in summary: "WN"}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{y desc replicated x times}
#'  \item{obj.desc}{Depth of Parameters e.g. list(1)}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' WN()
#' WN(sigma=3.4)
WN = function(sigma2 = NULL) {
  starting = FALSE
  if(is.null(sigma2)){
    sigma2 = 3
    starting = TRUE
  }
  if(length(sigma2) != 1){
    stop("Bad WN model submitted. Must be a double that indicates the standard deviation.")
  }
  out = structure(list(process.desc = "WN",
                       theta = sigma2,
                       plength = 1,
                       desc = "WN",
                       obj.desc = list(1),
                       starting = starting), class = "ts.model")
  invisible(out)
}

#' @title Create an Random Walk (RW) Process
#' @description Sets up the necessary backend for the RW process.
#' @param sigma2 A \code{double} value for the variance, \eqn{\sigma ^2}{sigma^2}, of a WN process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{Used in summary: "RW"}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{y desc replicated x times}
#'  \item{obj.desc}{Depth of Parameters e.g. list(1)}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' RW()
#' RW(sigma=3.4)
RW = function(sigma2 = NULL) {
  starting = FALSE
  if(is.null(sigma2)){
    sigma2 = 4
    starting = TRUE
  }
  if(length(sigma2) != 1){
    stop("Bad RW model submitted. Must be a double that indicates the standard deviation.")
  }
  out = structure(list(process.desc = "RW",
                       theta = sigma2,
                       plength = 1,
                       desc = "RW",
                       obj.desc = list(1),
                       starting = starting), class = "ts.model")
  invisible(out)
}

#' @title Create an Drift (DR) Process
#' @description Sets up the necessary backend for the DR process.
#' @param slope A \code{double} value for the slope of a DR process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{Used in summary: "DR"}
#'  \item{theta}{slope}
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(1)}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @examples
#' DR()
#' DR(slope=3.4)
DR = function(slope = NULL) {
  starting = FALSE
  if(is.null(slope)){
    slope = 5
    starting = TRUE
  }
  if(length(slope) != 1){
    stop("Bad Drift model submitted. Must be a double that indicates a slope.")
  }
  out = structure(list(process.desc = "DR",
                       theta = slope,
                       plength = 1,
                       desc = "DR",
                       obj.desc = list(1),
                       starting = starting), class = "ts.model")
  invisible(out)
}

#' @title Create an Autoregressive Moving Average (ARMA) Process
#' @description Sets up the necessary backend for the ARMA process.
#' @param ar A \code{vector} or \code{integer} containing either the coefficients for \eqn{\phi}{phi}'s or the process number \eqn{p} for the Autoregressive (AR) term.
#' @param ma A \code{vector} or \code{integer} containing either the coefficients for \eqn{\theta}{theta}'s or the process number \eqn{q} for the Moving Average (MA) term.
#' @param sigma2 A \code{double} value for the standard deviation, \eqn{\sigma}{sigma}, of the ARMA process.
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{\eqn{AR*p}{AR x p}, \eqn{MA*q}{MA x q}}
#'  \item{theta}{\eqn{\sigma}{sigma}}
#'  \item{plength}{Number of Parameters}
#'  \item{obj.desc}{y desc replicated x times}
#'  \item{obj}{Depth of Parameters e.g. list(c(length(ar),length(ma),1) )}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @details
#' A standard deviation is required since the model generation statements utilize 
#' randomization functions expecting a standard deviation instead of a variance.
#' @author JJB
#' @examples
#' # Create an ARMA(1,2) process
#' ARMA(ar=1,2)
#' # Creates an ARMA(3,2) process with predefined coefficients.
#' ARMA(ar=c(0.23,.43, .59), ma=c(0.4,.3))
#' 
#' # Creates an ARMA(3,2) process with predefined coefficients and standard deviation
#' ARMA(ar=c(0.23,.43, .59), ma=c(0.4,.3), sigma2 = 1.5)
ARMA = function(ar = 1, ma = 1, sigma2 = 1.0) {
  # Assume the user specified data
  starting = FALSE
  
  # Get initial parameters
  p = length(ar)
  q = length(ma)
  
  # If P or Q == 1, this implies we might have a starting guess. 
  if( p == 1 || q == 1 ){
    if(p == 1){
      if(is.whole(ar) & ar != 0){
        ar = rep(-1, ar)
        starting = TRUE
      }else if(ar == 0){
        ar = numeric(0) # creates a size 0 vector
      }
    }
    
    if(q == 1){
      if(is.whole(ma) & ma != 0){
        ma = rep(-2, ma)
        starting = TRUE
      }else if(ma == 0){
        ma = numeric(0) 
      }
    }
  }
  
  # Update the values.
  p = length(ar)
  q = length(ma)
  
  out = structure(list(process.desc = c(rep("AR", p), rep("MA",q), "SIGMA2"),
                       theta = c(ar, ma, sigma2),
                       plength = p + q + 1,
                       desc = "ARMA",
                       obj.desc = list(c(p,q,1)),
                       starting = starting), class = "ts.model")
  invisible(out)
}

#' @title Multiple a ts.model by constant
#' @description Sets up the necessary backend for creating multiple model objects.
#' @method * ts.model
#' @param x A \code{numeric} value
#' @param y A \code{ts.model} object
#' @return An S3 object with called ts.model with the following structure:
#' \describe{
#'  \item{process.desc}{y desc replicated x times}
#'  \item{theta}{y theta replicated x times}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{y desc replicated x times}
#'  \item{obj.desc}{Depth of Parameters e.g. list(c(1,1),c(1,1))}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @keywords internal
#' @export
#' @examples
#' 4*DR()+2*WN()
#' DR()*4 + WN()*2
#' AR1(phi=.3,sigma=.2)*3
"*.ts.model" = function(x, y) {
  # Handles the ts.model()*c case
  if(!is.numeric(x)){
    temp = x
    x = y
    y = temp
  }
  out = structure(list(process.desc = rep(y$process.desc,x),
                       theta = rep(y$theta,x),
                       plength = y$plength*x,
                       desc = rep(y$desc,x),
                       obj.desc = rep(y$obj.desc,x),
                       starting = y$starting), class = "ts.model")
  invisible(out)
}

#' @title Add ts.model objects together
#' @description Sets up the necessary backend for combining ts.model objects.
#' @method + ts.model
#' @param x A \code{ts.model} object
#' @param y A \code{ts.model} object
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{process.desc}{combined x, y desc}
#'  \item{theta}{combined x, y theta}
#'  \item{plength}{Number of Parameters}
#'  \item{desc}{Add process to queue e.g. c("AR1","WN")}
#'  \item{obj.desc}{Depth of Parameters e.g. list(1, c(1,1), c(length(ar),length(ma),1) )}
#'  \item{starting}{Guess Starting values? TRUE or FALSE (e.g. specified value)}
#' }
#' @author JJB
#' @export
#' @keywords internal
#' @examples
#' DR()+WN()
#' AR1(phi=.3,sigma=.2)
"+.ts.model" = function(x, y) {
  starting = FALSE
  if(y$starting & x$starting){
    starting = TRUE
  }
  out = structure(list(process.desc = c(x$process.desc, y$process.desc),
                       theta = c(x$theta,y$theta),
                       plength = x$plength + y$plength,
                       desc = c(x$desc, y$desc),
                       obj.desc = c(x$obj.desc, y$obj.desc),
                       starting = starting), class = "ts.model")
  invisible(out)
}

#' @title Multiple a ts.model by constant
#' @description Sets up the necessary backend for creating multiple model objects.
#' @method print ts.model
#' @export
#' @param x A \code{numeric} value
#' @param ... further arguments passed to or from other methods.
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}
#'  \item{theta}
#' }
#' @keywords internal
#' @author JJB
#' @examples
#' QN() + DR() + WN() + RW() + AR1() + ARMA(1,2)
#' AR1(phi=.9,sigma2=.1) + WN(sigma2=1) + 
#' RW(sigma2=.3) + DR(slope=.5) + QN(q2=.9) + ARMA(ar=c(.3,.1),ma=c(.3,.2), sigma2= .99)
#' 
#' AR1(.9,.1) + WN(1) + RW(.3) + DR(.5) + QN(.9) + ARMA(c(.3,.1),c(.3,.2), .99)
print.ts.model = function(x, ...){

  desctable = data.frame("Terms" = x$process.desc, "Initial Values" = x$theta, stringsAsFactors = FALSE);
  cat("\nGuess Starting Values:", x$starting, "\n")
  if(x$starting){
    cat("The program will attempt to guess starting values for...\n")
    print(desctable[,1], row.names = FALSE)
    cat("To have the option of using your own starting values, please supply values for each parameter.\n")
  }else{
    print(desctable, row.names = FALSE)
    cat("The model will be initiated using the initial values you supplied.\n")
  }
}

#' @title Create a ts.model from desc string
#' @description Sets up the necessary backend for using Cpp functions to build R ts.model objects
#' @param desc A \code{character} vector containing: \code{"AR1"},\code{"DR"},\code{"WN"},\code{"RW"},\code{"QN"}
#' @return An S3 object with called ts.model with the following structure:
#' \itemize{
#'  \item{desc}
#'  \item{theta}
#' }
#' @author JJB
#' @keywords internal
#' @examples
#' desc.to.ts.model(c("AR1","WN"))
desc.to.ts.model = function(desc){
  theta = .Call('gmwm_model_theta', PACKAGE = 'gmwm', desc)
  
  out = structure(list(process.desc = .Call('gmwm_model_process_desc', PACKAGE='gmwm', desc),
                       theta = theta,
                       plength = length(theta),
                       desc = desc,
                       obj.desc = .Call('gmwm_model_objdesc', PACKAGE = 'gmwm', desc),
                       starting = TRUE), class = "ts.model")
  invisible(out)
}
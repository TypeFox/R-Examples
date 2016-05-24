#'@title  Fitting Linear Models with Endogenous Regressors using Lewbel's Higher Moments Approach
#'@aliases hmlewbel
# Description
#'@description  Fits linear models with one endogenous regressor using internal instruments built using the approach described in 
#' Lewbel A. (1997). This is a statistical technique to address the endogeneity problem where no external instrumental
#' variables are needed. The implementation allows the incorporation of external instruments if available. 
#' An important assumption for identification is that the endogenous variable has a skewed distribution.  
#
# Arguments
#'@param    y   the vector or matrix containing the dependent variable.
#'@param    X   the data frame or matrix containing the exogenous regressors of the model.
#'@param    P   the endogenous variables of the model as columns of a matrix or dataframe.  
#'@param    G   the functional form of G. It can take four values, \code{x2}, \code{x3},\code{lnx} or \code{1/x}. The last two forms are conditional on the values of the exogenous variables: greater than 1 or different from 0 respectively.
#'@param    IIV  stands for "internal instrumental variable". It can take six values: \code{g,gp,gy,yp,p2} or \code{y2}. Tells the function 
#'which internal instruments to be constructed from the data. See "Details" for further explanations.
#'@param EIV  stands for "external instrumental variable". It is an optional argument that lets the user specify any external variable(s) to be used as instrument(s).
#'@param data  optional data frame or list containing the variables in the model.
#

#'@details 
#'Consider the model below: 
#'\deqn{ Y_{t} = \beta_{0}+ \gamma^{'}X_{t} + \alpha P_{t}+\epsilon_{t} \hspace{0.3cm} (1) }{Y_t = b0 + gamma' * X_t + alpha * P_t + eps_t}
#' \deqn{ P_{t} = Z_{t}+\nu_{t} \hspace{2.5 cm} (2)}{P_t = Z + nu_t} 
#' 
#'The observed data consist of \eqn{Y_{t}}{Y}, \eqn{X_{t}}{X} and \eqn{P_{t}}{P_t}, while \eqn{Z_{t}}{Z_t}, \eqn{\epsilon_{t}}{eps_t}, and \eqn{\nu_{t}}{nu_t} 
#'are unobserved. The endogeneity problem arises from the correlation of \eqn{P_{t}}{P} with the structural error, \eqn{\epsilon_{t}}{eps_t}, 
#'since \eqn{E(\epsilon \nu)\neq 0}{E(eps * nu) > 0}. 
#'The requirement for the structural and measurement error is to have mean zero, but no restriction is imposed on their distribution. 
#'
#'Let \eqn{\bar{S}}{mean(S)} be the sample mean of a variable \eqn{S_{t}}{S_t} and \eqn{G_{t} = G(X_{t})}{G = G(X)} for any given function \eqn{G}{G} that 
#'has finite third own and cross moments. Lewbel(1997) proves that the following instruments can be constructed and used with 2SLS to obtain consistent estimates:
#' \deqn{ q_{1t}=(G_{t} - \bar{G})  \hspace{1.6 cm}(3a)}{q1 = G_t - mean(G)}
#' \deqn{ q_{2t}=(G_{t} - \bar{G})(P_{t}-\bar{P}) \hspace{0.3cm} (3b) }{q2 = (G_t - mean(G))(P_t -mean(P))}
#' \deqn{ q_{3t}=(G_{t} - \bar{G})(Y_{t}-\bar{Y}) \hspace{0.3cm} (3c)}{q3 = (G_t - mean(G))(Y_t - mean(Y))}
#' \deqn{ q_{4t}=(Y_{t} - \bar{Y})(P_{t}-\bar{P}) \hspace{0.3cm} (3d)}{q4 = (Y_t - mean(Y))(P_t - mean(P))}
#' \deqn{ q_{5t}=(P_{t}-\bar{P})^{2} \hspace{1.5 cm} (3e) }{q5 = (P_t - mean(P))^2}
#' \deqn{ q_{6t}=(Y_{t} - \bar{Y})^{2}\hspace{1.5 cm} (3f)}{q6 = (Y_t - mean(Y))^2}
#' 
#'Instruments in equations \code{3e} and \code{3f} can be used only when the measurement and the structural errors are symmetrically distributed.
#'Otherwise, the use of the instruments does not require any distributional assumptions for the errors. Given that the regressors \eqn{G(X) = X} 
#'are included as instruments, \eqn{G(X)} should not be linear in \eqn{X} in equation \code{3a}.
#'
#'Let small letter denote deviation from the sample mean: \eqn{s_{i} = S_{i}-\bar{S}}{s_i = S_i - mean(S)}. Then, using as instruments the variables presented in
#'equations \code{3} together with \code{1} and \eqn{X_{t}}{X}, the two-stage-least-squares estimation will provide consistent estimates for the parameters
#'in equation \code{1} under the assumptions exposed in Lewbel(1997).
#'

#Return Value
#'@return Returns an object of class \code{ivreg}, with the following components:
#'\item{coefficients}{ parameters estimates.}
#'\item{residulas}{ a vector of residuals.}
#'\item{fitted.values}{ a vector of predicted means.}
#'\item{n}{ number of observations.}
#'\item{df.residual}{ residual degrees of freedom for the fitted model.}
#'\item{cov.unscaled}{ unscaled covariance matrix for coefficients.}
#'\item{sigma}{ residual standard error.}
#'\item{call}{ the original function call.}
#'\item{formula}{ the model formula.}
#'\item{terms}{ a list with elements "regressors" and "instruments" containing the terms objects for the respective components.}
#'\item{levels}{ levels of the categorical regressors.}
#'\item{contrasts}{ the contrasts used foe categorical regressors.}
#'\item{x}{ a list with elements "regressors", "instruments", "projected", containing the model matrices from the respective components.
#' "projected" is the matrix of regressors projected on the image of the instruments.}
#'@keywords endogenous
#'@keywords latent
#'@keywords instruments
#'@author The implementation of the model formula by Raluca Gui based on the paper of Lewbel (1997).
#'@references  Lewbel, A. (1997). Constructing Instruments for Regressions with Measurement Error when No Additional Data Are Available,
#' with An Application to Patents and R&D. \emph{Econometrica}, \bold{65(5)}, 1201-1213.
#'@seealso \code{\link{internalIV}}, \code{\link[AER]{ivreg}}, \code{\link{liv}}
#'@examples 
#'#load data 
#'data(dataHMLewbel)
#'y <- dataHMLewbel$y
#'X <- cbind(dataHMLewbel$X1,dataHMLewbel$X2)
#'colnames(X) <- c("X1","X2")
#'P <- dataHMLewbel$P
#'
#'# call hmlewbel with internal instrument yp = (Y - mean(Y))(P - mean(P))
#'hmlewbel(y,X,P, G = "x2", IIV = "yp")  
#'
#'# build an additional instrument p2 = (P - mean(P))^2  using the internalIV() function 
#'eiv <- internalIV(y,X,P, G="x2", IIV = "p2")
#'
#'# use the additional variable as external instrument in hmlewbel()
#'h <- hmlewbel(y,X,P,G = "x2",IIV = "yp", EIV=eiv) 
#'h$coefficients
#'
#'# get the robust standard errors using robust.se() function from package ivpack
#'# library(ivpack)
#'# sder <- robust.se(h)
#'
# make availble to the package users
#'@export


hmlewbel <- function(y,X,P, G = c("x2","x3","lnx","1/x"), IIV = c("g","gp","gy","yp","p2","y2"), EIV=NULL, data=NULL){
  
  # check to see if any external instruments were provided
  if (!is.null(EIV)) {
        EIV <- as.matrix(EIV)
    }
  
 # computes the internal IVs proposed by Lewbel 97 - g, gp, gy, yp, y2, p2, depending on 
# the values provided by the user in IIV
  IV <- internalIV(y,X,P,G,IIV, data)
  
  # check if external instruments were provided. If, yes, add them to the IIV
  if (is.null(EIV)){IV1 <- IV} else {
    IV1 <- cbind(IV,EIV)
  }  
  
  # function that checks whether the model error is symmetrically distributed
  # if yes, then IIV6 can be used as instrument
  # checkAssumptions(y,X,P,IIV,EIV, data)
  
  # uses ivreg function from \pkg{AER} for users to be able to use afterwards the package ivpack 
  # hmlewbel should return an object of class "ivreg"
  dataHM <- data.frame(cbind(X,P))
  res <- AER::ivreg(y ~.|X + IV1, data = dataHM, x=TRUE )
  
  res$call <- match.call()
  class(res) <- c("ivreg")
  return(res)
    
}

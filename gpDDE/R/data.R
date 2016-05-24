#' Nicholson's Blowflies Data.
#'
#'  A set of control experiments in which the adult food supply was restricted to 0.4g per day. The adult blowfly population exhibits oscillation due to constant resource level.
#'
#' @format A data frame with two variables: \code{day} and \code{y}.
#' \code{day} indicates the time observations are made. The experiment last from the 40th day to the 315th day. Observations are made daily.
#' \code{y} is the blowfly population count.
"blowflydata"

#' Simulated dataset from a delay SIR model
#'
#' A simulated dataset from a delay SIR model specified as following:
#' \deqn{\dot{S} = \rho(t)-\beta I_{\tau}S\\
#'  \dot{I} = \beta I_{\tau}S-\gamma I}
#' The parameters are set at: \eqn{\rho(t)=4000\times(\sin(t)+2)}, \eqn{\gamma = 5},
#' \eqn{\beta = 0.0012}, and \eqn{\tau = 0.5}.  We simulated the numerical solution from time \eqn{t=0} to \eqn{t=50} and the process is then sampled regularly ten times per unit time.  Independent normal observational noise with sd = 100 is added to the numerical solution.
#'
#' @format A matrix of two columns. The first column is the observed state \code{S} and the second column is the observed state \code{I}.
"DSIRdata"

#' Fitted Delay SIR model
#'
#' Fitted results of the Delay SIR model with data of the \code{DSIR.data} and fitted using the code in the example section of \code{Profile.LS.DDE}.
#' @format A list returned by \code{Profile.LS.DDE}
"DSIR.fit"

#' Initial values for the spline coefficients for the Delay SIR model
#'
#' Initial values for the spline coefficients for the Delay SIR model used in the example of \code{Profile.LS.DDE} with data of the \code{DSIR.data}.
#' @format A matrix. Two columns for state "S" and "I" respectively.
"DSIRInitCoefs"

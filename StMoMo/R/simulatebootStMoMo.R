
#' Simulate future sample paths from a Bootstrapped Stochastic 
#' Mortality Model
#'
#' Simulate future sample paths from a Bootstrapped Stochastic Mortality Model.
#' The period indexes \eqn{\kappa_t^{(i)}, i = 1,..N,} are modelled
#' using a Multivariate Random Walk with Drift. The cohort index 
#' \eqn{\gamma_{t-x}} is modelled using an ARIMA\eqn{(p, d, q)}. By default
#' an ARIMA\eqn{(1, 1, 0)} with a constant is used.
#' 
#' @param object an object of class \code{"bootStMoMo"} with the bootstrapped 
#' parameters of a stochastic mortality model.
#' @param nsim number of sample paths to simulate from each bootstrapped 
#' sample. Thus if there are \code{nBoot} bootstrapped samples the total 
#' number of paths will be \code{nsim * nBoot}.
#' @inheritParams simulate.mrwd
#' @param oxt optional array/matrix/vector or scalar of known offset to be 
#' added in the simulations. This can be used to specify any a priori known 
#' component to be added to the simulated predictor.  
#' @inheritParams forecast.fitStMoMo
#' 
#' @return A list of class \code{"simStMoMo"} with components
#' 
#' \item{rates}{ a three dimensional array with the future simulated rates.}
#' 
#' \item{ages}{ vector of ages corresponding to the first dimension of 
#' \code{rates}.}
#' 
#' \item{years}{ vector of years for which a simulations has been produced. 
#' This corresponds to the second dimension of \code{rates}.}  
#'  
#' \item{kt.s}{ information on the simulated paths of the period indices of 
#' the model. This is a list with the simulated paths of \eqn{\kappa_t}
#'  (\code{sim}) and the \code{years} for which simulations were produced. 
#'  If the mortality model does not have any age-period terms (i.e. \eqn{N=0})
#'  this is set to \code{NULL}.}
#'   
#' \item{gc.s}{ information on the simulated paths of the cohort index of the
#'  model. This is a list with the simulated paths of \eqn{\gamma_c} 
#'  (\code{sim}) and the \code{cohorts} for which simulations were produced. 
#'  If the mortality model does not have a cohort effect this is set to 
#'  \code{NULL}.} 
#'  
#' \item{oxt.s}{ a three dimensional array with the offset used in the 
#' simulations.}
#' 
#' \item{fitted}{ a three dimensional array with the in-sample rates of the 
#' model for the years for which the mortality model was fitted 
#' (and bootstrapped).}
#'  
#'  \item{jumpchoice}{Jump-off method used in the simulation.}
#'  
#'  \item{model}{the bootstrapped model from which the simulations were 
#'  produced.}
#'  
#' @details
#' For further details see \code{\link{simulate.fitStMoMo}}. 
#' 
#' @seealso \code{\link{bootstrap.fitStMoMo}}, \code{\link{simulate.fitStMoMo}}
#' 
#' @examples
#' #Long computing times
#' \dontrun{
#' #Lee-Carter: Compare projection with and without parameter uncertainty
#' library(fanplot)
#' LCfit <- fit(lc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'              ages = EWMaleData$ages, years = EWMaleData$years)
#' LCResBoot <- bootstrap(LCfit, nBoot = 500)
#' LCResBootsim <- simulate(LCResBoot)
#' LCsim <- simulate(LCfit, nsim = 500)#' 
#' plot(LCfit$years, log(LCfit$Dxt / LCfit$Ext)["10", ], 
#'      xlim = range(LCfit$years, LCsim$years),
#'      ylim = range(log(LCfit$Dxt / LCfit$Ext)["10", ], 
#'                   log(LCsim$rates["10", , ])),
#'      type = "l", xlab = "year", ylab = "log rate", 
#'      main = "Mortality rate projection at age 10 with and without parameter uncertainty") 
#' fan(t(log(LCResBootsim$rates["10", , ])),start = LCResBootsim$years[1], 
#'     probs = c(2.5, 10, 25, 50, 75, 90, 97.5), n.fan = 4, 
#'     fan.col = colorRampPalette(c(rgb(0, 0, 1), rgb(1, 1, 1))), ln = NULL)
#' fan(t(log(LCsim$rates["10", 1:(length(LCsim$years) - 3), ])),
#'     start = LCsim$years[1], probs = c(2.5, 10, 25, 50, 75, 90, 97.5),
#'     n.fan = 4, fan.col = colorRampPalette(c(rgb(1, 0, 0), rgb(1, 1, 1))), 
#'     ln = NULL)
#' }
#' 
#' @export
simulate.bootStMoMo <-function(object, nsim = 1, seed = NULL, h = 50, 
                               oxt = NULL, gc.order = c(1, 1, 0), 
                               gc.include.constant = TRUE,
                               jumpchoice = c("fit", "actual"), 
                               kt.lookback = NULL, gc.lookback = NULL,
                               ...) {  
  jumpchoice <- match.arg(jumpchoice)
  ages <- object$model$ages
  nAges <- length(ages)
  nBoot <- length(object$bootParameters)
  nPath <- nsim * nBoot
  
  ## Handle generator seed
  if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- .Random.seed
  else {
    R.seed <- .Random.seed
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  ## Offset  
  if (is.null(oxt)) 
    oxt <- 0  
  oxt.s <- array(oxt, dim = c(nAges, h, nBoot))
  
  
  ## Do simulations  
  
  #Initialse output variables
  modelFit <- object$bootParameters[[1]]
  modelFit$Dxt <- object$model$Dxt
  modelFit$Ext <- object$model$Ext
  tempSim <- simulate.fitStMoMo(object = modelFit, 
                                nsim = nsim, seed = NULL, h = h, 
                                oxt = oxt.s[, , 1], gc.order = gc.order, 
                                gc.include.constant = gc.include.constant, 
                                jumpchoice = jumpchoice, 
                                kt.lookback = kt.lookback, 
                                gc.lookback = gc.lookback)
  rates <- array(NA, c(dim(tempSim$rates)[1:2], nPath), 
                 list(dimnames(tempSim$rates)[[1]], 
                      dimnames(tempSim$rates)[[2]], 1:nPath))
  fitted <- array(NA, c(dim(tempSim$fitted)[1:2], nPath), 
                  list(dimnames(tempSim$fitted)[[1]], 
                       dimnames(tempSim$fitted)[[2]], 1:nPath))  
  if (is.null(tempSim$kt.s)) {
    kt.s <- NULL
  } else {
    kt.s <- list(sim = NULL, years = tempSim$kt.s$years)
    kt.s$sim <- array(NA, c(dim(tempSim$kt.s$sim)[1:2], nPath), 
                      list(dimnames(tempSim$kt.s$sim)[[1]], 
                           dimnames(tempSim$kt.s$sim)[[2]], 1:nPath))
  }
  if (is.null(tempSim$gc.s)) {
    gc.s <- NULL
  } else {
    gc.s <- list(sim = NULL, cohorts = tempSim$gc.s$cohorts)
    gc.s$sim <- array(NA, c(dim(tempSim$gc.s$sim)[1], nPath), 
                      list(dimnames(tempSim$gc.s$sim)[[1]], 1:nPath))
  }
  
  # Do simulations for each bootstrap sample
  i <- 1
  while (i <= nBoot) {    
    rates[, , ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$rates
    fitted[, , ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$fitted
    if (!is.null(kt.s)) 
      kt.s$sim[, , ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$kt.s$sim  
    if (!is.null(gc.s)) 
      gc.s$sim[, ((i - 1) * nsim + 1):(i * nsim)] <- tempSim$gc.s$sim  
    i <- i + 1
    if ( i <= nBoot) {
      modelFit <- object$bootParameters[[i]]
      modelFit$Dxt <- object$model$Dxt
      modelFit$Ext <- object$model$Ext
      tempSim <- simulate.fitStMoMo(object = modelFit, 
                              nsim = nsim, seed = NULL, h = h, 
                              oxt = oxt.s[, , i], gc.order = gc.order, 
                              gc.include.constant = gc.include.constant,                               
                              jumpchoice = jumpchoice,
                              kt.lookback = kt.lookback, 
                              gc.lookback = gc.lookback)
    
    }    
  }  
  structure(list(rates = rates, ages = tempSim$ages, 
                 years = tempSim$years, kt.s = kt.s, gc.s = gc.s, 
                 oxt.s = oxt.s, fitted = fitted, jumpchoice = jumpchoice,
                 model = object, call = match.call()), 
            class = "simStMoMo")
}


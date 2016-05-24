#' Simulate future sample paths from a Stochastic Mortality Model
#'
#' Simulate future sample paths from a Stochastic Mortality Model.
#' The period indexes \eqn{\kappa_t^{(i)}, i = 1,..N,} are modelled
#' using a Multivariate Random Walk with Drift. The cohort index 
#' \eqn{\gamma_{t-x}} is modelled using an ARIMA\eqn{(p, d, q)}. By default
#' an ARIMA\eqn{(1, 1, 0)} with a constant is used.
#' 
#' @param object an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' @param nsim number of sample paths to simulate.
#' @inheritParams simulate.mrwd
#' @param oxt optional matrix/vector or scalar of known offset to be 
#' added in the simulations. This can be used to specify any a priori 
#' known component to be added to the simulated predictor.
#' @inheritParams forecast.fitStMoMo
#'
#' @return A list of class \code{"simStMoMo"} with components:
#' 
#' \item{rates}{ a three dimensional array with the future simulated rates.}
#' 
#' \item{ages}{ vector of ages corresponding to the first dimension of 
#' \code{rates}.}
#' 
#' \item{years}{vector of years for which a simulations has been produced. 
#' This corresponds to the second dimension of \code{rates}.}  
#' 
#' \item{kt.s}{ information on the simulated paths of the period indexes 
#' of the model. This is a list with the  \code{model} fitted to 
#' \eqn{\kappa_t}; the simulated paths (\code{sim}); and the \code{years} 
#' for which simulations were produced. If the mortality model does not 
#' have any age-period terms (i.e. \eqn{N=0}) this is set to \code{NULL}.}
#'   
#' \item{gc.s}{ information on the simulated paths of the cohort index of 
#' the model. This is a list with the \code{model} fitted to \eqn{\gamma_c}; 
#' the simulated paths (\code{sim}); and the \code{cohorts} for which 
#' simulations were produced. If the mortality model does not have a cohort 
#' effect this is set to \code{NULL}.} 
#' 
#' \item{oxt.s}{ a three dimensional array with the offset used in the 
#' simulations.}
#' 
#' \item{fitted}{ a three dimensional array with the in-sample rates of 
#' the model for the years for which the mortality model was fitted.}
#'  
#' \item{jumpchoice}{Jump-off method used in the simulation.}
#'  
#' \item{model}{ the model fit from which the simulations were produced.}
#' 
#' @details
#' Fitting and simulation of the Multivariate Random Walk with Drift
#' for the period indexes is done using the function \code{\link{mrwd}}.
#' Fitting and simulation of the ARIMA model for the cohort index
#' is done with function \code{\link[forecast]{Arima}} from package 
#' \pkg{forecast}. See the latter function for further details on 
#' input arguments \code{gc.order} and \code{gc.include.constant}. 
#' 
#' Note that in some cases simulations of the 
#' cohort effects may be needed for a horizon longer than \code{h}.
#' This is the case when in the fitted model the most recent cohorts 
#' have been zero weighted. The simulated cohorts can be seen in 
#' \code{gc.s$cohorts}. 
#'
#' @seealso \code{\link{forecast.fitStMoMo}}
#'
#'@examples
#' #Lee-Carter
#' LCfit <- fit(lc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'              ages = EWMaleData$ages, years = EWMaleData$years,
#'              ages.fit = 55:89)
#' LCsim <- simulate(LCfit, nsim = 100)
#' 
#' par(mfrow=c(1, 2))
#' plot(LCfit$years, LCfit$kt[1, ], xlim = range(LCfit$years, LCsim$kt.s$years),
#'      ylim = range(LCfit$kt, LCsim$kt.s$sim), type = "l", 
#'      xlab = "year", ylab = "kt", 
#'      main = "Lee-Carter: Simulated paths of the period index kt")
#' matlines(LCsim$kt.s$years, LCsim$kt.s$sim[1, , ], type = "l", lty = 1)
#' 
#' plot(LCfit$years, (LCfit$Dxt / LCfit$Ext)["65", ], 
#'      xlim = range(LCfit$years, LCsim$years),
#'      ylim = range((LCfit$Dxt / LCfit$Ext)["65", ], LCsim$rates["65", , ]), 
#'      type = "l", xlab = "year", ylab = "rate", 
#'      main = "Lee-Carter: Simulated mortality rates at age 65")
#' matlines(LCsim$years, LCsim$rates["65", , ], type = "l", lty = 1)
#' 
#' #APC
#' par(mfrow=c(1, 3))
#' wxt <- genWeightMat(55:89,  EWMaleData$years, clip = 3)
#' APCfit <- fit(apc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years, 
#'               ages.fit = 55:89, wxt = wxt)
#' APCsim <- simulate(APCfit, nsim = 100, gc.order = c(1, 1, 0))
#' 
#' plot(APCfit$years, APCfit$kt[1, ], 
#'      xlim = range(APCfit$years, APCsim$kt.s$years),
#'      ylim = range(APCfit$kt, APCsim$kt.s$sim), type = "l",
#'      xlab = "year", ylab = "kt",
#'      main = "APC: Simulated paths of the period index kt")
#' matlines(APCsim$kt.s$years, APCsim$kt.s$sim[1, , ], type = "l", lty = 1)
#' 
#' plot(APCfit$cohorts, APCfit$gc, 
#'      xlim = range(APCfit$cohorts, APCsim$gc.s$cohorts),
#'      ylim = range(APCfit$gc, APCsim$gc.s$sim, na.rm = TRUE), type = "l",
#'      xlab = "year", ylab = "kt", 
#'      main = "APC: Simulated paths of the cohort index (ARIMA(1,1,0))")
#' matlines(APCsim$gc.s$cohorts, APCsim$gc.s$sim, type = "l", lty = 1)
#' 
#' plot(APCfit$years, (APCfit$Dxt / APCfit$Ext)["65", ], 
#'      xlim = range(APCfit$years, APCsim$years),
#'      ylim = range((APCfit$Dxt/APCfit$Ext)["65", ], APCsim$rates["65", , ]), 
#'      type = "l", xlab = "year", ylab = "rate", 
#'      main = "APC: Simulated of mortality rates at age 65")
#' matlines(APCsim$years, APCsim$rates["65", , ], type = "l", lty = 1)
#' 
#' #Compare LC and APC
#' library(fanplot)
#' par(mfrow=c(1, 1))
#' plot(LCfit$years, (LCfit$Dxt / LCfit$Ext)["65", ], 
#'      xlim = range(LCfit$years, LCsim$years),
#'      ylim = range((LCfit$Dxt / LCfit$Ext)["65", ], LCsim$rates["65", , ], 
#'      APCsim$rates["65", , ]), type = "l", xlab = "year", ylab = "rate", 
#'      main = "Fan chart of mortality rates at age 65 (LC vs. APC)")
#' fan(t(LCsim$rates["65", , ]), start = LCsim$years[1], 
#'     probs = c(2.5, 10, 25, 50, 75, 90, 97.5), n.fan = 4,
#'     fan.col = colorRampPalette(c(rgb(1, 0, 0), rgb(1, 1, 1))), ln = NULL)
#' fan(t(APCsim$rates["65", 1:(length(APCsim$years) - 3), ]), 
#'     start = APCsim$years[1], probs = c(2.5, 10, 25, 50, 75, 90, 97.5), 
#'     n.fan = 4, fan.col = colorRampPalette(c(rgb(0, 0, 1), rgb(1, 1, 1))), 
#'     ln = NULL) 
#'@export
simulate.fitStMoMo <- function(object, nsim = 1000, seed = NULL, h = 50, 
                               oxt = NULL, gc.order = c(1, 1, 0), 
                               gc.include.constant = TRUE,
                               jumpchoice = c("fit", "actual"),
                               kt.lookback = NULL, gc.lookback = NULL, 
                               ...) {
  
  jumpchoice <- match.arg(jumpchoice)
  #Handle generato seed
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
  
  #Fit model to kt
  N <- object$model$N
  kt <- object$kt
  years <- object$years
  nYears <- length(years)
  if (is.null(kt.lookback)) 
    kt.lookback <- nYears 
  if (kt.lookback <= 0)
    stop("kt.lookback must be positive")
  kt.lookback <- min(c(kt.lookback, nYears))
  yearsSim <- (years[nYears] + 1):(years[nYears] + h)
  ages <- object$ages
  nAges <- length(ages)
  kt.h <- kt
  kt.path <- NULL
  kt.model <- NULL
  years.h <- years
  years.s <- yearsSim
  if (N > 0) {    
    kt.nNA <- max(which(!is.na(kt[1, ])))
    kt.hNA <- nYears - kt.nNA
    kt.model <- mrwd(kt[, (1 + nYears - kt.lookback):kt.nNA]) 
    if (kt.hNA > 0) {
      years.h <- years[-((kt.nNA+1):nYears)]
      years.s <- c(years[(kt.nNA+1):nYears], years.s)
      kt.h <- array(kt.h[, 1:kt.nNA], c(nrow(kt), kt.nNA))
      dimnames(kt.h)[[2]] <- years.h      
    } 
    kt.sim <- array(NA,c(N, length(years.s), nsim), list(1:N, years.s, 1:nsim))
  }
    
  #fit model to gc
  gc <- object$gc
  cohorts <- object$cohorts
  nCohorts <- length(cohorts)
  if (is.null(gc.lookback)) gc.lookback <- nCohorts 
  if (gc.lookback <= 0)
    stop("gc.lookback must be positive")
  gc.lookback <- min(c(gc.lookback, nCohorts))
  gc.h <- gc
  cohorts.h <- cohorts
  gc.model <- NULL
  gc.path <- NULL
  cohorts.s <- (cohorts[nCohorts] + 1):(cohorts[nCohorts] + h)
  if (!is.null(object$model$cohortAgeFun)) {
    gc.nNA <- max(which(!is.na(gc)))
    gc.hNA <- nCohorts - gc.nNA
    gc.model <- forecast::Arima(gc[(1 + nCohorts - gc.lookback):gc.nNA],
                                order = gc.order, 
                                include.constant = gc.include.constant) 
    if (gc.hNA > 0) {      
      gc.h <- gc[-((gc.nNA+1):nCohorts)]
      cohorts.h <- cohorts[-((gc.nNA + 1):nCohorts)]
      cohorts.s <- c(cohorts[(gc.nNA + 1):nCohorts], cohorts.s)
    }    
    gc.sim <- array(NA, c(length(cohorts.s), nsim), list(cohorts.s, 1:nsim))
  }
  #Offset
  if (is.null(oxt)) 
    oxt <- 0
  oxt.s <- matrix(oxt, nrow = nAges, ncol = h)
  colnames(oxt.s) <- yearsSim
  rownames(oxt.s) <- ages
  #Do simulations  
  forcastRates <- array(NA, c(nAges, h, nsim), 
                        list(ages, yearsSim, 1:nsim))
  fittedRates  <- array(NA, c(nAges, nYears, nsim), 
                        list(ages, years, 1:nsim))
  
  if (jumpchoice == "actual") {
    jumpoffRates <- (object$Dxt / object$Ext)[, nYears]   
  }
  
  for (i in 1:nsim) {
    if (!is.null(kt.model)) {
      kt.path <- simulate(kt.model,h + kt.hNA)
      kt.sim[, , i] <- kt.path
    }
    if (!is.null(gc.model)) {
      gc.path <- as.vector(simulate(gc.model, h + gc.hNA))
      gc.sim[, i] <- gc.path
    }    
    ratesi <- predict(object, years = c(years.h, years.s),
                     kt = cbind(kt.h,kt.path), gc = c(gc.h, gc.path),
                     oxt = cbind(object$oxt,oxt.s), type = "rates")
    
    forcastRates[, , i] <- ratesi[, (nYears + 1):(nYears + h)]
    fittedRates[, , i] <- ratesi[, 1:nYears]
    if (jumpchoice == "actual") {     
      forcastRates[, , i] <- forcastRates[, , i] * jumpoffRates / ratesi[, nYears]
    }    
  }
  
  if (is.null(kt.model)) {
    kt.s <- NULL
  } else {
    kt.s <- list(sim = kt.sim, model = kt.model, years = years.s)
  }
  if (is.null(gc.model)) {
    gc.s <- NULL
  } else {
    gc.s <- list(sim = gc.sim, model = gc.model, cohorts = cohorts.s)
  }
    structure(list(rates = forcastRates, ages = ages, 
                 years = yearsSim, kt.s = kt.s, gc.s = gc.s, oxt.s = oxt.s, 
                 fitted = fittedRates, jumpchoice = jumpchoice,
                 model = object, call = match.call()), 
            class ="simStMoMo")  
  
}


#' Print an object of class \code{"simStMoMo"}
#' 
#' \code{print} method for class \code{"simStMoMo"}. 
#' @usage 
#' \method{print}{simStMoMo}(x, ...)
#' @param x an object of class \code{"simStMoMo"}.
#' @param ... arguments to be passed to or from other methods.
#' @export 
#' @method print simStMoMo
print.simStMoMo <- function(x, ...) {
  cat("Simulations of Stochastic Mortality Model")
  cat(paste("\nCall:", deparse(x$call)))
  cat("\n\nSimulation based on")
  
  cat(paste("\nCall:", deparse(x$model$call)))  
  cat(paste("\n\nJump-off method:", x$jumpchoice))
  cat(paste("\nYears in simulation:", min(x$years), "-", max(x$years)))
  cat(paste("\nAges in simulation:", min(x$ages), "-", max(x$ages), "\n"))  
  cat(paste("\nNumber of paths:", dim(x$rates)[3], "\n"))  
}

  
  
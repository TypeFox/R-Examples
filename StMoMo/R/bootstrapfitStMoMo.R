
#' Generic method for bootstrapping a fitted Stochastic Mortality Model
#'
#' \code{bootstrap} is a generic function for bootstrapping Stochastic 
#' Mortality Models. The function invokes particular methods which depend on 
#' the class of the first argument. 
#'
#' \code{bootstrap} is a generic function which means that new fitting 
#' strategies can be added for particular stochastic mortality models. See 
#' for instance \code{\link{bootstrap.fitStMoMo}}.
#'
#' @param object an object used to select a method. Typically of class 
#' \code{fitStMoMo} or an extension of this class.
#' 
#' @param nBoot number of bootstrap samples to produce.
#'
#' @param ... arguments to be passed to or from other methods.
#'
#' @export
bootstrap =  function(object, nBoot, ...)
  UseMethod("bootstrap")


#' Bootstrap a fitted Stochastic Mortality Model
#' 
#' Produce bootstrap parameters of a Stochastic Mortality Model to account for 
#' parameter uncertainty.
#' 
#' @param object an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' @inheritParams bootstrap
#' @param type type of bootstrapping approach to be applied. 
#' \code{"semiparametric"}(default) uses the assumed distribution of the 
#' deaths to generate bootstrap samples. \code{"residual"} resamples the 
#' deviance residuals of the model to generate bootstrap samples. 
#' @param deathType type of deaths to sample in the semiparametric bootstrap. 
#' \code{"observed"} (default) resamples the observed deaths. \code{"fitted"} 
#' resamples the fitted deaths. This parameter is only used if \code{type} is 
#' \code{"semiparametric"}.
#' 
#' @return A list with class \code{"bootStMoMo"} with components:
#' 
#' \item{bootParameters}{ a list of of length \code{nBoot} with the fitted 
#' parameters for each bootstrap replication.} 
#' \item{model}{ the model fit that has been bootstrapped.}
#' \item{type}{ type of bootstrapping approach applied.}
#' \item{deathType}{ type of deaths sampled in case of semiparametric 
#' bootstrap.}
#' 
#' @details
#' When \code{type} is \code{"residual"} the residual bootstrapping approach 
#' described in Renshaw and Haberman (2008) is applied, which is an 
#' adaptation of the approach of Koissi et al (2006). In the case of a 
#' \code{"logit"} link with Binomial responses the adaptation described in 
#' Debon et al, (2010, section 3) is used.
#' 
#' When \code{type} is \code{"semiparametric"} the semiparametric approach 
#' described in Brouhns et al.(2005) is used. In the case of a \code{"logit"} 
#' link with Binomial responses a suitable adaptation is applied. If 
#' \code{deathType} is \code{"observed"} then the observed deaths are used in 
#' the sampling as in Brouhns et al. (2005) while if \code{deathType} is 
#' \code{"fitted"} the fitted deaths are used in the sampling as in 
#' Renshaw and Haberman (2008).
#' 
#' @seealso \code{\link{simulate.bootStMoMo}}, \code{\link{plot.bootStMoMo}}
#' 
#' @references 
#' Brouhns, N., Denuit M., & Van Keilegom, I. (2005). Bootstrapping the 
#' Poisson log-bilinear model for mortality forecasting. 
#' Scandinavian Actuarial Journal, 2005(3), 212-224. 
#' 
#' Debon, A., Martinez-Ruiz, F., & Montes, F. (2010). A geostatistical 
#' approach for dynamic life tables: The effect of mortality on remaining 
#' lifetime and annuities. Insurance: Mathematics and Economics, 47(3), 
#' 327-336. 
#' 
#' Renshaw, A. E., & Haberman, S. (2008). On simulation-based approaches to 
#' risk measurement in mortality with specific reference to Poisson Lee-Carter 
#' modelling. Insurance: Mathematics and Economics, 42(2), 797-816. 
#' 
#' @examples 
#' #Long computing times
#' \dontrun{
#' LCfit <- fit(lc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'              ages = EWMaleData$ages, years = EWMaleData$years)
#' 
#' LCResBoot <- bootstrap(LCfit, nBoot = 500, type = "residual")
#' plot(LCResBoot)
#' LCSemiObsBoot <- bootstrap(LCfit, nBoot = 500, type = "semiparametric")
#' plot(LCSemiObsBoot)
#' LCSemiFitBoot <- bootstrap(LCfit, nBoot = 500, type = "semiparametric",  
#'                            deathType = "fitted")
#' plot(LCSemiFitBoot)
#' }
#' 
#' @export
bootstrap.fitStMoMo <- function(object, nBoot = 1, 
                                type = c("semiparametric", "residual"),
                                deathType = c("observed", "fitted"), ...) {
  
  type <- match.arg(type)  
  deathType <- match.arg(deathType)
  link <- object$model$link
  
  #Generate death samples
  
  if (link == "log") {
    if (type == "residual") {
      bootSamples <- genPoissonResBootSamples(
        devRes = residuals(object, scale = FALSE), 
        dhat = fitted(object, type = "deaths"), nBoot)
    } else if (type == "semiparametric") {     
        if (deathType == "observed") {
          D <- object$Dxt
        } else { 
          D <- fitted(object, type = "deaths")
        }
      bootSamples <- genPoissonSemiparametricBootSamples(D = D, nBoot)
      
    }
  } else if (link == "logit") {
    if (type == "residual") {
      bootSamples <- genBinomialResBootSamples(
        devRes = residuals(object, scale = FALSE),
        dhat = fitted(object, type = "deaths"), E = object$Ext, nBoot)
    } else if (type == "semiparametric") {
      if (deathType == "observed") {
        D <- object$Dxt
      } else { 
        D <- fitted(object, type = "deaths")
      }
      bootSamples <- genBinomialSemiparametricBootSamples(D = D, 
                                                          E = object$Ext,
                                                          nBoot)      
    }
  }

  #Fit the model to each of 
  refit <- function(Dxt) {
    getMinimalFitStMoMo(fit(object = object$model, Dxt = Dxt, Ext = object$Ext, 
                            ages = object$ages, years = object$years,
                            oxt = object$oxt, wxt = object$wxt, 
                            start.ax = object$ax, start.bx = object$bx, 
                            start.kt = object$kt, start.b0x = object$b0x, 
                            start.gc = object$gc, verbose = FALSE))
  }
  bootParameters <- lapply(bootSamples, FUN = refit)      
  structure(list(bootParameters = bootParameters, model = object, type = type, 
                 deathType = deathType, call = match.call()), 
            class = "bootStMoMo")  
}

#' Print an object of class \code{"bootStMoMo"}
#' 
#' \code{print} method for class \code{"bootStMoMo"}. 
#' @usage 
#' \method{print}{bootStMoMo}(x, ...)
#' @param x an object of class \code{"bootStMoMo"}.
#' @param ... arguments to be passed to or from other methods.
#' @export 
#' @method print bootStMoMo
 print.bootStMoMo <- function(x,...) {
  cat("Bootstrapped Stochastic Mortality Model")
  cat(paste("\nCall:", deparse(x$call)))
  cat("\n\n")  
  if (x$type == "residual") {
    cat("Residual bootstrap based on\n")    
  } else {
    cat(paste("Semiparametric bootstrap using", x$deathType, 
              "deaths and based on\n"))        
  }
  print(x$model$model)  
  cat(paste("\n\nNumber of bootstrap samples:", length(x$bootParameters), "\n"))  
}


#' Map poisson deviance residuals into deaths
#' 
#' This funciont uses the procedure described in Renshaw and Haberman (2008)
#' 
#' @param dhat fitted number of deaths (scalar)
#' @param res deviance residual (scalar)
#' 
#' @return matching observed deaths
#' 
#' @references
#' Renshaw, A. E., & Haberman, S. (2008). On simulation-based approaches to 
#' risk measurement in mortality with specific reference to Poisson Lee-Carter 
#' modelling. Insurance: Mathematics and Economics, 42(2), 797-816. 
#' 
#' @keywords internal
poissonRes2death <- function(dhat, res) {
  
  if (is.na(dhat)) 
    return(NA)
  a <- log(dhat)
  c <- res^2 / 2 - dhat
  g <- function(d){
    d * log(d) - d * (1 + a) - c  
  }
  dg <- function(d){
    log(d) - a  
  }
  if ( (res < 0 & c > 0) | (res == -sqrt(2 * dhat) & c == 0)){
    d <- 0
  } else {
    start <- max(0.01, dhat + sign(res) * 0.5 * dhat)
    d <- rootSolve::multiroot(f = g, start = start, 
                              jacfunc = dg, positive = TRUE)$root    
  }
  d  
}

#' Generate poisson residual bootstrap sample
#' 
#' Generate poisson bootstrap samples using the procedure
#' described in Renshaw and Haberman (2008)
#' 
#' @param devRes Matrix of deviance residuals
#' @param dhat Marix of fitted deaths
#' @inheritParams bootstrap
#' 
#' @return a list of length \code{nBoot} of marices with matching sampled 
#' deaths
#' 
#' @references
#' Renshaw, A. E., & Haberman, S. (2008). On simulation-based approaches to 
#' risk measurement in mortality with specific reference to Poisson Lee-Carter 
#' modelling. Insurance: Mathematics and Economics, 42(2), 797-816. 
#' 
#' @keywords internal
genPoissonResBootSamples <- function(devRes, dhat, nBoot) {    
  nYears <- ncol(dhat)
  nAges <- nrow(dhat)
  nObs <- nYears * nAges
  sampleRes <- as.vector(devRes$residuals)
  sampleRes <- sampleRes[!is.na(sampleRes)]
  vecdHat <- as.vector(dhat)
  
  funSampleDeath <- function(x = 0) {
    resStar <- sample(sampleRes, nObs, replace = TRUE)
    dStar <- mapply(poissonRes2death, vecdHat, resStar)    
  }  
  bootSamples <- lapply(X = 1:nBoot,FUN = funSampleDeath)
  lapply(X = bootSamples, FUN = matrix, nrow = nAges, ncol = nYears)  
}

#' Map Binomial deviance residuals into deaths
#' 
#' This funciont uses the procedure described in 
#' Debon et al. (2010, section 3)
#' 
#' @param dhat fitted number of deaths (scalar)
#' @param e number of exposures (scalar)
#' @param res deviance residual (scalar)
#' 
#' @return matching observed deaths
#' 
#' @references
#' Debon, A., Martinez-Ruiz, F., & Montes, F. (2010). A geostatistical 
#' approach for dynamic life tables: The effect of mortality on remaining 
#' lifetime and annuities. Insurance: Mathematics and Economics, 
#' 47(3), 327-336. 
#' 
#' @keywords internal
binomRes2q <- function(dhat, e, res) {
  
  if (is.na(dhat) || is.na(e)) 
    return(NA)
  if (e <= 0) return(0L)
  a <- log(dhat / (e - dhat))
  c <- res^2 / 2 + e * log(e - dhat)
  g <- function(d) {
    d * log(d / (e - d)) + e * log(e - d) - a * d - c  
  }
  dg <- function(d) {
    log(d / (e - d)) - a  
  }
  
  if (c == 0) {
    d <- 0
  } else {
    start <- min(max(0.01, dhat + sign(res) * 0.5 * dhat), e - 0.01)
    d <- rootSolve::multiroot(f = g, start = start, jacfunc = dg, 
                              positive = TRUE)$root    
  }
  min(d, e)
}

#' Generate Binomial residual bootstrap sample
#' 
#' Generate Binomial bootstrap samples using the procedure
#' described in Debon et al. (2010, section 3)
#' 
#' @param devRes Matrix of deviance residuals
#' @param dhat Matrix of fitted deaths
#' @param E Matrix of observed exposures
#' @inheritParams bootstrap
#' 
#' @return a list of length \code{nBoot} of matrices with matching sampled 
#' deaths
#' 
#' @references
#' Debon, A., Martinez-Ruiz, F., & Montes, F. (2010). A geostatistical 
#' approach for dynamic life tables: The effect of mortality on remaining 
#' lifetime and annuities. Insurance: Mathematics and Economics, 47(3), 
#' 327-336.  
#' 
#' @keywords internal
genBinomialResBootSamples <- function(devRes, dhat, E, nBoot) {
  
  nYears <- ncol(dhat)
  nAges <- nrow(dhat)
  nObs <- nYears * nAges
  sampleRes <- as.vector(devRes$residuals)
  sampleRes <- sampleRes[!is.na(sampleRes)]
  vecdHat <- as.vector(dhat)
  vecE <- as.vector(E)
  
  funSampleDeath <- function(x = 0) {
    resStar <- sample(sampleRes, nObs, replace = TRUE)
    dStar <- mapply(FUN = binomRes2q, vecdHat, vecE, resStar)    
  }
  
  bootSamples <- lapply(X = 1:nBoot, FUN = funSampleDeath)
  lapply(X = bootSamples, FUN = matrix, nrow = nAges, ncol = nYears)
  
}

#' Generate Poisson semiparametric bootstrap samples
#' 
#' Generate Poisson semiparametric bootstrap samples using the procedure
#' described in Brouhns et al (2005).
#' 
#' @param D matrix of deaths
#' @inheritParams bootstrap
#' 
#' @return a list of length \code{nBoot} of matrices with sampled deaths
#' 
#' @references
#' Brouhns, N., Denuit M., & Van Keilegom, I. (2005). Bootstrapping the 
#' Poisson log-bilinear model for mortality forecasting. 
#' Scandinavian Actuarial Journal, 2005(3), 212-224. 
#' 
#' @keywords internal
genPoissonSemiparametricBootSamples <- function(D, nBoot){
  nYears <- ncol(D)
  nAges <- nrow(D)
  nObs <- nYears * nAges
  vecD <- as.vector(D)  
  
  funSampleDeath <- function(x = 0) {    
    dStar <- sapply(X = vecD, FUN = function(lambda) rpois(1, lambda))    
  }  
  bootSamples <- lapply(X = 1:nBoot, FUN = funSampleDeath)
  lapply(X = bootSamples, FUN = matrix, nrow = nAges, ncol = nYears)
}

#' Generate Binomial semiparametric bootstrap samples
#' 
#' Generate Binomial semiparametric bootstrap samples using a suitable 
#' adaptation of the Poisson procedure described in Brouhns et al (2005).
#' 
#' @param D matrix of deaths
#' @param E matrix of exposures
#' @inheritParams bootstrap
#' 
#' @return a list of length \code{nBoot} of matrices with sampled deaths
#' 
#' @references
#' Brouhns, N., Denuit M., & Van Keilegom, I. (2005). Bootstrapping the 
#' Poisson log-bilinear model for mortality forecasting. 
#' Scandinavian Actuarial Journal, 2005(3), 212-224. 
#' 
#' @keywords internal
genBinomialSemiparametricBootSamples <- function(D, E, nBoot) {
  nYears <- ncol(D)
  nAges <- nrow(D)
  nObs <- nYears * nAges
  vecRate <- as.vector(D / E)  
  vecE <- round(as.vector(E))
  
  funSampleDeath <- function(x = 0) {    
    dStar <- mapply(FUN = function(e, q) rbinom(1, e, q), vecE, vecRate)    
  }  
  bootSamples <- lapply(X = 1:nBoot, FUN = funSampleDeath)
  lapply(X = bootSamples, FUN = matrix, nrow = nAges, ncol = nYears)
}
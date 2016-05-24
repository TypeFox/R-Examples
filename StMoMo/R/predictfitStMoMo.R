#' Predict method for Stochastic Mortality Models fits
#' 
#' Obtain predictions from a Stochastic Mortality Model
#' fit. 
#' 
#' This function evaluates 
#' \deqn{\hat{\eta}_{xt} = o_{xt} + \alpha_x + 
#' \sum_{i=1}^N \beta_x^{(i)}\kappa_t^{(i)} + \beta_x^{(0)}\gamma_{t-x}}
#' for a fitted Stochastic Mortality model.
#' In producing a prediction the static age function, \eqn{\alpha_x}, and the
#' age-modulating parameters, \eqn{\beta_x^{(i)}, i=0, ..., N}, are taken from
#' the fitted model in \code{object} while the period indexes, 
#' \eqn{\kappa_t^{(i)}, i=1,..., N}, and cohort index, \eqn{\gamma_{t-x}}, 
#' are taken from the function arguments.  
#' 
#' This function can be useful, for instance, in producing forecasts of 
#' mortality rates using time series models different to those available 
#' in \code{\link{forecast.fitStMoMo}} (See examples below). 
#' 
#' @param object an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' 
#' @param years vector of years for which a prediction is required.
#' 
#' @param kt matrix of values of the period indexes to use for the prediction. 
#' If the model has any age-period term this argument needs to be provided and 
#' the number of rows in \code{kt} must be equal to the number of age-period 
#' terms in the model and the number of columns in \code{kt} must correspond 
#' to the length of \code{years}. If the Stochastic Mortality Model doesn't 
#' have any age-period terms this argument is ignored and needs not be 
#' provided.
#'  
#' @param gc vector of values of the cohort indexes to use for the prediction. 
#' If the model has a cohort effect this argument needs to be provided. 
#' In this case the length of \code{gc} must be equal to the number of cohorts 
#' for which a prediction is being produced, namely, 
#' \code{length(object$ages) + length(years) - 1}. If the Stochastic Mortality 
#' Model doesn't have a cohort effect  this argument is ignored and needs not 
#' be provided. 
#' 
#' @param oxt optional matrix/vector or scalar of known offset to be used in 
#' the prediction.
#' @param type the type of the predicted values that should be returned. The 
#' alternatives are \code{"link"}(default) and \code{"rates"}.
#' @param ... other arguments.
#'   
#' @return A matrix with the predicted values.
#'   
#' @seealso \code{\link{forecast.fitStMoMo}}   
#'   
#' @examples
#' library(forecast)
#' #Lee-Carter forecast using auto.arima
#' LCfit <- fit(lc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'              ages = EWMaleData$ages, years = EWMaleData$years,
#'              ages.fit = 55:89)
#' ktForLC <- forecast(auto.arima(as.vector(LCfit$kt)), h = 30) 
#' mxtForLC <- predict(LCfit, years = 2012:2041, kt = ktForLC$mean, 
#'                     type = "rates")
#' mxthatLC <- fitted(LCfit, type = "rates")
#' mxt <- LCfit$Dxt / LCfit$Ext
#' plot(1961:2041, (cbind(mxthatLC, mxtForLC))["80", ], type = "l", 
#'      xlab = "year", ylab = "death rate", 
#'      main = "Fitted vs. Observed rates at age 80")
#' points(1961:2011, mxt["80", ])
#' 
#' #Age-Period-Cohort forecast using auto.arima
#' APCfit <- fit(apc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'               ages = EWMaleData$ages, years = EWMaleData$years,
#'               ages.fit = 55:89)
#' ktForAPC <- forecast(auto.arima(as.vector(APCfit$kt)), h = 30)
#' gcForAPC <- forecast(auto.arima(as.vector(APCfit$gc), max.d = 1), h = 30)
#' mxtForAPC <- predict(APCfit, years = 2012:2041, kt = ktForAPC$mean, 
#'                      gc = c(tail(APCfit$gc, 34), gcForAPC$mean), 
#'                      type = "rates")
#' mxthatAPC <- fitted(APCfit, type = "rates")
#' lines(1961:2041 , (cbind(mxthatAPC, mxtForAPC))["80", ], type = "l", 
#'       col = "blue")
#'
#' @export 
predict.fitStMoMo <- function(object, years, kt = NULL, gc = NULL, oxt = NULL, 
                              type = c("link", "rates"), ...) {
  type <- match.arg(type)
  ages <- object$ages
  nAges <- length(ages)
  nYears <- length(years)
  cohorts <- (years[1] - ages[nAges]):(years[nYears] - ages[1])
  nCohorts <- length(cohorts)
  
  #Check inputs
  if (object$model$N > 0) {
    kt <- as.matrix(kt)
    if (ncol(kt) == 1) 
      kt <- t(as.matrix(kt))
    if (ncol(kt) != nYears) {
      stop( "Mismatch between the dimension of kt and the number of years.")
    }
    if (nrow(kt) != object$model$N) {
      stop("Mismatch between the dimension of kt and the number of age-period 
           terms in the model.")
    }
  } else {
    if (!is.null(kt)) {
      warning("kt argument ignored as the model doesn't have any age-period
              terms.")
    }
    kt <- NULL
  }
  if (!is.null(object$model$cohortAgeFun)) {
    gc <- as.vector(gc)
    if (length(gc) != nCohorts) {
      stop(paste("Mismatch between the length of gc and the number of cohorts 
                  required. Required number of cohorts:", nCohorts))
    }
  } else {
    if (!is.null(gc)) {
      warning("gc argument ignored as the model doesn't have an age-cohort
              term.")
     }
    gc <- NULL
  }
  if (is.null(oxt)) {
    oxt <- matrix(0, nrow = nAges, ncol = nYears)
    rownames(oxt) <- ages
    colnames(oxt) <- years        
  } else {
    oxt <- matrix(oxt, nrow = nAges, ncol = nYears)
    rownames(oxt) <- ages
    colnames(oxt) <- years             
  }   
  
  #compute prediction
  link <- predictLink(ax = object$ax, bx = object$bx, kt = kt, 
                      b0x = object$b0x, gc = gc, oxt = oxt, 
                      ages = ages, years = years)
  rates <- switch(object$model$link, log = exp(link), logit = invlogit(link))
  switch(type, rates = rates, link = link)  
}
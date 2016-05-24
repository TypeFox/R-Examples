#' Compute fitted values for a Stochastic Mortality Model
#' 
#' Returns fitted values for the data used in fitting a Stochastic Mortality 
#' Model.
#' 
#' @param object an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' @param type the type of the fitted values that should be returned. The 
#' alternatives are \code{"link"}(default), \code{"rates"}, and 
#' \code{"deaths"}.
#' @param ... other arguments.
#'
#' @return A matrix with the fitted values.
#' 
#' @examples
#' LCfit <- fit(lc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'              ages = EWMaleData$ages, years = EWMaleData$years,
#'              ages.fit = 55:89)
#' matplot(LCfit$ages, fitted(LCfit), type = "l", lty = 1, 
#'         col = rainbow(length(LCfit$years)), xlab = "year", 
#'         ylab = "log death rate", main = "Fitted rates")
#' 
#' qxthat <- fitted(LCfit, type = "rates")
#' qxt <- LCfit$Dxt / LCfit$Ext
#' plot(LCfit$years, qxt["65", ], xlab = "year", ylab = "death rate",
#'      main = "fitted vs. observed rates at age 65")
#' lines(LCfit$years, qxthat["65", ])
#' @export 
fitted.fitStMoMo<-function(object, type = c("link", "rates", "deaths"), ...) {
  
  type <- match.arg(type)
  link <- with(object, predictLink(ax = ax, bx = bx, kt = kt, b0x = b0x, 
                                   gc = gc, oxt = oxt, ages = ages, 
                                   years = years))
  rates <- switch(object$model$link, log = exp(link), logit = invlogit(link))
  deaths <- object$Ext * rates
  switch(type, rates = rates, deaths = deaths, link = link)  
}

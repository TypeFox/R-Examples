if(!isGeneric("logLik"))
  setGeneric(name = "logLik", def = function(object, ...) standardGeneric("logLik"))

setMethod(f = "logLik", signature = "glimML", definition = function(object, ...){
  val <- object@logL
  attr(val, "df") <- object@nbpar
  attr(val, "nobs") <- df.residual(object) + object@nbpar
  class(val) <- "logLik"
  val
  })

if(!isGeneric("AIC"))
  setGeneric("AIC", function(object, ..., k = 2) standardGeneric("AIC"))

## Removed to comply with recommandations on avoiding to redefine generic functions from the stats package (24 Aug 2008)
##setMethod(f = "AIC", signature = "logLik", definition = function(object, ..., k = 2){
##  npar <- attr(object, "df")
##  nobs <- attr(object, "nobs")
##  c(AIC = -2 * c(object) + k * npar, AICc = -2 * c(object) + k * npar + 2 * npar * (npar + 1) / (nobs - npar - 1))
##  })

setMethod(f = "AIC", signature = "glimML", definition = function(object, ..., k = 2){
  ## local function to compute AIC and AICc
  AIC1 <- function(x, k = k){
  npar <- attr(x, "df")
  nobs <- attr(x, "nobs")
  c(AIC = -2 * as.vector(x) + k * npar, AICc = -2 * as.vector(x) + k * npar + 2 * npar * (npar + 1) / (nobs - npar - 1))
  }
  ## Actual computation
  object <- list(object, ...)
  val <- lapply(object, logLik)
  val <- as.data.frame(t(sapply(val, function(el) c(attr(el, "df"), AIC1(el, k = k)))))
  names(val) <- c("df", "AIC", "AICc")
  Call <- match.call()
  Call$k <- NULL
  row.names(val) <- as.character(Call[-1])
  new(Class = "aic", istats = val)
  })

setMethod(f = "show", signature = "aic", definition = function(object) print(object@istats))

if(!isGeneric("summary"))
  setGeneric(name = "summary", def = function(object, ...) standardGeneric("summary"))

setMethod(f = "summary", signature = "aic", definition = function(object, which = c("AIC", "AICc")){
  which <- match.arg(which)
  if(!is.element(which, c("AIC", "AICc")))
    stop(which, " is not a valid choice: must be either ", dQuote("AIC"), " or ", dQuote("AICc"), ".\n", sep = "")
  x <- object@istats
  res <- x[, match(c("df", which), table = names(x))]
  aic <- res[ , 2]
  res <- res[order(aic), ]
  aic <- sort(aic)
  res$Delta <- aic - min(aic)
  res$W <- exp(-.5 * res$Delta) / sum(exp(-.5 * res$Delta))
  res$Cum.W <- cumsum(res$W)
  res
  })

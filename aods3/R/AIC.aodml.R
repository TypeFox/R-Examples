AIC.aodml <- function(object, ..., k = 2) {
  
  # AIC = -2 x logL + ...
  
  # local function (x is a logLik object)
	zaic <- function(x, k = k){
  	df.model <- attr(x, "df")
  	nobs <- attr(x, "nobs")
    ## NB for AICc, k = 2 by definition
  	c(AIC = -2 * as.vector(x) + k * df.model,
  		AICc = -2 * as.vector(x) + 2 * df.model + 2 * df.model * (df.model + 1) / (nobs - df.model - 1))
    }
	
  # actual computation
  object <- list(object, ...)
  val <- lapply(object, logLik)
  val <- as.data.frame(t(sapply(val, function(z) c(attr(z, "nobs"), attr(z, "df"), zaic(z, k = k)))))
  names(val) <- c("nobs", "df.model", "AIC", "AICc")
  call <- match.call()
  row.names(val) <- as.character(call[-1])
	val
	
}

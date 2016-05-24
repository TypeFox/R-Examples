waldci <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = NULL, ...)
{
  UseMethod("waldci")
}

waldci.default <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = NULL, ...)
{
  ## use S4 methods if loaded
  coef0 <- if("stats4" %in% loadedNamespaces()) stats4::coef else coef
  vcov0 <- if("stats4" %in% loadedNamespaces()) stats4::vcov else vcov

  ## extract coefficients and standard errors
  est <- coef0(x)
  if(is.null(vcov.)) se <- vcov0(x) else {
      if(is.function(vcov.)) se <- vcov.(x)
        else se <- vcov.
  }
  se <- sqrt(diag(se))

  ## match using names and compute t/z statistics
  if(!is.null(names(est)) && !is.null(names(se))) {
    anames <- names(est)[names(est) %in% names(se)]
    est <- est[anames]
    se <- se[anames]
  }
  
  ## process level
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  
  ## get quantile from central limit theorem
  if(is.null(df)) {
    df <- try(df.residual(x), silent = TRUE)
    if(inherits(df, "try-error")) df <- NULL
  }
  if(is.null(df)) df <- 0
  fac <- if(is.finite(df) && df > 0) qt(a, df = df) else qnorm(a)

  ## set up confidence intervals
  ci <- cbind(est + fac[1] * se, est + fac[2] * se)
  colnames(ci) <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3L), "%")
  
  ## process parm
  if(is.null(parm)) parm <- seq_along(est)
  if(is.character(parm)) parm <- which(names(est) %in% parm)
  ci <- ci[parm, , drop = FALSE]
  return(ci)
} 

waldci.glm <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = Inf, ...)
  waldci.default(x, parm = parm, level = level, vcov. = vcov., df = df, ...)  

waldci.mlm <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = NULL, ...)
{
  ## obtain vcov
  v <- if(is.null(vcov.)) vcov(x) else if(is.function(vcov.)) vcov.(x) else vcov.

  ## nasty hack: replace coefficients so that their names match the vcov() method
  x$coefficients <- structure(as.vector(x$coefficients), .Names = colnames(vcov(x)))

  ## call default method
  waldci.default(x, parm = parm, level = level, vcov. = v, df = df, ...)
}

waldci.survreg <- function(x, parm = NULL, level = 0.95, vcov. = NULL, df = Inf, ...)
{
  if(is.null(vcov.)) v <- vcov(x) else {
      if(is.function(vcov.)) v <- vcov.(x)
  	else v <- vcov.
  }
  if(length(x$coefficients) < NROW(x$var)) {
    x$coefficients <- c(x$coefficients, "Log(scale)" = log(x$scale))
  }
  waldci.default(x, parm = parm, level = level, vcov. = v, df = df, ...)  
} 

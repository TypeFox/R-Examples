ArchTest <- function (x, lags=12, demean = FALSE) 
{
# Capture name of x for documentation in the output  
  xName <- deparse(substitute(x))
# 
  x <- as.vector(x)
  if(demean) x <- scale(x, center = TRUE, scale = FALSE)
#  
  lags <- lags + 1
  mat <- embed(x^2, lags)
  arch.lm <- summary(lm(mat[, 1] ~ mat[, -1]))
  STATISTIC <- arch.lm$r.squared * length(resid(arch.lm))
  names(STATISTIC) <- "Chi-squared"
  PARAMETER <- lags - 1
  names(PARAMETER) <- "df"
  PVAL <- 1 - pchisq(STATISTIC, df = PARAMETER)
  METHOD <- "ARCH LM-test;  Null hypothesis:  no ARCH effects"
  result <- list(statistic = STATISTIC, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name =
                 xName)
  class(result) <- "htest"
  return(result)
}

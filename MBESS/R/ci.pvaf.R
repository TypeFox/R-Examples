ci.pvaf <- function (F.value = NULL, df.1 = NULL, df.2 = NULL, N = NULL, 
          conf.level = 0.95, alpha.lower = NULL, alpha.upper = NULL, 
          ...) 
{
  if (is.null(alpha.lower) & is.null(alpha.upper)) {
    alpha.lower <- (1 - conf.level)/2
    alpha.upper <- (1 - conf.level)/2
  }
  if (!is.null(alpha.lower) && !is.null(alpha.upper)) {
    conf.level <- 1 - alpha.upper - alpha.lower
  }
  if (!is.null(alpha.lower) & is.null(alpha.upper)) 
    stop("This is a problem with the desired confidence level ('alpha.lower' specified but not 'alpha.upper').")
  if (is.null(alpha.lower) & !is.null(alpha.upper)) 
    stop("This is a problem with the desired confidence level ('alpha.upper' specified but not 'alpha.lower').")
  if (is.null(df.1) | is.null(df.2) | is.null(N)) 
    stop("You need to specify 'df.1', 'df.2', and 'N'.")
  if (is.null(F.value)) 
    stop("You must specify the observed F-value ('F.value') from the analysis of variance.")
  if (alpha.lower > 0.5 || alpha.lower < 0) 
    stop(" 'alpha.lower' must be smaller than .5 and nonnegative.")
  if (conf.level > 1 || conf.level < 0) 
    stop(" 'conf.level' must be larger than 0 and smaller than 1. ")
  if (F.value <= 0) 
    stop(" 'F.value' must be larger than 0. ")
  if (N <= 0 || N <= df.1 + df.2) 
    stop("N must be larger than df.1+df.2")
  
  Lims <- conf.limits.ncf(F.value = F.value, conf.level = NULL, 
                          df.1 = df.1, df.2 = df.2, alpha.lower = alpha.lower, 
                          alpha.upper = alpha.upper, ...)
  
  # The following checks are when the lower limit is theoreteically under zero; it is set to zero.
  if(is.na(Lims$Lower.Limit)) {Lims$Lower.Limit <- 0}
  if(is.na(Lims$Prob.Less.Lower)) {Lims$Prob.Less.Lower <- 0}
  
  Actual.Coverage <- 1 - (Lims$Prob.Greater.Upper + Lims$Prob.Less.Lower)
  
  Lower.lim <- Lims$Lower.Limit/(Lims$Lower.Limit + N)
  Upper.lim <- Lims$Upper.Limit/(Lims$Upper.Limit + N)
  
  print(paste("The", 1 - (alpha.lower + alpha.upper), "confidence limits (and the actual confidence interval coverage) for the proportion of variance of the dependent variable accounted for by knowing group status are given as:"))
  
  return(list(Lower.Limit.Proportion.of.Variance.Accounted.for = Lower.lim,
              Probability.Less.Lower.Limit = Lims$Prob.Less.Lower,
              Upper.Limit.Proportion.of.Variance.Accounted.for = Upper.lim,
              Probability.Greater.Upper.Limit = Lims$Prob.Greater.Upper,
              Actual.Coverage=Actual.Coverage))
}
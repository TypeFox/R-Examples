##  recog.R
##
## $Revision: 1.4 $ $Date: 2014/06/24 02:13:35 $
##

recogniseCdf <- function(s="punif") {
  if(!is.character(s) || length(s) != 1) return(NULL)
  if(nchar(s) <= 1 || substr(s,1,1) != "p") return(NULL)
  root <- substr(s, 2, nchar(s))
  a <- switch(root,
              beta     = "beta",
              binom    = "binomial",
              birthday = "birthday coincidence",
              cauchy   = "Cauchy",
              chisq    = "chi-squared",
              exp      = "exponential",
              f        = "F",
              gamma    = "Gamma",
              geom     = "geometric",
              hyper    = "hypergeometric",
              lnorm    = "log-normal",
              logis    = "logistic",
              nbinom   = "negative binomial",
              norm     = "Normal",
              pois     = "Poisson",
              t        = "Student's t",
              tukey    = "Tukey (Studentized range)",
              unif     = "uniform",
              weibull  = "Weibull",
              NULL)
  if(!is.null(a))
    return(paste(a, "distribution"))
  b <- switch(root,
              AD     = "Anderson-Darling",
              CvM    = "Cramer-von Mises",
              wilcox = "Wilcoxon Rank Sum",
              NULL)
  if(!is.null(b))
    return(paste("null distribution of", b, "Test Statistic"))
  return(NULL)
}

         

Unitroot <- function(x, trend=c("c", "nc", "ct"), method=c('adf', 'McKinnon'),
                     lags= 2 ){
  type <- match.arg(trend)
  mtd <- match.arg(method)
  lag <- max(1, lags-1)
#
  {
    if(require('fUnitRoots'))
      urtest <- {
        if(mtd=="adf")
          adfTest(x, lags=lag, type=type)
        else
          unitrootTest(x, lags=lag, type=type)
      }
    else
      stop('require(fUnitRoots);  not available.')
  }
  urtest@call <- match.call()
#  urtest@data.name <- deparse(substitute(x))
#
  urtest
}

summary.fHTEST <- function(object, ...){
  if(!require('fUnitRoots'))stop('Need package fUnitRoots')
#
  test <- object@test
#
  stat <- test$statistic
  type <- names(stat)
  cat("Test for Unit Root:",
      c("McKinnon's", "Augmented DF")[1+(type=="Dickey-Fuller")],
      "test\n")
  cat("Null Hypothesis:  There is a unit root.\n")
  cat("   Type of Test:  t test\n")
#
  cat(" Test Statistic:  ", round(stat, 3), "\n")
#
  Call <- object@call
#
  p.value <- padf(stat, trend=Call[["trend"]], statistic="t")
  cat("        P-value: ", round(p.value, 4), "\n")
  cat("NOTE:  p.value by linear interpolation in a table;\n")
  cat("NOTE:  In the example on Tsay, p. 70, it differed\n")
  cat("NOTE:  from the S-PLUS Finmetrics answer by 10%\n")
#
  print(summary(test$lm))
  cat("NOTE:  The order of terms and labeling here are different\n")
  cat("NOTE:  from S-PLUS Finmetics but seem more consistent\n")
  cat("NOTE:  with the ADF definition per Tsay, expression (2.40).\n")
  cat("NOTE:  S-PLUS Finmetrics reports 7 more degrees of freedom\n")
  cat("NOTE:  than adfTest{fUnitRoots} in the example on Tsay, p. 70.\n")
#
  invisible(object)
}

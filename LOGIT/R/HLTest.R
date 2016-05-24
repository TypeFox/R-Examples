#HLTest.r  Bilder & Loughin, Analysis of Categorical Data with R, 2015, Chapman & HAll/CRC
#' @title  Display Hosmer-Lemeshow statistic and table of probabilities following logistic
#' regression using glm with binomial family.
#' @description Provides a Hosmer-Lemeshow statistic and table following logistic regression.
#' @aliases HLTest
#' @usage HLTest(obj, g)
#' @importFrom stats family quantile xtabs pchisq
#' @format \describe{
#' \item{x}{
#' The function has two arguments: model name, number of groups}
#' }
#' @details HLTest is a post-estimation function for logistic regression, following the use
#' of glm(). Usage displays a table of observed vs predicted groups and an overall
#' H-L goodness-of-fit statistic.
#' @param obj model name
#' @param g number of groups
#' @return list
#' @note HLTest must be loaded into memory in order to be effectve. As a function in LOGIT,
#' it is immediately available to a user. My thanks to Bilger and Loughlin for the
#' use of their function.
#'@examples
#' library(MASS)
#' library(LOGIT)
#' data(medpar)
#' mylogit <- glm( died ~  los + white + hmo, family=binomial, data=medpar)
#' grp10 <- HLTest(obj=mylogit, g=10)
#' cbind(grp10$observed, round(grp10$expect, digits = 1))
#' grp10
#'
#' @seealso \code{\link{glm}}
#' @author Adapted from Loughlin, T.M. in Bilder and Loughlin, 2015
#'
#' @references Hilbe, J. M. (2015), Practical Guide to Logistic Regression, Chapman & Hall/CRC.
#'
#' Bilder, C.R. and Loughlin, T.M. (2015), Analysis of Categorical Data with R, Chapman & Hall/CRC.
#'
#' Hilbe, J. M. (2009), Logistic Regression Models, Chapman & Hall/CRC.
#'
#' Hosmer, D.W., Lemeshow, S, and Sturdivant, R.X (2013), Applied Logistic Regression, 3rd ed, Wiley.
#' @keywords models
#' @export
#'
HLTest<-function(obj, g) {
 # first, check to see if we fed in the right kind of object
 stopifnot(family(obj)$family == "binomial" && family(obj)$link == "logit")
 y = obj$model[[1]]
 trials = rep(1, times = nrow(obj$model))
 if(any(colnames(obj$model) == "(weights)"))
  trials <- obj$model[[ncol(obj$model)]]
 # the double bracket (above) gets the index of items within an object
 if (is.factor(y))
  y = as.numeric(y) == 2  # Converts 1-2 factor levels to logical 0/1 values
 yhat = obj$fitted.values
 interval = cut(yhat, quantile(yhat, 0:g/g), include.lowest = TRUE)  # Creates factor with levels 1,2,...,g
 Y1 <- trials*y
 Y0 <- trials - Y1
 Y1hat <- trials*yhat
 Y0hat <- trials - Y1hat
 obs = xtabs(formula = cbind(Y0, Y1) ~ interval)
 expect = xtabs(formula = cbind(Y0hat, Y1hat) ~ interval)
 if (any(expect < 5))
  warning("Some expected counts are less than 5. Use smaller number of groups")
 pear <- (obs - expect)/sqrt(expect)
 chisq = sum(pear^2)
 P = 1 - pchisq(chisq, g - 2)
 # by returning an object of class "htest", the function will perform like the
 # built-in hypothesis tests
 return(structure(list(
  method = c(paste("Hosmer and Lemeshow goodness-of-fit test with", g, "bins", sep = " ")),
  data.name = deparse(substitute(obj)),
  statistic = c(X2 = chisq),
  parameter = c(df = g-2),
  p.value = P,
  pear.resid = pear,
  expect = expect,
  observed = obs
 ), class = 'htest'))
}

###
### mombf.R
###

mombf <- function(lm1,coef,g,prior.mode,baseDensity='normal',nu=3,theta0,logbf=FALSE,B=10^5) {
  if (!(baseDensity %in% c('normal','t'))) stop('baseDensity must either be normal or t')
  if (baseDensity=='t') {
    if (length(g)>1) stop("For baseDensity=='t' g is only allowed to have length 1")
    if (nu <= 2) stop('For t base density nu must be >=3, otherwise prior is improper')
    if (nu >= length(lm1$residuals)-length(coef)+2) warning(paste('Based on the amount of data available, nu should be <',length(lm1$residuals)-length(coef)+2,'in order to guarantee finite sample consistency'))
  }
  UseMethod("mombf")
}


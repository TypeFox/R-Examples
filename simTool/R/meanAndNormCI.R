#' A convenient function to calculate the mean and a 95\% confidence interval
#'
#' The 95\% confidence interval is based on a normal approximation.
#'
#'
#'@param results  a numeric or logical vector
#'@author  Marsel Scheer
#'@examples
#'
#'meanAndNormCI(rexp(10^4, rate=2))
#'
#'@export
meanAndNormCI <-
function(results){
    m = mean(results)
    s = sd(results)/sqrt(length(results))
    c(mean=m, lower=m-1.96*s, upper=m+1.96*s)
}

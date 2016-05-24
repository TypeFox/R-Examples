#'Computing P-values with a Simulated Sample from the Null Distribution
#'
#'This function computes the p-value of a statistic using a simulated sample
#'from its therotical null distribution.
#'
#'
#'@param s 
#'the observation whose p-value is computed.  For instance a
#'Likelihood Ratio.
#'
#'@param theo_dist 
#'the sample of the distribution under the null hypothesis.
#'
#'@return The p-value associated to the observation \code{s}.
#'
#'@author Boris P. Hejblum
#'
#'@seealso \code{\link{rmixchisq}}
#'
#'@keywords internal
#'
#'
#'@examples
#'
#'\dontrun{ 
#'theo_dist <- rnorm(n=10000, mean=0, sd=1)
#'TcGSA:::pval_simu(s=1.96, theo_dist)
#'1-pnorm(q=1.96, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
#'}
#'
pval_simu <-
function(s,theo_dist){
  1-length(which(theo_dist<s))/length(theo_dist)
}

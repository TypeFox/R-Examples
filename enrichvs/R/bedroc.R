# Calculate a BEDROC for the highly-ranked compounds
# x: a vector for scores
# y: a vector for labels
#
# e.g.)
# > x <- rnorm(100001) - 1:100001 * 0.00005
# > y <- c(rep(1,501), rep(0,length(x)-501))
# > bedroc(x, y, alpha=20.0, decreasing=TRUE)
#
# Ref.: Truchon et al. Evaluating Virtual Screening Methods: 
#   Good and Bad Metrics for the "Early Recognition" Problem.
# J. Chem. Inf. Model. (2007) 47, 488-508.
#

bedroc <- function(x, y, decreasing=TRUE, alpha=20.0) {
  if ( length(x) != length(y) ){
    stop(paste("The number of scores must be equal to the number of labels."))
  }
  N <- length(y)
  n <- length( which(y==1) )
  ord <- order(x, decreasing=decreasing)
  m_rank <- which( y[ord] == 1 )
  s <- sum( exp(-alpha * m_rank / N ) )
  ra <- n / N
  ri <- (N - n) / N
  random_sum <- ra * exp( -alpha / N )*(1.0 - exp( -alpha ) )/ ( 1.0 - exp( -alpha / N ) )
  fac <- ra * sinh( alpha / 2.0 ) / ( cosh(alpha / 2.0 ) - cosh(alpha / 2.0 - alpha * ra ) )
  cte = 1.0/ ( 1 - exp( alpha * ri) )
  return( s / random_sum * fac + cte )
}

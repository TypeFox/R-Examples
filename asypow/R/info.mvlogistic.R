info.mvlogistic <- function( coef, design, rss=1 ){
#     Information Matrix for a Multivariate Logistic

# MODEL IS:
#        Prob event(u) = exp(u) / ( 1 + exp(u) )
# where u = Sum( coef[i]*x[i] )
# and coef is the argument so named.  The x are rows of the 
# design matrix.


# coef - Vector of length p (number of covariates) giving coefficients
#        of variables

# design - Matrix of length n X p each row of which gives values of
#          covariates at one of the n design points
# Note: Most models will include a constant term and the column of
#       design corresponding to this term will be identically 1

# rss - The relative sample size at each design point.  If changed
#       from its default value should be a vector of length n.
# Note:  The information matrix is calculated for one observation
#        spread over the n design points in proportions determined
#        by rss.

# NOTE: The primary use of this routine will be for tables analysed
#       using anova techniques with a logistic model.

#----------------------------------------------------------------------

p <- length( coef )

if (is.vector(design)) {
   design <- matrix(design,nrow=1)
}

tmp <- dim( design )

if ( tmp[2] != p ) 
   stop('Number of variables in coef and design differ')

n <- tmp[1]

lrss <- length(rss)
if ( lrss == 1 ) {
   wt <- rep(1,n)
   tot <- n } else if ( lrss == n ) {
   wt <- rss
   tot <- sum(rss) } else 
         stop('Length of rss does not match number of design points.')

info <- matrix( 0, ncol=p, nrow=p )

for (i in 1:n) {
  
  row <- design[i,]

  u <- crossprod( row, coef )
  dim(u) <- NULL

  p <- invlogit( u )

  mat <- outer( row, row )

  info <- info + wt[i] * p * (1-p) * mat     }

return( info/tot )     }

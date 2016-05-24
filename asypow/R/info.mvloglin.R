info.mvloglin <- function( coef, design, rss=1 ){
#     Information Matrix for a Log-Linear model

# coef - Vector of length p (number of covariates) giving coefficients
#        of variables

# MODEL IS:
#        Prob event(u) = exp(u)
# where u = Sum(log( coef[i] )*x[i] )

# design - Matrix of length n X p each row of which gives values of
#          covariates (x values) at one of the n design points
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

model.info <- matrix( 0, ncol=p, nrow=p )

lcoef <- log(coef)

for (i in 1:n) {
  
  row <- design[i,]

  u <- crossprod( row, lcoef )
  dim(u) <- NULL
  p <- exp( u )
  if (p >= 1) 
     stop(paste("Probability p generated with from design point",i,
                "is not an element of (0,1)"))

  mat <- outer( row/coef, row/coef )

  mult <- (wt[i] * (p / (1-p)))

  model.info <- model.info + mult * mat     }

return( model.info/tot )     }

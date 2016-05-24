ci.cum <-
function( obj,
      ctr.mat = NULL,
       subset = NULL,
         intl = 1,
        alpha = 0.05,
          Exp = TRUE,
       ci.Exp = FALSE,
       sample = FALSE )
{
# First extract all the coefficients and the variance-covariance matrix
cf  <- COEF( obj )
vcv <- VCOV( obj )

# Check if the intervals matches ctr.mat
if( length( intl ) == 1 ) intl <- rep( intl, nrow( ctr.mat ) )
if( length( intl ) != nrow( ctr.mat ) ) stop( "intl must match ctr.mat" )

# Workaround to expand the vcov matrix with 0s so that it matches
# the coefficients vector in case of (extrinsic) aliasing.
if( any( is.na( cf ) ) )
  {
vM <- matrix( 0, length( cf ), length( cf ) )
dimnames( vM ) <- list( names( cf ), names( cf ) )
vM[!is.na(cf),!is.na(cf)] <- vcv
cf[is.na(cf)] <- 0
vcv <- vM
   }

if( is.character( subset ) ) {
  sb <- numeric(0)
  for( i in 1:length( subset ) ) sb <- c(sb,grep( subset[i], names( cf )  ))
  subset <- sb # unique( sb )
  }
# If subset is not given, make it the entire set
if( is.null( subset ) ) subset <- 1:length( cf )
# Exclude units where aliasing has produced NAs.
# Not needed after replacement with 0s
# subset <- subset[!is.na( cf[subset] )]
 cf <-  cf[subset]
vcv <- vcv[subset,subset]
if( is.null( ctr.mat ) )
  {
  ctr.mat <- diag( length( cf ) )
  rownames( ctr.mat ) <- names( cf )
  }
if( dim( ctr.mat )[2] != length(cf) )
    stop( paste("\n Dimension of ", deparse(substitute(ctr.mat)),
          ": ", paste(dim(ctr.mat), collapse = "x"),
          ", not compatible with no of parameters in ",
          deparse(substitute(obj)), ": ", length(cf), sep = ""))

# Finally, here is the actual computation of the estimates
    ct <- ctr.mat %*% cf
    vc <- ctr.mat %*% vcv %*% t( ctr.mat )
# If a sample is requested replace the eatimate by a sample
    if( sample ) ct <- t( mvrnorm( sample, ct, vc ) )
# If Exp was requested, we take the exponential of the estimates
# before we cumulate the sum
if( Exp )
  {
  ct <- exp( ct )
  vc <- ( ct[,1] %*% t(ct[,1]) ) * vc
  }
# Here is the cumulation matrix
  cum.mat <- 1 - upper.tri( diag(ct[,1]) )
  # Multiply columns of the matrix with interval lengths
  cum.mat <- t( intl * t( cum.mat ) )
# This is then multiplied to the coefficients
    ct <- cum.mat %*% ct
    if( sample ) return( ct )
    else
      {
      vc <- cum.mat %*% vc %*% t( cum.mat )
      se <- sqrt( diag( vc ) )
      if( !ci.Exp )
        {
        cum <- cbind( ct, se ) %*% ci.mat( alpha=alpha )
        return( cbind( cum, "StdErr"=se ) )
        }       
      else
        {  
        cum <- exp( cbind( log(ct), se/ct ) %*% ci.mat( alpha=alpha ) )
        return( cbind( cum, "Erf"=exp( qnorm(1-alpha/2)*se/as.vector(ct) ) ) )
        } 
      }
}

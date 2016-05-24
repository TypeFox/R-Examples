# The coef() methods in nlme and lme4 do something different
# so we make a workaround by specifying our own generic methods
COEF          <- function( x, ... ) UseMethod("COEF")
COEF.default  <- function( x, ... ) coef( x, ... )
VCOV          <- function( x, ... ) UseMethod("VCOV")
VCOV.default  <- function( x, ... ) vcov( x, ... )

# Then we can get from these methods what we want from lme, mer etc.
COEF.lme      <- function( x, ... ) nlme::fixed.effects( x )
COEF.mer      <- function( x, ... ) lme4::fixef( x )
COEF.lmerMod  <- function( x, ... ) lme4::fixef( x )
# The vcov returns a matrix with the wrong class so we strip that:
VCOV.lme      <- function( x, ... ) as.matrix(vcov( x ))
VCOV.mer      <- function( x, ... ) as.matrix(vcov( x ))
VCOV.lmerMod  <- function( x, ... ) as.matrix(vcov( x ))

# For the rest of the non-conforming classes we then just need the methods not defined
COEF.crr      <- function( object, ... ) object$coef
VCOV.crr      <- function( object, ... ) object$var
COEF.MIresult <- function( object, ... ) object$coefficients
VCOV.MIresult <- function( object, ... ) object$variance
COEF.mipo     <- function( object, ... ) object$qbar
VCOV.mipo     <- function( object, ... ) object$t
COEF.polr     <- function( object, ... ) summary(object)$coefficients
VCOV.gnlm     <- function( object, ... ) object$cov
VCOV.rq       <- function( object, ... ) summary(object, cov=TRUE)$cov

ci.lin <-
function( obj,
      ctr.mat = NULL,
       subset = NULL,
       subint = NULL,
        diffs = FALSE,
         fnam = !diffs,
         vcov = FALSE,
        alpha = 0.05,
           df = Inf,
          Exp = FALSE,
       sample = FALSE )
{
# First extract all the coefficients and the variance-covariance matrix
cf  <- COEF( obj )
vcv <- VCOV( obj )
# Workaround to expand the vcov matrix with 0s so that it matches
# the coefficients vector in case of (extrinsic) aliasing.
if( any( is.na( cf ) ) )
  {
if( inherits( obj, c("coxph") ) )
  { # aliased parameters are only NAs in coef, but omitted from vcov
  wh <- !is.na(cf)
  cf <- cf[wh]
  vcv <- vcv[wh,wh]
  }
else
if( inherits( obj, c("clogistic") ) )
  {
  cf[is.na(cf)] <- 0
  }
else
  {
vM <- matrix( 0, length( cf ), length( cf ) )
dimnames( vM ) <- list( names( cf ), names( cf ) )
vM[!is.na(cf),!is.na(cf)] <- vcv
cf[is.na(cf)] <- 0
vcv <- vM
  }
  }

# Function for computing a contrast matrix for all possible
# differences between a set of parameters.
all.dif <-
function( cf, pre=FALSE )
{
nn <- length( cf )
nr <- nn * ( nn - 1 ) / 2
nam <- names( cf )

# Work out the indexes of parameter pairs to compare
#
xx <- numeric( 0 )
for( j in 2:nn ) xx <- c(xx, j:nn )
ctr <- cbind( rep( 1:(nn-1), (nn-1):1 ), xx )

# Now for the annotation:
# Find out how large a proportion of rownames are identical
i <- 1
while( all( substr( nam, 1, i ) == substr( nam[1], 1, i ) ) ) i <- i+1

# If a factor name is given, then use this, otherwise the identical part
# of the parameter names
if( is.character( pre ) ) {
prefix <- pre
  pre <- TRUE
} else {
prefix <- substr( nam[1], 1, i-1 )
}
rn <- paste( if( pre ) prefix else "",
             substring( nam[ctr[,1]], i ), "vs.",
             substring( nam[ctr[,2]], i ) )

# Finally, construct the contrast matrix and attach the rownames
cm <- matrix( 0, nr, nn )
cm[cbind(1:nr,ctr[,1])] <- 1
cm[cbind(1:nr,ctr[,2])] <- -1
rownames( cm ) <- rn
cm
# end of function for all differences 
}

# Were all differences requested?
#
if( diffs )
  {
  if( is.character( subset ) )
    {
    if ( inherits( obj, "lm" ) &
         length( grep( subset, names( obj$xlevels ) ) )>0 )
       { # The case of factor level differences we find the relevant
         # subset of parameters by reconstructing names of parameters
         wf <- grep( subset, af <- names( obj$xlevels ) )
         # All factor levels
         fn <- obj$xlevels[[af[wf]]]
         # Reconstruct names of relevant parameter names
         pnam <- paste( af[wf], fn, sep="" )
         # Find them in the parameter vector
         wh <- match( pnam, names( coef( obj ) ) )
         # Get the relevant subset, and stick in 0s for NAs
         cf <- coef( obj )[wh]
         cf[is.na( cf )] <- 0
         vcv <- vcov( obj )[wh,wh]
         vcv[is.na( vcv )] <- 0
         names( cf ) <- rownames( vcv ) <- colnames( vcv ) <-
             paste( subset, ": ", fn, sep="" )
       } else {
         subset <- grep( subset, names( cf ) )
             cf <-  cf[subset]
            vcv <- vcv[subset,subset]
       }
    } else {
     cf <-  cf[subset]
    vcv <- vcv[subset,subset]
   }
   ctr.mat <- all.dif( cf, pre=fnam )
  }

if( !diffs )
  {
  if( is.character( subset ) ) {
    sb <- numeric(0)
    for( i in 1:length( subset ) ) sb <- c(sb,grep( subset[i], names( cf )  ))
    subset <- sb # unique( sb )
    }
  if( is.character( subint ) ) {
    sb <- 1:length(cf)
    for( i in 1:length(subint) ) sb <- intersect( sb, grep(subint[i],names(cf)) )
    subset <- sb # unique( sb )
    }
  if( is.null( subset ) & is.null( subint ) ) subset <- 1:length( cf )
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
  }
# Finally, here is the actual computation
  if( sample )
    {
    # mvrnorm() returns a vector if sample=1, otherwise a sample x
    # length(cf) matrix - hence the rbind so we always get a row 
    # matrix and res then becomes an nrow(ctr.mat) x sample matrix
    res <- ctr.mat %*% t( rbind(mvrnorm( sample, cf, vcv )) )   
    }  
  else
    {
    ct <- ctr.mat %*% cf
    vc <- ctr.mat %*% vcv %*% t( ctr.mat )
    se <- sqrt( diag( vc ) )
    ci <- cbind( ct, se ) %*% ci.mat( alpha=alpha, df=df )
    t0 <- cbind( se, ct/se, 2 * ( pnorm( -abs( ct / se ) ) ) )
    colnames(t0) <- c("StdErr", "z", "P")
    res <- cbind(ci, t0)[, c(1, 4:6, 2:3), drop=FALSE]
    if( Exp )
      {
      res <- cbind(      res[,1:4     ,drop=FALSE],
                    exp( res[,c(1,5,6),drop=FALSE] ) )
      colnames( res )[5] <- "exp(Est.)"
      }
    }
# Return the requested structure
if( sample ) invisible( res ) else
if( vcov ) invisible( list( est=ct, vcov=vc ) ) else res
}

# Handy wrapper
ci.exp <-
function( ..., Exp=TRUE, pval=FALSE )
{
if( Exp )
ci.lin( ..., Exp=TRUE  )[,if(pval) c(5:7,4)   else 5:7     ,drop=FALSE]
else
ci.lin( ..., Exp=FALSE )[,if(pval) c(1,5,6,4) else c(1,5,6),drop=FALSE]
}

# Wrapper for predict.glm
ci.pred <-
function( obj, newdata,
         Exp = NULL,
       alpha = 0.05,
          df = Inf )
{
if( !inherits( obj, "glm" ) ) stop("Not usable for non-glm objects")
# get prediction and se on the link scale
zz <- predict( obj, newdata=newdata, se.fit=TRUE, type="link" )
# compute ci on link scale
zz <- cbind( zz$fit, zz$se.fit ) %*% ci.mat( alpha=alpha, df=df )
# transform as requested
if( missing(Exp) ) {   return( obj$family$linkinv(zz) )
} else {  if(  Exp ) { return(                exp(zz) ) 
   } else if( !Exp )   return(                    zz  )
       }  
}

Deming <-
function( x, y, vr=sdr^2, sdr=sqrt(vr), boot=FALSE, keep.boot=FALSE, alpha=0.05 )
{
if( missing( vr ) & missing( sdr ) ) var.ratio <- 1
else var.ratio <- vr
vn <- c( deparse( substitute( x ) ),
         deparse( substitute( y ) ) )
pn <- c( "Intercept", "Slope", paste( "sigma", vn, sep="." ) )

alfa <- alpha
dfr <- data.frame( x=x, y=y )
dfr <- dfr[complete.cases(dfr),]
x <- dfr$x
y <- dfr$y
n <- nrow( dfr )
SSDy <- var( y )*(n-1)
SSDx <- var( x )*(n-1)
SPDxy <- cov( x, y )*(n-1)
beta <- ( SSDy - var.ratio*SSDx +
          sqrt( ( SSDy - var.ratio*SSDx )^2 +
                4*var.ratio*SPDxy^2 ) ) / ( 2*SPDxy)
alpha <- mean( y ) - mean( x ) * beta
ksi <- ( var.ratio*x + beta*(y-alpha) )/(var.ratio+beta^2)
sigma.x <- ( var.ratio*sum( (x-ksi)^2 ) +
                       sum( (y-alpha-beta*ksi)^2 ) ) /
# The ML-estiamtes requires 2*n at this point bu we do not think we have that
# many observations so we stick to (n-2). Any corroboation from litterature?
           ( (n-2)*var.ratio )
sigma.y <- var.ratio*sigma.x
sigma.x <- sqrt( sigma.x )
sigma.y <- sqrt( sigma.y )
if( !boot ){
res <- c( alpha, beta, sigma.x, sigma.y )
names( res ) <- pn
res
}
else
{
if( is.numeric( boot ) ) N <- boot else N <- 1000
res <- matrix( NA, N, 4 )
for( i in 1:N )
   {
   wh <- sample( 1:n, n, replace=TRUE )
   res[i,] <- Deming( x[wh], y[wh], vr=var.ratio, boot=FALSE )
   }
ests <- cbind( c(alpha,beta,sigma.x, sigma.y),
               se <- sqrt( diag( cov( res ) ) ),
               t( apply( res, 2, quantile, probs=c(0.5,alfa/2,1-alfa/2 ), na.rm=T ) ) )
colnames( res ) <- rownames( ests ) <- pn
colnames( ests )<- c("Estimate", "S.e.(boot)", colnames(ests)[3:5] )
if(keep.boot)
  {
  print( ests )
  invisible( res )
  }
else
  {
  cat( vn[2], " = alpha + beta*", vn[1], "\n" )
  ests
  }
}
}
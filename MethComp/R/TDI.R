TDI <-
function( y1, y2, p = 0.05, boot=1000, alpha = 0.05 )
{
  if( length(y1) != length(y2) )
    stop( "Lengths of y1 and y2 must be the same!" )

  # Analytical approximation
  mu  <- mean(y1-y2)
  sigma <- sd(y1-y2)
  TDI.appr <- qnorm(1-p/2)*sqrt(mu^2+sigma^2)
  names( TDI.appr ) <- "(approximate)"

  # Limits of agreement
  LoA <- mu + c(-1,1)*qnorm(1-p/2)*sigma
  names( LoA ) <- c("lower","upper")

  # Numerical calculation
  FF <- function( x ) pnorm( (x-mu)/sigma ) - pnorm( (-x-mu)/sigma ) - (1-p)
  int <- 2 * max(abs(LoA)) * c(-1,1)
  TDI.num <- uniroot( FF, int )$root
  names( TDI.num ) <- "(numeric)"

  # Bootstrap c.i.
  if( is.numeric(boot) )
    {
    nn <- length(y1)
    tdi <- numeric(boot)
    for( i in 1:boot )
       {
       wh <- (1:nn)[sample(1:nn,nn,replace=T)]
       tdi[i] <- TDI( y1[wh], y2[wh], p=p, boot=FALSE )[[1]]
       }
    TDI.num <- c(TDI.num, quantile( tdi, c(0.5,alpha/2,1-alpha/2) ))
    }

  # Put it all together
  res <- list( c(TDI.appr,TDI.num), LoA )
  names( res ) <- c( paste((1-p)*100,"% TDI",sep=""),
                     paste((1-p)*100,"% Limits of Agreement",sep="") )
  return( res )
}
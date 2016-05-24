bond.cir <-
function( alpha, beta, sigma,  q1, q2, r0, n, maturities, days = 360){
  # Simulation of the values and yields of zero-coupon bonds
  # when the (annualized ) spot rate (in percent) is modeled by a 
  # Feller process satisfying
  #
  #                dr = alpha(beta-r)dt + sigma sqrt(r) dW,   
  #
  #  with market price of risk q = q1/sqrt(r) +q2 sqrt(r). The maturities are 1,3,6 ad 12
  #  months.
  #
  #  Input
  #    r0: initial value
  #    n: number of values
  #   days: days in a year convention (360 default)
  #   maturities: maturities in years (row vector).
  #
  #
  #   Output
  #      P: bond values
  #      R: Annual rate for the bond
  #      tau: maturities in years.
  #
  # Example:
  # out = bond.cir(0.5,2.55,0.365,0.3,0,3.55,1080,c(1/12, 3/12, 6/12, 1),365)
  
  nbBonds <- length(maturities)
  tau <- matrix( 1, n+1, 1) %*% maturities # maturity in years
  
  scalingFact   <- 100
  
  # simulate short rate (%) on an annual time scale 
  r <- sim.cir( alpha, beta, sigma, r0, n ,1/days) 
  
  
  # get pricing coefficients
  params <- get.cir.param( c(alpha,beta,sigma,q1,q2), tau, scalingFact )
  A <- params$A
  B <- params$B
  
  P <- 100 * exp( A - B * t((matrix(1, nbBonds,1) %*% r)) / scalingFact )
  R <- ( B * t((matrix( 1,nbBonds,1) %*% r)) - scalingFact * A ) / tau

  month <- 12*maturities
  library('ggplot2')
  library('reshape')
  tmp <- data.frame(x=(1:nrow(P)), P=P)
  framed.data = melt(tmp, id='x')
  # To remove the notes by R CMD check!
  value <- NULL 
  variable <- NULL
  x <- NULL
  Days <- NULL
  values.graph <- ggplot(framed.data,
                         aes(x=x, y=value, group=variable, colour=variable))
  
  values.graph <- values.graph + geom_line() + 
    ggtitle(sprintf('Daily values of a zero-coupon bond with $100 face value for different maturities'))
  
  print(values.graph)
  
  tmp <- data.frame(Days=(1:nrow(R)), '%'=R)
  
  names(tmp) <- c('Days',sprintf('\tau <-  %d months' ,month))
  framed.data = melt(tmp, id='Days')
  
  values.graph <- ggplot(framed.data,
                         aes(x=Days, y=value, group=variable, colour=variable))
  
  values.graph <- values.graph + geom_line() + 
    ggtitle(sprintf('Annual yield of a zero-coupon bond for different maturities'))
  
  print(values.graph)
  
  return(list(P=P, R=R, tau=tau, r=r))
}

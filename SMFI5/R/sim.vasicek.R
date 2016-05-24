sim.vasicek <-
function( alpha, beta, sigma, r0, n, h){
  # Simulation of a Ornstein-Uhlenbeck process satisfying
  #
  #    dr <- alpha(beta-r)dt + sigma  dW,
  #
  # observed at periods h, 2h, ..., nh.
  #
  #  Input
  #    r0: initial value
  #     n: number of values
  #     h: time step between observations.
  #
  #   Output
  #      r: annual rate
  #
  # r = sim.vasicek(0.5,2.55,0.365,2.55,360,1/360)
  
  phi <- exp( - alpha * h)
  vol <- sigma * sqrt( ( 1 - phi^2 ) / ( 2 * alpha ) )
  
  r <- mat.or.vec(n+1,1)
  
  r[1] <- r0
  
  for(i in 2:(n+1)){
  r[i]  <- beta + phi * ( r[i-1] - beta ) + vol * rnorm(1)
  }
  
  t <- (1:n)*h
  
  t <- c(0,t)
  
  #library('ggplot2')
  #library('reshape')
  tmp <- data.frame(t=t, r=r)
  names(tmp) <- c('t','Annual rate in %')
  framed.data = melt(tmp, id='t')
  
  # To remove the notes by R CMD check!
  value <- NULL 
  variable <- NULL
  
  values.graph <- ggplot(framed.data,
                         aes(x=t, y=value, group=variable, colour=variable))
  
  values.graph <- values.graph + geom_line() + 
    ggtitle(sprintf('OU model for the spot rate with alpha = %f, beta = %f, sigma = %f',alpha,beta,sigma))
  
  print(values.graph)
  
  return(r)

}

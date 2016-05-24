sim.cir <-
function( alpha, beta, sigma, r0, n, h){
  # Simulation of a Feller process satisfying
  # dr <- alpha(beta-r)dt + sigma sqrt(r) dW   
  #
  #  Input
  #    r0: initial value
  #     n: number of values
  #     h: time step between observations.
  #
  #
  #   Output
  #      r: annual rate in percent
  #
  # r = sim.cir(0.5,2.55,0.365,2.55,720,1/360)
  
  # precomputations
  
  sigmaSquared <- sigma^2
  nu           <- 4 * alpha * beta / sigmaSquared 
  phi          <- exp( - alpha * h )
  omega        <- sigmaSquared * ( 1-phi )  / ( 4 * alpha )
  
  r    <- mat.or.vec(n+1,1)
  r[1] <- r0
  
  for(t in 2:(n+1)){
    x <- r[t-1] / omega
    D <- x * phi  # non-centrality parameter
    tt <- sim.n.chi2( nu, D )
    r[t] <- omega * tt
     #r[t] <- omega * rchisq(1, nu, D )
  }
  
  
  #library('ggplot2')
  #library('reshape')
  tmp <- data.frame(x=(1:length(r)+1), r=r)
  names(tmp) <- c('x','Annual rate in %')
  framed.data = melt(tmp, id='x')
  
  # To remove the notes by R CMD check!
  value <- NULL 
  variable <- NULL
  x <- NULL
  
  values.graph <- ggplot(framed.data,
                         aes(x=x, y=value, group=variable, colour=variable))
  
  values.graph <- values.graph + geom_line() + 
    ggtitle(sprintf('CIR model for the spot rate with alpha = %f, beta = %f, sigma = %f, nu = %f',alpha,beta,sigma,nu))
  
  print(values.graph)
  
  return(r)
}

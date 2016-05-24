# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
# ** Copyright UCAR (c) 1992 - 2004 
# ** University Corporation for Atmospheric Research(UCAR) 
# ** National Center for Atmospheric Research(NCAR) 
# ** Research Applications Program(RAP) 
# ** P.O.Box 3000, Boulder, Colorado, 80307-3000, USA 
# ** 2004/9/1 14:13:55 
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=* 
crps <- function(obs, pred, ...)
  ## Tilmann Gneiting's crps code, assumes pred is either a vector of length
  ## 2 (mu, sig) or a matrix of mu and sig if each forcast is different
  {
  if(is.null( dim(pred)) & length(pred)==2){mu <- pred[1];
                                            sigma <- pred[2]} else {
    mu<- as.numeric( pred[,1] ); sigma <- as.numeric( pred[,2]) }
    
  z <- (obs-mu)/sigma ## center and scale
  
  crps<- sigma * (z*(2*pnorm(z,0,1)-1) + 2*dnorm(z,0,1) - 1/sqrt(pi))
  ign <-  0.5*log(2*pi*sigma^2) + (obs - mu)^2/(2*sigma^2)
  pit <- pnorm(obs, mu,sigma )

return(list(crps = crps, CRPS = mean(crps), ign = ign, IGN = mean(ign), pit = pit) )

}


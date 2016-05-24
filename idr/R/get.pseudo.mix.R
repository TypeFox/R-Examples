get.pseudo.mix <-
function(x, mu, sigma, rho, p){

  
  # first compute cdf for points on the grid
  # generate 200 points between [-3, mu+3*sigma]
  nw <- 1000
  w <- seq(min(-3, mu-3*sigma), max(mu+3*sigma, 3), length=nw) 
  w.cdf <- p*pnorm(w, mean=mu, sd=sigma) + (1-p)*pnorm(w, mean=0, sd=1)

  i <- 1

  quan.x <- rep(NA, length(x))

  for(i in c(1:nw)){
    index <- which(x >= w.cdf[i] & x < w.cdf[i+1])
    quan.x[index] <- (x[index]-w.cdf[i])*(w[i+1]-w[i])/(w.cdf[i+1]-w.cdf[i]) +w[i]
  }

  index <- which(x < w.cdf[1])
  if(length(index)>0)
    quan.x[index] <- w[1]

  index <- which(x > w.cdf[nw])
  if(length(index)>0)
    quan.x[index] <- w[nw]  
  
  invisible(quan.x)
}


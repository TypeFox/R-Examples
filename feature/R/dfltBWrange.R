######### R-function dfltBWrange  #########

# Obtain default set of grid counts from a
# multivariate point cloud 'x'.

# Last changed: 22 JUL 2005

dfltBWrange <- function(x,tau) {
  d <- ncol(x)
  if (d==1) x <- as.matrix(x)

  r <- 2
  cmb.fac.upp <- (4/((d+2*r+2)*nrow(x)))^(1/(d+2*r+4))
  r <- 0
  cmb.fac.low <- (4/((d+2*r+2)*nrow(x)))^(1/(d+2*r+4))

  ## Compute the scale in each direction

  st.devs <- apply(x,2,sd)
  IQR.vals <- apply(x, 2, IQR)/(qnorm(3/4) - qnorm(1/4))
  sig.hats <- apply(cbind(st.devs,IQR.vals),1,min)
  ##range.vals <- apply(x,2,max) - apply(x,2,min)

  range.h <- list()
  for (id in 1:d)
  {
    h.upp <- cmb.fac.upp*sig.hats[id]
    h.low <- 0.1*cmb.fac.low*sig.hats[id] ##3*(range.vals[id] + 2*tau*h.upp)/((gridsize[id]-1)*tau)
    range.h[[id]] <- c(h.low,h.upp)
  }

  return(range.h)
}

######## End of dfltBWrange ########

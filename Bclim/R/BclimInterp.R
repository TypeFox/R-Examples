
BclimInterp <- function(Bclim.data, Bclim.res, tgrid=seq(0,14,by=0.1)) {
	# Get relevant parts from data given in arguments
  chron <- Bclim.res$chron.store  # n.s by n where n.s is the numer of stored MCMC iterations
  alpha <- sqrt(Bclim.res$phi2.store/Bclim.res$phi1.store) # Need this parameterisation for bridging both are n.s by m
  delta <- sqrt(Bclim.res$phi1.store*Bclim.res$phi2.store)
	clim <- Bclim.res$c.store   # n.s by n by m
	v <- Bclim.res$v.store  # n.s by n-1 by m
  eps = 1e-5 # A small number
  
  # Get the number of interpolated samples we need and the number of sections
  n.s = dim(clim)[1] # Number of samples
  n = dim(clim)[2] # Number of layers
  m = dim(clim)[3] # Number of climate dimensions
  n.g = length(tgrid) # Number of grid points
  
  # Create some arrays to hold the interpolated climates and squared volatilities v
  # Note that the storage here is longer as we have to interpolate onto a finer grid which includes the chronologies
  clim.interp <- array(NA, dim = c(n.s, n.g, m)) 
  v.interp.all <- array(NA, dim = c(n.s , n.g-1+n, m))
  v.interp <- array(NA, dim = c(n.s, n.g-1, m))
  
  # Extrapolation function
  NIGextrap = function(c.s,t.s,t.out,phi1,phi2,future=FALSE) {
    #c.s=curr.clim[1];t.s=curr.chron[1];t.out=t.select;phi1=Bclim.res$phi1[1];phi2=Bclim.res$phi2[1]
    t.diff = abs(diff(sort(c(t.out,t.s))))
    mu = phi1*t.diff # Re-parameterise on to the R version of the IG distribution
    lambda = phi2*phi1*t.diff^2
    v.out = statmod::rinvgauss(length(t.diff),mu,lambda)
    if(future==TRUE) {
      # extrapolating into the future
      return(list(NIG=rev(cumsum(rev(rnorm(length(t.diff),mean=0,sd=sqrt(v.out)))))+c.s,IG=v.out))
    } else {
      return(list(NIG=cumsum(rnorm(length(t.diff),mean=0,sd=sqrt(v.out)))+c.s,IG=v.out))
    }
  }
    
  # Bridging function - taken from Ribeiro 2008
  NIGB <- function(delta, IG.sum, tb.points, c.start, c.end){
    total.points <- length(tb.points)
    NIG.bridge <- rep(NA,total.points)
    NIG.bridge[1] <- c.start
    NIG.bridge[total.points] <- c.end
    IG.increment <- rep(NA, total.points-2)
    z <- IG.sum
    eps = 1e-5
    
    #curr.cstart = c.start
    
    l <- 2
    while(l < total.points){
      
      if((tb.points[l]-tb.points[l-1])<eps) {
        IG.increment[l-1] = 0 
      } else {   
        
        # Generate chi-square random variate
        q <- rnorm(1,0,1)^2  
        
        # Reparameterise
        mu <- (tb.points[total.points]-tb.points[l]) / (tb.points[l]-tb.points[l-1])   
        #lambda <- (1/delta^2 * (tb.points[total.points]-tb.points[l])^2) / z
        lambda <- (delta^2 * (tb.points[total.points]-tb.points[l])^2) / z
        
        if(z==0) {
          print(c(delta,IG.sum,tb.points,c.start,c.end))
          stop("Problem in Inverse Gaussian Bridge - IG sum is zero")
        }
        
        # Compute the roots of the chi-square random variate
        s1 <- mu + (mu^2*q)/(2*lambda) - (mu/(2*lambda)) * sqrt(4*mu*lambda*q + mu^2*q^2)
        if(lambda<eps) s1 <- mu
        s2 <- mu^2 / s1
        
        # Acceptance/rejection of root
        s <- ifelse(runif(1,0,1) < mu*(1+s1)/((1+mu)*(mu+s1)), s1, s2)
        
        # Compute the IG incrrement
        IG.increment[l-1] <- z / (1 + s)
        if(any(IG.increment<0,na.rm=TRUE)) stop("Inverse Gaussian bridge producing negative variances")
        
      # End of if statement
      }
           
      # Rescale the sum of left-over distance of the IG bridge
      z <- z - IG.increment[l-1]
              
      #Compute the current value of the NIG bridge
      #NIG.bridge[l] <- (c.start*z + c.end*(IG.sum-z)) / IG.sum + rnorm(1, 0, 1) * sqrt((IG.sum-z)*z / IG.sum)
      #NIG.bridge[l] <- curr.cstart+ (c.end-curr.cstart)*IG.increment[l-1]/(IG.increment[l-1]+z) + rnorm(1) * sqrt(IG.increment[l-1]*z/(z+IG.increment[l-1]))
      #curr.cstart = NIG.bridge[l]
      
      l <- l + 1
      
      # End of while loop
    }
    
    W1 = 1/IG.increment
    W2 = cumsum(1/(IG.sum-IG.increment))
    bridge.var = 1/(W1+W2)
    bridge.mean = c.start+bridge.var*W2*(c.end-c.start) 
    NIG.bridge[2:((total.points)-1)] = rnorm(total.points-2,mean=bridge.mean,sd=sqrt(bridge.var))
    
    return(list(IGB = c(IG.increment, (IG.sum - sum(IG.increment))), NIGB = NIG.bridge))
  }
  
  # Now loop through each climate dimension 
  for (j in 1:Bclim.data$m){  
    cat("Perform interpolation for climate dimension", j, "\n")
    # Now loop through each sample n.s
    for (i in 1:n.s) { #loop through each sample
      # Report on completion every 10 samples
      if(i%%10==0) {
        cat("\r")
        cat("Completed:",format(round(100*i/n.s,2), nsmall = 2),"%")
        if(i<n.s) {
          cat("\r")
        } else {
          cat("\n")
        }
      }
      
      # Get small vectors of the current values
      curr.clim = clim[i,,j] # length n
      curr.v = v[i,,j] # length n-1
      curr.chron = chron[i,] # length n
      tgrid.all = sort(c(curr.chron,tgrid)) # length n + n.g
      
      # If there are any chronology times which exactly match the tgrid then fill in those climates and variances as zero
      if(any(diff(tgrid.all)<eps)) {
        # Just shift the chronology by a little bit
        curr.chron = curr.chron+eps
        tgrid.all = sort(c(curr.chron,tgrid))
      }
      
      # See if extrapolation into the future is required
      if(tgrid[1]<min(curr.chron)) {
        # Get all the bits of tgrid that are less than the smallest value of the chronology
        t.select.index = which(tgrid<min(curr.chron))
        t.all.select.index = which(tgrid.all<min(curr.chron))
        t.select = tgrid[t.select.index]
        
        # Perform extrapolation
        clim.temp <- NIGextrap(curr.clim[1], curr.chron[1], t.select ,Bclim.res$phi1.store[i,j], Bclim.res$phi2.store[i,j],future=TRUE)
        clim.interp[i,t.select.index,j] <- clim.temp$NIG
        v.interp.all[i,t.all.select.index,j] = clim.temp$IG
      }
      
      # See if extrapolation into the past is required
      if(tgrid[length(tgrid)]>max(curr.chron)) {
        t.select.index = which(tgrid>max(curr.chron))
        t.all.select.index = which(tgrid.all>max(curr.chron))
        t.select = tgrid[t.select.index]
        
        # Perform extrapolation
        clim.temp = NIGextrap(curr.clim[length(curr.clim)],curr.chron[length(curr.clim)],t.select,Bclim.res$phi1.store[i,j],Bclim.res$phi2.store[i,j])
        clim.interp[i,t.select.index,j] = clim.temp$NIG
        v.interp.all[i,t.all.select.index-1,j] = clim.temp$IG
      }
  
      # Now look at how many gaps their are
      for(k in 1:(n-1)) {
        # Find out if in this section there are any tgrid points inside
        t.select.index = which(tgrid>=curr.chron[k] & tgrid<curr.chron[k+1])
        t.all.select.index = which(tgrid.all>=curr.chron[k] & tgrid.all<curr.chron[k+1])
        
        # If the length of t.select.index is positive then there are points inside and we should use the NIG bridge
        if(length(t.select.index)>0) {
          # Select which bits of the grid are inbetween the current section of the chronology
          t.select = tgrid[t.select.index]
          
          # Now interpolate using the NIGB code
          clim.temp <- NIGB(delta[i,j], curr.v[k], c(curr.chron[k],t.select,curr.chron[k+1]), curr.clim[k], curr.clim[k+1])
          
          clim.interp[i,t.select.index,j] = clim.temp$NIGB[-c(1,length(clim.temp$NIGB))]
          v.interp.all[i,t.all.select.index,j] = clim.temp$IGB
        } 
        
        # If there's no grid points just store the v value from the original data
        if(length(t.select.index)==0 & length(t.all.select.index)>0) {
          v.interp.all[i,t.all.select.index,j] = curr.v[k]
        }
      }    
      
      # Now sum over the interpolated v.interp.all values so as to get the correct interpolated v values
      tgrid.select.lower = match(tgrid,tgrid.all)
      tgrid.select.upper = c((tgrid.select.lower-1)[2:length(tgrid.select.lower)],length(tgrid.all))
      for(l in 1:(n.g-1)) {
        v.interp[i,l,j] = sum(v.interp.all[i,tgrid.select.lower[l]:tgrid.select.upper[l],j])    
      }
      
      if(any(is.na(clim.interp[i,,j]))) stop("Some interpolated climates have been given NAs")
      if(any(is.na(v.interp[i,,j]))) stop("Some interpolated v values have been given NAs")
    
    } # End of samples loop
  } # End of climate dimension loop
  
  for(j in 1:Bclim.data$m) {
    clim.interp[,, j] <- clim.interp[,, j] * sqrt(Bclim.data$ScVar)[j] + Bclim.data$ScMean[j]
    v.interp[,, j] <- sqrt(v.interp[,, j]) * sqrt(Bclim.data$ScVar)[j]
  }
  
  cat("Completed! \n")
  return(list(clim.interp = clim.interp, v.interp = v.interp, time.grid=tgrid))
}

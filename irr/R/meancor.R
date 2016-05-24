"meancor" <-
function(ratings, fisher=TRUE) {
	ratings <- as.matrix(na.omit(ratings))

	ns <- nrow(ratings)
	nr <- ncol(ratings)

	for (i in 1:(nr-1)) for (j in (i+1):nr) {
    if ((i==1) & (j==(i+1))) r <- cor(ratings[,i],ratings[,j])
    else r <- c(r, cor(ratings[,i],ratings[,j]))
	}

  delr <- 0
  
  if (fisher) {
    delr <- length(r) - length(r[(r<1) & (r>-1)])
    #Eliminate perfect correlations (r=1, r=-1)
    r <- r[(r<1) & (r>-1)]
    
 	  rz  <- 1/2*log((1+r)/(1-r))
 	  mrz <- mean(rz)

 	  coeff <- (exp(2*mrz)-1)/(exp(2*mrz)+1)
  	SE    <- sqrt(1/(ns-3))

  	u       <- coeff/SE
  	p.value <- 2 * (1 - pnorm(abs(u)))
  }
  else {
    coeff <- mean(r)
  }

  rval <- list(method = "Mean of bivariate correlations R",
               subjects = ns, raters = nr,
               irr.name = "R", value = coeff)

  if (fisher) rval <- c(rval, stat.name = "z", statistic = u, p.value = p.value)
  if (delr>0) rval <- c(rval, error = paste(delr, ifelse(delr==1, "perfect correlation was", "perfect correlations were"), "dropped before averaging"))
  class(rval) <- "irrlist"

  return(rval)
}


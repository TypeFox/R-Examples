"meanrho" <-
function(ratings, fisher=TRUE) {
	ratings <- as.matrix(na.omit(ratings))

	ns <- nrow(ratings)
	nr <- ncol(ratings)

	#Test for ties
	TIES = FALSE
	testties <- apply(ratings, 2, unique)
	if (!is.matrix(testties)) TIES <- TRUE
	else { if (length(testties) < length(ratings)) TIES <- TRUE }

	ratings.rank <- apply(ratings,2,rank)

	for (i in 1:(nr-1)) for (j in (i+1):nr) {
    if ((i==1) & (j==(i+1))) r <- cor(ratings[,i], ratings[,j], method="spearman")
    else r <- c(r, cor(ratings[,i], ratings[,j], method="spearman"))
	}

  delr <- 0

  if (fisher) {
    delr <- length(r) - length(r[(r<1) & (r>-1)])
    #Eliminate perfect correlations (r=1, r=-1)
    r <-  r[(r<1) & (r>-1)]

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

  rval <- list(method = "Mean of bivariate rank correlations Rho",
               subjects = ns, raters = nr,
               irr.name = "Rho", value = coeff)

  if (fisher) rval <- c(rval, stat.name = "z", statistic = u, p.value = p.value)
  if (delr>0) {
    if (TIES) {
      rval <- c(rval, error = paste(delr, ifelse(delr==1, "perfect correlation was", "perfect correlations were"), "dropped before averaging",
                                    "\n Coefficient may be incorrect due to ties"))
    }
    else {
      rval <- c(rval, error = paste(delr, ifelse(delr==1, "perfect correlation was", "perfect correlations were"), "dropped before averaging"))
    }
  }
  
  if ((delr==0) & (TIES)) rval <- c(rval, error = "Coefficient may be incorrect due to ties")
  class(rval) <- "irrlist"

  return(rval)
}


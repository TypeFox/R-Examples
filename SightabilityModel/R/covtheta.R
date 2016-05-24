covtheta <-
function(total, srates, stratum, subunit, covars, betas, varbetas, nboots){

  #  Arguments (the first four should have one observation for each individual in the obs. dataset
  #  total = total number of animals in the observed group
  #  srates = sampling rate for the stratum (for each observed animal)
  #  subunit = subunit identifier 
  #  covars = matrix of covariates used in the sightability model
  #  beta = vector of parameter estimates from logistic model
  #  varbeta = var/cov matrix for above parameter estimates
  # nboots = number of bootstrap replicates
  
  # Check on length of vector arguments
    n <- length(total)
    if(length(srates) != n) {stop("Srates vector needs to be same dimension as total vector")}
    if(length(stratum) != n) {stop("Stratum vector needs to be same dimension as total vector")}
    if(length(srates) != n) {stop("Subunit vector needs to be same dimension as total vector")}

  # Make sure strata are numbered 1 to h
    stratum <- as.numeric(factor(stratum))

  # Form xmatrix
    xdata <- as.matrix(cbind(rep(1, n), covars))
    ncovars <- dim(xdata)
    nxdata <- nrow(xdata)
    if(ncovars[1] != n) {stop("Covariate matrix needs to be same length as total vector")}  
      
  # Make sure beta is a matrix of dim ncovar x 1 
    thetas <- matrix(NA, nboots, n)
    for(j in 1:nboots){
      tempbet <- betas[j, ]
      tempvarbeta <- varbetas[j, , ]
      
  #  Estimate theta for each animal = 1/prob(detection)
      thetas[j, ] <- 1+exp(-xdata%*%tempbet-diag(xdata%*%tempvarbeta%*%t(xdata))/2)  
  } 
    
    
  # Use empirical variance for Cov(thetas)
    smat <- cov(thetas) 
    return(smat)
}

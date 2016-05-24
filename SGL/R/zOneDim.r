oneDim <-
function(data, index, thresh = 0.0001, nlam = 20, lambdas = NULL, beta.naught = rep(0,ncol(data$x)), inner.iter = 100, outer.iter = 100, outer.thresh = 0.0001, gamma = 0.8, step = 1, reset = 10, alpha = 0.95, min.frac = 0.05, verbose = FALSE){

  if(is.null(lambdas)){
  lambdas <- betterPathCalc(data = data, index = index, alpha=alpha, min.frac = min.frac, nlam = nlam, type = "linear")
  }

  X <- data$x
  y <- data$y
  n <- nrow(X)
  p <- ncol(X)

  ## Setting up group lasso stuff ##
  
  ord <- order(index)
  index <- index[ord]
  X <- X[,ord]
  unOrd <- match(1:length(ord),ord)

    ## Coming up with other C++ info ##

  groups <- unique(index)
  num.groups <- length(groups)
  range.group.ind <- rep(0,(num.groups+1))
  for(i in 1:num.groups){
    range.group.ind[i] <- min(which(index == groups[i])) - 1
  }
  range.group.ind[num.groups+1] <- ncol(X)

  group.length <- diff(range.group.ind)

  ## DONE SETTING UP C STUFF ##

  #alpha <- sqrt(2*log(p))/(1+sqrt(2*log(num.groups)/min(group.length)) + sqrt(2*log(p)))

  nlam = length(lambdas)
  beta.old <- rep(0,ncol(X))
  beta.is.zero <- rep(1,num.groups)
  beta <- array(0, c(ncol(X),nlam))

  eta <- rep(0,n)

  for(k in 1:nlam){
    beta.is.zero <- rep(1, num.groups)
    beta.old <- rep(0, ncol(X))
    eta <- rep(0,n)
      
    junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(lambdas[k]*alpha), lambda2 = as.double(lambdas[k]*(1-alpha)), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset))
   
    beta.new <- junk$beta
    beta[,k] <- beta.new
    beta.is.zero <- junk$betaIsZero
    eta <- junk$eta
    beta.old <- beta.new
    if(verbose == TRUE){
      write(paste("***Lambda", k, "***"),"")
    }
    }
  return(list(beta = beta[unOrd,], lambdas = lambdas))
}

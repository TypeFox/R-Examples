findNum <- function(data, alpha = 0.95, min.frac = 0.05, nlam = 20, type = "linear", num = 5, del = 0.9){

reset <- 10
step <- 1
gamma <- 0.8

inner.iter <- 1000
outer.iter <- 1000
thresh = 10^(-3)
outer.thresh = thresh
  
  n <- nrow(data$x)

  if(type == "linear"){

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
  range.group.ind[num.groups + 1] <- ncol(X)
  
  group.length <- diff(range.group.ind)
  beta.naught <- rep(0,ncol(X))
  beta <- beta.naught
  
  beta.is.zero <- rep(1, num.groups)
  beta.old <- rep(0, ncol(X))

  beta <- array(0, c(ncol(X),nlam,nlam))
                 #matrix(0, nrow = ncol(X), ncol = nlam)

  eta <- rep(0,n)

  max.lam <- max(t(X)%*%y)/n

  is.nonzero <- 0
  above <- 0
  change <- 0

  move <- 0.99

  while(is.nonzero != 5){

  junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(alpha*max.lam), lambda2 = as.double((1-alpha)*max.lam), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset))

  is.nonzero <- sum(abs(junk$beta))
  if(is.nonzero < num){
  change <- above
  max.lam <- max.lam * move
  above <- 0
  }
  if(is.nonzero > num){
  change <- 1 - above
  max.lam <- max.lam /move
  above <- 1
  }
  if(change == 1){
  move <- move * del
  }
}
}
  return(junk$beta)
}

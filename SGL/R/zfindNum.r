findNum <- function(data, index, alpha = 0.95, type = "linear", num = 5, del = 0.9, thresh = 10^(-3)){

reset <- 10
step <- 1
gamma <- 0.8
inner.iter <- 1000
outer.iter <- 1000

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

  eta <- rep(0,n)

  max.lam <- max(t(X)%*%y)/n

  is.nonzero <- 0
  above <- 0
  change <- 0
  minBelow <- 100

  move <- 0.99

  moveOn <- 0
  count <- 0

  while(moveOn == 0){

  junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(alpha*max.lam), lambda2 = as.double((1-alpha)*max.lam), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset))

    beta.new <- junk$beta
    beta.is.zero <- junk$betaIsZero
    eta <- junk$eta
    beta.old <- beta.new

  is.nonzero <- sum(junk$beta!=0)

  if(is.nonzero < num){
  count <- 0
  minBelow <- max.lam
  max.lam <- max.lam * move
  }
  if(is.nonzero > num){
  count <- count + 1
  if(count > 5){
  moveOn <- 1
  }
  max.lam <- mean(c(max.lam,minBelow))
  }

  if(is.nonzero == num){
  moveOn <- 1
  }

}
}
write(is.nonzero,"")
  return(junk$beta)
}

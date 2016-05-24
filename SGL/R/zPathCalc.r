pathCalc <- function(data, alpha = 0.95, min.frac = 0.05, nlam = 20, type = "linear"){

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

  beta <- matrix(0, nrow = ncol(X), ncol = nlam)

  eta <- rep(0,n)

  max.lam <- max(abs(t(X)%*%y))/n

  is.nonzero <- 0

  while(is.nonzero == 0){

  junk <- .C("linNest", X = as.double(as.vector(X)), y = as.double(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(alpha*max.lam), lambda2 = as.double((1-alpha)*max.lam), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), step = as.double(step), reset = as.integer(reset))

  is.nonzero <- sum(abs(junk$beta))
  max.lam <- max.lam * 0.99
}

  max.lam <- max.lam / 0.99

  min.lam <- min.frac*max.lam

  lambdas <- exp(seq(log(max.lam),log(min.lam), (log(min.lam) - log(max.lam))/(nlam-1)))
    
 #  lambdas <- seq(max.lam,min.lam, (min.lam - max.lam)/(nlam-1))

  }

  if(type == "logit"){

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
  betas <- matrix(0, nrow = ncol(X), ncol = nlam)

  eta <- rep(0,n)

  intercepts <- log(sum(y)) - log(n-sum(y))

  eta = eta + intercepts

  resp <- 4*(y-1/2)

  max.lam <- max(abs(t(X)%*%resp))/n

  is.nonzero <- 0

  while(is.nonzero == 0){

    junk <- .C("logitNest", X = as.double(as.vector(X)), y = as.integer(y), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(alpha*max.lam), lambda2 = as.double(max.lam * (1-alpha)), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), eta = as.double(eta), gamma = as.double(gamma), betaIsZero = as.integer(beta.is.zero), betaZero = as.double(intercepts), step = as.double(step))

  is.nonzero <- sum(abs(junk$beta))
  max.lam <- max.lam * 0.99
}

  max.lam <- max.lam / 0.99

  min.lam <- min.frac*max.lam

  lambdas <- exp(seq(log(max.lam),log(min.lam), (log(min.lam) - log(max.lam))/(nlam-1)))

  }

  if(type == "cox"){

  covariates <- data$x
  n <- nrow(covariates)
  p <- ncol(covariates)  
  time <- data$time
  status <- data$status

  ## Ordering Response and Removing any Censored obs before first death ##
  death.order <- order(time)
  ordered.time <- sort(time)  

  X <- covariates[death.order,]  
  ordered.status <- status[death.order]

  first.blood <- min(which(ordered.status == 1))

  X <- X[first.blood:n,]
  ordered.status <- ordered.status[first.blood:n]
  ordered.time <- ordered.time[first.blood:n]
  death.order <- death.order[first.blood:n]
  n <- n-first.blood+1
  
 death.times <- unique(ordered.time[which(ordered.status == 1)])  ## Increasing list of times when someone died (censored ends not included) ##

  ## Calculating Risk Sets ##
  
  risk.set <- rep(0,n)
  for(i in 1:n){
    risk.set[i] <- max(which(death.times <= ordered.time[i]))
  }

  ## Calculating risk set beginning/ending indices ##
  
  risk.set.ind <- rep(0,(length(death.times)+1))  
  for(i in 1:length(death.times)){
    risk.set.ind[i] <- min(which(ordered.time >= death.times[i]))
  }
  risk.set.ind[length(risk.set.ind)] <- length(ordered.time) + 1

  ## Calculating number of deaths at each death time ##
  num.deaths <- rep(0,length(death.times))
  for(i in 1:length(ordered.time)){
    if(ordered.status[i] == 1){
      num.deaths[which(death.times == ordered.time[i])] <-  num.deaths[which(death.times == ordered.time[i])] + 1
    }
  }
    
    ## Finding death indices and number of deaths ##
    death.index <- which(ordered.status == 1)
    total.deaths <- length(death.index)
    
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

    beta.is.zero <- rep(1, num.groups)
  
  eta <- rep(0,n)

  ## DONE SETTING UP COX MODEL STUFF


junk1 <- .C("Cox", riskSetInd = as.integer(risk.set.ind), riskSet = as.integer(risk.set), numDeath = as.integer(num.deaths), status = as.integer(ordered.status), ndeath = as.integer(length(death.times)), nrow = as.integer(n), ncol = as.integer(p), beta = as.double(rep(0,p)), eta = as.double(rep(0,n)), y = as.double(rep(0,n)), weights = as.double(rep(0,n)))

  resp <- junk1$y * junk1$weights

  max.lam <- max(abs(t(X)%*%resp))/n

  is.nonzero <- 0

   while(is.nonzero == 0){

    junk <- .C("coxSolver", X = as.double(as.vector(X)), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(alpha*max.lam), lambda2 = as.double((1-alpha)*max.lam), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), riskSetInd = as.integer(risk.set.ind), riskSet = as.integer(risk.set), numDeath = as.integer(num.deaths), status = as.integer(ordered.status), ndeath = as.integer(length(death.times)), eta = as.double(eta), gamma = as.double(gamma), deathInd = as.integer(death.index), totDeath = as.integer(total.deaths), betaIsZero = as.integer(beta.is.zero), step = as.double(step))

  is.nonzero <- sum(abs(junk$beta))
  max.lam <- max.lam * 0.99
}

  max.lam <- max.lam / 0.99

  min.lam <- min.frac*max.lam

  lambdas <- exp(seq(log(max.lam),log(min.lam), (log(min.lam) - log(max.lam))/(nlam-1)))

  }

  return(lambdas)
}
               

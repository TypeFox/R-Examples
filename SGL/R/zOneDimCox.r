oneDimCox <-
function(data, index, thresh = 0.0001, lambdas= NULL, beta.naught = rep(0,ncol(data$x)), inner.iter = 100, outer.iter = 100, outer.thresh = 0.0001, gamma = 0.8, step = 1, alpha = 0.95, min.frac = 0.05, nlam = 20, reset = 20, verbose = FALSE){

  if(is.null(lambdas)){
  lambdas <- betterPathCalc(data = data, index = index, alpha=alpha, min.frac = min.frac, nlam = nlam, type = "cox")
  }

    ## SETTING UP COX MODEL STUFF ##
  
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
  eta <- rep(0,n)

  ## Finding death indices and number of deaths ##
  death.index <- which(ordered.status == 1)
  total.deaths <- length(death.index)

  ## DONE SETTING UP COX MODEL STUFF

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

  nlam = length(lambdas)
  beta.old <- rep(0,ncol(X))
  beta.is.zero <- rep(1,num.groups)
  beta <- matrix(0, nrow = ncol(X), ncol = nlam)
  
    beta.is.zero <- rep(1, num.groups)
    beta.old <- rep(0, ncol(X))
    eta <- rep(0,n)
    for(i in 1:nlam){

      junk <- .C("coxSolver", X = as.double(as.vector(X)), index = as.integer(index), nrow = as.integer(nrow(X)), ncol = as.integer(ncol(X)), numGroup = as.integer(num.groups), rangeGroupInd = as.integer(range.group.ind), groupLen = as.integer(group.length), lambda1 = as.double(lambdas[i]*alpha), lambda2 = as.double((1-alpha)*lambdas[i]), beta = as.double(beta.old), innerIter = as.integer(inner.iter), outerIter = as.integer(outer.iter), thresh = as.double(thresh), outerThresh = as.double(outer.thresh), riskSetInd = as.integer(risk.set.ind), riskSet = as.integer(risk.set), numDeath = as.integer(num.deaths), status = as.integer(ordered.status), ndeath = as.integer(length(death.times)), eta = as.double(eta), gamma = as.double(gamma), deathInd = as.integer(death.index), totDeath = as.integer(total.deaths), betaIsZero = as.integer(beta.is.zero), step = as.double(step))
    
    beta.new <- junk$beta
    beta[,i] <- beta.new
    beta.is.zero <- junk$betaIsZero
    eta <- junk$eta
    beta.old <- beta.new
    if(verbose == TRUE){
      write(paste("***Lambda", i, "***"),"")
      }
  }
  return(list(beta = beta[unOrd,], death.times = death.times, ordered.time = ordered.time, X = X, lambdas = lambdas))
}


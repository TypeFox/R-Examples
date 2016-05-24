coxCrossVal <-
function(data, index, nfold = 10, nlam = 20, lambdas = lambdas, min.frac = 0.05, alpha = 0.95, maxit = 10000, gamma = 0.8, thresh = 0.0001, verbose = TRUE, step = 1, reset = 10){

  ## Setting up basic stuff

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
  beta.naught <- rep(0,ncol(X))
  beta <- beta.naught

  ## Done with group stuff ##
  
  y <- rep(0,n)
  weights <- rep(1,n)

  ## finding the path

MainSol <- oneDimCox(data, index, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, step = step, reset = reset, alpha = alpha)
  
  lambdas <- MainSol$lambdas

  lldiff <- rep(0, nlam)
  lldiffFold <- matrix(0, nrow = nlam, ncol = nfold)

  size <- floor(nrow(X)/nfold)
  o_flow <- c(rep(1,nrow(X) - size * nfold), rep(0, nfold - (nrow(X) - size * nfold)))
  sizes <- size + o_flow
  ind.split <- c(1,cumsum(sizes))
  
  ind <- sample(1:nrow(data$x), replace = FALSE)
  for(i in 1:nfold){
    ind.out <- ind[ind.split[i]:ind.split[i+1]]
    ind.in <- ind[-(ind.split[i]:ind.split[i+1])]
    new.data <- list(x = data$x[ind.in,], y = data$y[ind.in])
    
    new.data <- list(x = data$x[ind.in,], time = data$time[ind.in], status = data$status[ind.in])

    new.sol <- oneDimCox(new.data, index, thresh = thresh, inner.iter = maxit, lambdas = lambdas, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, gamma = gamma, step = step, reset = reset, alpha = alpha)

	for(k in 1:nlam){
      lldiffFold[k,i] = log.likelihood.calc(X, new.sol$beta[,k], death.times, ordered.time) - log.likelihood.calc(new.sol$X, new.sol$beta[,k], new.sol$death.times, new.sol$ordered.time)

      lldiff[k] <- lldiff[k] + log.likelihood.calc(X, new.sol$beta[,k], death.times, ordered.time) - log.likelihood.calc(new.sol$X, new.sol$beta[,k], new.sol$death.times, new.sol$ordered.time)
      }
    if(verbose == TRUE){
   write(paste("*** NFOLD ", i, "***"),"")
  }
  }
  lldiffSD <- apply(lldiffFold,1,sd) * sqrt(nfold)
  obj <- list(lambdas = lambdas, lldiff = lldiff, llSD = lldiffSD, fit = MainSol)
  class(obj)="cv.SGL"
  return(obj)
}


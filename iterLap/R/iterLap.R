nlm.control <- function(d, typsize = rep(1, d), fscale = 1, print.level = 0,
                ndigit = 12, gradtol = 1e-6, stepmax = max(1000 * sqrt(sum((d/typsize)^2)), 1000),
                steptol = 1e-6, iterlim = 200, check.analyticals = TRUE){
  list(typsize = typsize, fscale = fscale, print.level = print.level,
       ndigit = ndigit, gradtol = gradtol, stepmax = stepmax,
       steptol = steptol, iterlim = iterlim,
       check.analyticals = check.analyticals)
}

nlopt <- function(par, fn, ...,   method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                  control = list(), hesscontrol = list(), lower = -Inf, upper = Inf){
  ## optimizer
  method <- match.arg(method)
  if(method == "nlm"){
    if (!is.null(control)) {
      if (!is.list(control)) {
        stop("when specified, 'nlmcontrol' must be a list")
      }
      pardim <- length(par)
      nlmctrl <- do.call("nlm.control", c(d = pardim, control))
    }
    nlmcall <- list(f = fn, p = par, ..., hessian = TRUE)
    nlmcall <- c(nlmcall, nlmctrl)
    optres <- do.call(nlm, nlmcall)
    par <- optres$estimate
    hess <- optres$hessian
    minimum <- optres$minimum
    if(optres$code < 3){
      conv <- 0
    } else {
      conv <- 1
    }
  } else {
    if(method == "nlminb"){
      optres <- nlminb(start = par, objective = fn, ..., control=control,
                       lower=lower, upper=upper)
      par <- optres$par
      hess <- optimHess(par, fn, gr=NULL, ..., control = hesscontrol)
      hess <- 0.5 * (hess + t(hess))
      minimum <- optres$objective
    } else {
      optres <- optim(par, fn, ..., method = method,
                      control = control, hessian = TRUE)
      par <- optres$par
      hess <- optres$hessian
      minimum <- optres$value
    }
    conv <- optres$convergence
  }
  list(par=par, hessian=hess, minimum=minimum, conv = conv)
}

iterLap.control <- function(pardim, gridSize = ceiling(50*pardim^1.25),
                            delta = 0.01, maxDim = 20, info = 0, eps = 0.005){
  ## gridSize - number of grid points for new components
  ## delta - stopping criterion based on maximum error on the grid
  ## maxDim - Maximum dimension for iterative process
  ## info - Plot information during function execution (possible values 0,1,2)
  ## eps - percentage for stopping criterion based on normalization constant
  list(gridSize=gridSize, delta=delta, maxDim=maxDim, eps=eps, info=info)
}


GRApprox <- function(post, start, grad, 
                     method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                     control = list(), ...){
  ## post - (non-normalized) log-posterior density
  ## start - vector of starting values if dimension=1
  ##         otherwise matrix of starting values with
  ##         the starting values in the rows
  ## grad - gradient (currently not implemented)
  ## control - control parameters for nlopt
  ##
  ## ... - additional arguments for post function
  method <- match.arg(method)
  if(!is.matrix(start)){
    start <- as.matrix(start, ncol=1)
  }
  nStart <- nrow(start)
  d <- ncol(start)
  
  post2 <- function(pars, ...){
    -post(pars, ...)
  }
  obj <- list()
  for(j in 1:nStart){
    nlcall <- list(par=start[j,], fn = post2, ..., method = method)
    nlcall <- c(nlcall, control)
    obj[[j]] <- do.call(nlopt, nlcall)
  }
  means <- matrix(0, nrow = nStart, ncol = d)
  sigmainv <- sigmas <- eigenHess <- list()
  weights <- dets <- numeric(nStart)
  z <- 1
  for(j in 1:nStart){
    eg <- eigen(obj[[j]]$hessian, symmetric = TRUE)
    detSig <- 1/prod(eg$values)
    if(all(eg$values > 0)){ # only include maxima of posterior
      m <- obj[[j]]$par
      sw <- sweep(means, 2, m) # check whether found optimum coincides with already found optimum
      mabs <- max(sum(abs(m)), 1)
      an <- apply(sw, 1, function(x) sum(abs(x)))
      ind1 <- j == 1
      ind2 <- all(an[1:z] > 0.01*mabs)
      if( ind1 | (ind2 & ind2) ){
        means[z,] <- m
        sigmas[[z]] <- solve(obj[[j]]$hessian)
        sigmainv[[z]] <- obj[[j]]$hessian
        eigenHess[[z]] <- eg
        dets[z] <- detSig
        weights[z] <- detSig^0.5*exp(-obj[[j]]$minimum)
        z <- z + 1
      }
    }
  }
  means <- means[1:(z-1), , drop=FALSE]
  weights <- weights[1:(z-1)]/sum(weights[1:(z-1)])
  out <- list(weights=weights, means=means,
              sigmas=sigmas, eigenHess=eigenHess,
              dets = dets[1:(z-1)], sigmainv = sigmainv)
  class(out) <- "mixDist"
  out
}

print.mixDist <- function(x, digits = 4, ...){
  dims <- dim(x$sigmas[[1]])[1]
  ncomp <- length(x$weights)
  cat("Mixture Approximation\n\n")
  cat("Dimensions:", dims, "\n")
  cat("Number of Components:", ncomp, "\n")
  cat("Means:\n")
  mns <- x$means
  nams <- paste("Component ", 1:ncomp, ":", sep="")
  dimnames(mns) <- list(nams, 1:dims)  
  print(round(mns, digits))
}

summary.mixDist <- function(object, digits = 4, ...){
  class(object) <- "summary.mixDist"
  object
}

print.summary.mixDist <- function(x, digits = 4, ...){
  dims <- dim(x$sigmas[[1]])[1]
  ncomp <- length(x$weights)
  cat("Mixture Approximation\n\n")
  cat("Number of Components:", ncomp, "\n")
  for(m in 1:ncomp){
    cat("------------------------------\n")    
    cat("Component", m, "\n")
    cat("Weight:", signif(x$weights[m], digits), "\n")
    cat("Means and Standard Variations:\n")
    mns <- x$means[m,]
    sdev <- sqrt(diag(x$sigmas[[m]]))
    pmat <- rbind(means=mns, sdevs = sdev)
    dimnames(pmat) <- list(c("means", "sdev"), 1:dims)
    print(round(pmat, digits))
    if(dims > 1){
      cat("Correlation Matrix:\n")
      corMat <- cov2cor(x$sigmas[[m]])
      dimnames(corMat) <- list(1:dims, 1:dims)
      corMat <- round(corMat, digits)
      cf <- format(corMat, digits = digits)
      cf[row(cf) > col(cf)] <- ""
      print(cf, quote = FALSE)
    }
    cat("\n")
  }
}


getGrid <- function (n, mean, eigenSigma, hess = TRUE,
                     method = c("Pseudo-Rand", "Quasi-Rand"),
                     df = Inf, ctrl = list(seed=4711, scrambling=1)){
    if (length(mean) != nrow(eigenSigma$vectors)) {
        stop("mean and sigma have non-conforming size")
    }
    if(hess){
      ord <- order(eigenSigma$values)
      eg <- eigenSigma
      eigenSigma$values <- 1/eg$values[ord]
      eigenSigma$vectors <- eg$vectors[,ord]
    }
    dim <- length(mean)
    method <- match.arg(method)
    if(method == "Pseudo-Rand"){
      if(df == Inf){
        vals <- rnorm(n * dim)
      } else {
        vals <- rt(n * dim, df = df)
      }
      Mat <- matrix(vals, nrow = n)
    }
    if(method == "Quasi-Rand"){
      if(method == "Quasi-Rand")
        Mat <- sobol(n, dim = dim, 
                     scrambling = ctrl$scrambling,
                     seed = ctrl$seed)
      if(df == Inf){      
        Mat <- qnorm(Mat)
      } else {
        Mat <- qt(Mat, df = df)
      }
    }
    ev <- eigenSigma
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors)
    retval <- Mat %*% retval
    retval <- sweep(retval, 2, mean, "+")
    retval
}

dmvtAdapted <- function (x, obj, ind = NULL, mat = FALSE, df = 4){
  if (is.vector(x)) {
    x <- matrix(x, nrow = length(x))
  }
  dm <- length(obj$sigmas)
  if(is.null(ind)){
    ind <- 1:dm
  } else {
    if(!mat)
      stop("Cannot calculate density.")
  }
  out <- matrix(nrow = nrow(x), ncol = length(ind))
  z <- 1
  for(i in ind){
    sigmaInv <- obj$sigmainv[[i]]
    detSigma <- obj$dets[i]
    mean <- obj$means[i,]
    m <- NCOL(sigmaInv)
    centr <- sweep(x, 2, mean)
    distval <- rowSums((centr %*% sigmaInv) * centr)
    logdet <- log(detSigma)
    logretval <- lgamma((m + df)/2) - (lgamma(df/2) + 0.5 * (logdet + 
        m * logb(pi * df))) - 0.5 * (df + m) * logb(1 + distval/df)
    out[,z] <- exp(logretval)
    z <- z + 1
  }
  if(mat){
    return(out)
  } else {
    out%*%obj$weights
  }
}

dmvtAdaptedC <- function(x, obj, df){
  mns <- as.double(t(obj$means))
  sigmainv <- as.double(unlist(obj$sigmainv))
  dets <- as.double(obj$dets)
  weights <- as.double(obj$weights)
  ncomp <- as.integer(length(weights))
  dm <- as.integer(length(mns)/ncomp)
  df <- as.integer(df)

  out <- .C("calcmixtdens", as.double(x), dm, ncomp, mns,
            sigmainv, dets, weights, df, double(dm), double(1))[[10]]
  out  
}

dmvnormAdapted <- function (x, obj, ind = NULL, mat = FALSE){
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  dm <- length(obj$sigmas)
  if(is.null(ind)){
    ind <- 1:dm
  } else {
    if(!mat)
      stop("Cannot calculate density.")
  }
  out <- matrix(nrow = nrow(x), ncol = length(ind))
  z <- 1
  for(i in ind){
    sigmaInv <- obj$sigmainv[[i]]
    detSigma <- obj$dets[i]
    mean <- obj$means[i,]
    m <- NCOL(sigmaInv)
    centr <- sweep(x, 2, mean)
    distval <- rowSums((centr %*% sigmaInv) * centr)
    logdet <- log(detSigma)
    logretval <- -m/2*log(2*pi)-0.5*logdet-0.5*distval
    out[,z] <- exp(logretval)
    z <- z + 1
  }
  if(mat){
    return(out)
  } else {
    out%*%obj$weights
  }
}

dmvnormAdaptedC <- function(x, obj){
  mns <- as.double(t(obj$means))
  sigmainv <- as.double(unlist(obj$sigmainv))
  dets <- as.double(obj$dets)
  weights <- as.double(obj$weights)
  ncomp <- as.integer(length(weights))
  dm <- as.integer(length(mns)/ncomp)

  out <- .C("calcmixmvdens", as.double(x), dm, ncomp, mns,
            sigmainv, dets, weights, double(dm), double(1))[[9]]
  out  
}


LSfitgeq0 <- function(y, X){
  D <- crossprod(X)
  d <- crossprod(y, X)
  Amat <- diag(dim(X)[2])
  bvec <- numeric(dim(X)[2])
  out <- try(solve.QP(D, d, Amat, bvec)$solution, silent = TRUE)
  if(!inherits(out, "try-error")){
    return(out)
  } else { # if solve.QP does not work use heuristic approach
    cf <- qr.coef(qr(X), y)
    ind <- cf < 0
    if(any(ind)){
      cf <- qr.coef(qr(X[,!ind]), y)
      out <- numeric(ncol(X))
      out[!ind] <- cf
      out[out < 0] <- 0
    }
    return(out)
  }
}

samplingMixDist <- function(nSim, obj, df = 4){
  ## nSim - number of simulations
  ## obj - mixDist object
  if(!inherits(obj, "mixDist"))
    stop("obj needs to be of class mixDist")
  weights <- obj$weights
  means <- obj$means
  egHess <- obj$eigenHess
  k <- length(weights)
  kmat <- rmultinom(1, nSim, weights)
  samples <- NULL
  for(i in 1:k){
    if(kmat[i,1] > 0){
      mat0 <- rmvtAdapt(kmat[i,1], means[i,], egHess[[i]], df = df)
      samples <- rbind(samples, mat0)
    }
  }
  ind <- sample(1:nSim)
  samples[ind,]
}

IS <- function(obj, nSim, df = 4, post, vectorized = FALSE,
               cores = 1, ...){
  if(!inherits(obj, "mixDist"))
    stop("obj needs to be of class mixDist")
  samp <- samplingMixDist(nSim, obj, df = df)
  if(vectorized){
    vals1 <- post(samp,...)
  } else {
    sampList <- split(samp, 1:nSim)
    if(cores == 1){
      vals1 <- lapply(sampList, function(x) post(x, ...))
    } else {
      vals1 <- mclapply(sampList, function(x) post(x, ...), mc.cores = cores)
    }
    vals1 <- unlist(vals1)
  }
  scal <- max(vals1)
  vals1 <- exp(vals1 - scal)
  if(df == Inf){
    vals2 <- dmvnormAdapted(samp, obj)    
  } else {
    vals2 <- dmvtAdapted(samp, obj, df = df)
  }
  w <- vals1/vals2
  const <- sum(w)
  w <- w/const
  out <- list(samp = samp, w = w, normconst = const/nSim*exp(scal),
              ESS = 1/sum(w^2))
  class(out) <- "IS"
  out
}

IMH <- function(obj, nSim, df = 4, post, vectorized = FALSE, cores = 1, ...){
  if(!inherits(obj, "mixDist"))
    stop("obj needs to be of class mixDist")
  samp <- samplingMixDist(nSim, obj, df = df)
  if(is.vector(samp)){
    samp <- matrix(samp, nrow = length(samp))
  }
  if(vectorized){
    vals1 <- post(samp,...)
  } else {
    sampList <- split(samp, 1:nSim)#as.list(data.frame(t(samp)))
    if(cores == 1){
      vals1 <- lapply(sampList, function(x) post(x, ...))
    } else {
      #vals1 <- mclapply(sampList, function(x) post(x, ...), mc.cores = cores)
      stop("currently only cores = 1 possible")      
    }
    vals1 <- unlist(vals1)
  }
  scal <- max(vals1)
  vals1 <- exp(vals1 - scal)
  vals2 <- dmvtAdapted(samp, obj, df = df)
  w <- vals1/vals2
  sampOut <- sampOut <- matrix(0, nrow = nSim, ncol = ncol(samp))
  wcur <- w[1]
  sampCur <- sampOut[1,] <- samp[1,]
  u <- runif(nSim)
  z <- 1
  for(i in 2:nSim){
    if(u[i] < w[i]/wcur){
      sampOut[i,] <- sampCur <- samp[i,]
      wcur <- w[i]
      z <- z + 1
    } else {
      sampOut[i,] <- sampCur
    }
  }
  const <- sum(w)
  w <- w/const
  out <- list(samp = sampOut, w = w, normconst = const/nSim*exp(scal),
              accept = z/nSim)
  class(out) <- "IMH"
  out
}

rmvtAdapt <- function(n, mean = NULL, eigenSigma = NULL, df){
  ## random numbers from multivariate t-distribution
  dm <- length(eigenSigma$values)
  grd <- getGrid(n, mean = rep(0, dm), eigenSigma, hess = TRUE,
                 method = "Pseudo-Rand")
  if(df != Inf){
    sqchsq <- sqrt(rchisq(n, df)/df)
  } else {
    sqchsq <- rep(1,n)
  }
  grd <- grd/sqchsq 
  mnmat <- matrix(mean, byrow = T, nrow = n, ncol = dm)
  grd + mnmat
}

resample <- function(n, obj){
  if(!inherits(obj, "IS"))
    stop("obj needs to be of class IS")
  samples <- obj$samp
  wgts <- obj$w
  ## performs residual resampling
  if(is.vector(samples)){
    samples <- matrix(samples, nrow = length(samples))
  }
  nS <- nrow(samples)
  reps <- floor(wgts*n)
  if(!all(reps==0)){
    inds <- rep(1:nS, reps)
    probAdjusted <- (wgts*n-reps)/(n-sum(reps))
    inds2 <- sample(1:nS, n-sum(reps), prob=probAdjusted, replace=TRUE)
    inds <- c(inds, inds2)
  } else {
    inds <- sample(1:nS, n, prob=wgts, replace=TRUE)
  }
  samples[inds,]
}

iterLap <- function(post, ..., GRobj = NULL, 
                 vectorized = FALSE, startVals = NULL,
                 method = c("nlminb", "nlm", "Nelder-Mead", "BFGS"),
                 control = NULL, nlcontrol = list()){
  ## GRobj - object of class mixDist
  ## post - log-posterior
  ## control - control parameters for iterated Laplace approximation
  ## nlcontrol - control parameters for nlm

  method <- match.arg(method)
  ## extract dimension and convert startVals to matrix
  if(!is.null(startVals)){
    if(!is.matrix(startVals)){
      startVals <- as.matrix(startVals, ncol=1)
    }
    pardim <- ncol(startVals)
  }

  ## calculcate (or extract) GR object
  if(is.null(GRobj)){
    if(is.null(startVals)){
      stop("need to specify starting matrix if GRobj is missing")
    }
    GRobj <- GRApprox(post, startVals, method = method, control=nlcontrol, ...)
  }
  obj <- GRobj # obj contains current approximation
  pardim <- ncol(obj$means)

  ## check control list
  if (!is.null(control)) {
    if(!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("iterLap.control", c(pardim=pardim,control))
  } else {
    ctrl <- iterLap.control(pardim = pardim)
  }
  
  propnorm <- numeric(ctrl$maxDim)
  
  ## calculate starting grid (from GR approximation)
  dm <- length(obj$weights)  
  grid <- NULL
  for(i in 1:dm){
    grd <- getGrid(ctrl$gridSize, obj$means[i,], obj$eigenHess[[i]], hess = TRUE,
                   method = "Quasi-Rand", df=Inf)
    grid <- rbind(grid, grd)
  }

  ## evaluate kernel and mixture of t-distributions at starting grid
  if(vectorized){
    vals <- post(grid, ...)
  } else {
    vals <- apply(grid, 1, function(x) post(x, ...))
  }
  lgMax <- max(vals)
  X <- dmvnormAdapted(grid, obj, mat = TRUE)
  z <- dm
  
  ## determine weights
  obj$weights <- LSfitgeq0(exp(vals-lgMax), X)
  propnorm[z] <- sum(obj$weights)*exp(lgMax)

  ## iterate number of components
  while(z < ctrl$maxDim){
    if(ctrl$info > 0){
      cat("Add Component: ", z+1, "\n")
    }

    pred <- X%*%obj$weights   # calculate predictions 
    if(checkStop(vals, lgMax, pred, ctrl, propnorm, z)){ # check stop criterion
      break
    } else { # calculate starting values for optimizer of residual
      startOpt <- getStartVal(vals, lgMax, pred, grid, obj$means[z,], z)
    }

    ## run optimizer
    optres <- opt(obj, z, startOpt, method, ctrl, nlcontrol, post, lgMax, ...)
    if(!optres$add){ # optimizer output does not contain local optimum
      if(ctrl$info > 0){
        cat("Cannot improve approximation.\n")
      }
      break
    } else {
      ## add component to approximation
      obj <- optres$obj      
      z <- z + 1

      ## evaluate new component on current grid (ind = z)
      Xgrd <- dmvnormAdapted(grid, obj, ind = z, mat = TRUE)
      X <- cbind(X, Xgrd)
      ## calculate grid of new component
      grdNew <- getGrid(ctrl$gridSize, obj$means[z,], obj$eigenHess[[z]],
                        hess = TRUE, method = "Quasi-Rand", df = Inf)
      ## evaluate all components on new grid
      Xnew2 <- dmvnormAdapted(grdNew, obj, mat = TRUE)
      X <- rbind(X, Xnew2)
      grid <- rbind(grid, grdNew)
      ## evaluate target function at new component
      if(vectorized){
        valsnew <- post(grdNew, ...)
      } else {
        valsnew <- apply(grdNew, 1, function(x) post(x, ...)) 
      }
      if(any(is.na(valsnew)) | any(abs(valsnew)==Inf)){
        warning("Inf/NA produced.")
      }
      vals <- c(vals, valsnew)
      lgMax <- max(vals)      
      obj$weights <-  LSfitgeq0(exp(vals-lgMax), X)
      propnorm[z] <- sum(obj$weights)*exp(lgMax)
      if(ctrl$info > 0){
        cat("Norm. Constant of Approx.: ", signif(propnorm[z], 5), "\n\n")
      }
    }
  }
  obj$weights <- abs(obj$weights/sum(abs(obj$weights)))
  class(obj) <- "mixDist"
  obj
}

opt <- function(obj, z, start, method, control, nlcontrol, post, lgMax, ...){
  ## optimizer
  post2 <- function(pars, obj, eps = 0.0001, ...){
    val <- exp(post(pars, ...)-lgMax)-dmvnormAdaptedC(pars, obj)
    if(val < eps | is.na(val)){
      val <- exp(val-eps)*eps
    }
    return(-log(val))
  }

  j <- 0
  while(j < nrow(start)){
    j <- j + 1
    nlcall <- list(par=start[j,], fn = post2, ...,
                   method = method, obj = obj, control = nlcontrol)
    optobj <- do.call(nlopt, nlcall)
    eg <- eigen(optobj$hessian, symmetric = TRUE)
    detSig <- 1/prod(eg$values)
    if(all(eg$values > 0)){ # only include maxima of posterior
      m <- optobj$par
      sw <- sweep(obj$means, 2, m) # check whether found optimum coincides with already found optimum
      mabs <- sum(abs(m))
      an <- apply(sw, 1, function(x) sum(abs(x)))
      ind2 <- all(an[1:z] > 0.005*mabs)
      if(ind2){
        if(control$info > 1){
          cat("Added Component:", signif(m,5), "\n")
        }
        obj$means <- rbind(obj$means, m)
        obj$sigmas[[z+1]] <- solve(optobj$hessian)
        obj$sigmainv[[z+1]] <- optobj$hessian
        obj$eigenHess[[z+1]] <- eg
        obj$dets <- c(obj$dets, detSig)
        add <- TRUE
        break
      } else {
        if(control$info > 1){        
          cat("Too close to already used components\n")
        }
        add <- FALSE
      }
    } else {
      if(control$info > 1){
        cat("Negative Eigenvalues\n")
      }
      add <- FALSE
    }
  }
  list(obj = obj, add = add)
}

calcDist <- function(start, oldVal){
  dfs <- sweep(start, 2, oldVal)
  rev(order(rowSums(dfs^2)))
}

extractMoments <- function(obj){
  mn <- obj$weights%*%obj$means
  sigm <- mn2 <- matrix(0, nrow = ncol(obj$means), ncol = ncol(obj$means))
  for(j in 1:length(obj$weights)){
    sigm <- sigm + obj$weights[j]*obj$sigmas[[j]]
    mn2 <- mn2 + obj$weights[j]*obj$means[j,]%*%t(obj$means[j,])
  }
  sigm <- sigm + mn2 - t(mn)%*%mn
  list(mean = as.numeric(mn), sigma = sigm)
}

checkStop <- function(vals, scal, pred, control, propnorm, z){
  out <- FALSE
  if(max(exp(vals-scal) - pred) < control$delta){
    out <- TRUE
    if(control$info > 0){
      cat("Maximum Error on grid smaller than delta\n")
    }
  }
  if(z >= 3){
    dif <- abs(propnorm[z]-mean(propnorm[(z-1):(z-2)]))
    std <- propnorm[z]
    if(dif/std < control$eps){
      if(control$info > 0){
        cat("Approximation of normalization constant does not improve\n")
      }
      out <- TRUE
    }
  }
  out
}

getStartVal <- function(vals, scal, pred, grid, oldVal, init){
  ## determine 10 grid point with largest importance weight
  ind <- rev(order(exp(vals-scal) / pred))[1:10]
  ## select starting values for k-means
  vals <- c(1,9,3,7,2,8,4,6,5,6,1,5,8,4,7,3,9,2,4,6,5,2,8,3,7,1,9)
  val <- vals[(init-1)%%27+1]        # emulate some randomness for ind0
  ind0 <- c(val, val+3, val+6)%%10+1 # this increases diversity
  startkm <- grid[ind[ind0],]        
  ## reduce this to 3 possibly distinct values
  start <- kmeans(grid[ind,], startkm)$centers
  ## sort those 3 values according to their distance to last point
  ind2 <- calcDist(start, oldVal)
  start[ind2, , drop=FALSE]
}

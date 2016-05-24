######################################
# Finds the cross-validated loglikelihood for a given penalty
######################################
cvl <- function(response, penalized, unpenalized, lambda1 = 0, lambda2= 0, positive = FALSE, 
    fusedl = FALSE, data, model = c("cox", "logistic", "linear", "poisson"), 
    startbeta, startgamma, fold, epsilon = 1e-10, maxiter, standardize = FALSE, 
    trace = TRUE, approximate = FALSE) {

  # Maximum number of iterations depends on the input
  if (missing(maxiter)) maxiter <- if (lambda1 == 0 && lambda2 == 0 && !positive) 25 else Inf

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # check for the presence of penalty parameters
  if (ncol(prep$X) >= nrow(prep$X) && all(lambda1 == 0) && all(lambda2 == 0))
    stop("High-dimensional data require a penalized model. Please supply lambda1 or lambda2.", call.=FALSE)

  # prepare the model
  fit <- .modelswitch(prep$model, prep$response, prep$offset, prep$strata)
  
  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  chr = prep$chr
  fusedl = prep$fusedl
  
  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)
                                 
  # choose the fastest way of computing in case of a linear model
  if(!fusedl){
  if(prep$model == "linear" && fold>1 && all(lambda1==0) && !any(positive)) approximate = TRUE }
  
  # make vectors of lambda1 and lambda2
  if (length(lambda1) == 1)
    uselambda1 <- lambda1 * prep$baselambda1
  else
    uselambda1 <- c(numeric(pu), lambda1) * prep$baselambda1
  if (length(lambda2) == 1)
    uselambda2 <- lambda2 * prep$baselambda2
  else
    uselambda2 <- c(numeric(pu), lambda2) * prep$baselambda2

  if(approximate && !fusedl)
  {
    # check if lambda1 = 0
    if (!all(lambda1 == 0))
    stop("Approximation method only works for ridge penalty, so lambda1 can not differ from 0", call.=FALSE)    
    # check if fold > 1
    if(fold <= 1)
    stop("Approximation method is not implemented for 1-fold cross-validation", call.=FALSE)
    # check if positive = FALSE everywhere
    if(any(positive))
    stop("Approximation method can not be used in combination with positivity constraints", call.=FALSE)     
            
    res <- .cvlapprox(prep$X, uselambda1, uselambda2,
    positive = prep$positive, beta = prep$beta, fit=fit, groups=groups, 
    epsilon=epsilon, maxiter=maxiter, trace = trace, quit.if.failed = FALSE)
    if(is.na(res$cvl) || res$cvl==-Inf)   #==? 
    {
      res$predictions <- NA
    }
    else
    {
      res$predictions <- .predictswitch(prep$model, res$predictions, groups)
    }
  }
  else
  {
    res <- .cvl(prep$X, uselambda1, uselambda2, fusedl = fusedl, chr = chr,
    positive = prep$positive, beta = prep$beta, fit=fit, groups=groups, 
    epsilon=epsilon, maxiter=maxiter, trace = trace, quit.if.failed = FALSE)
        if(is.na(res$cvl))
    {
      res$predictions <- NA
    }
    else
    {
      res$predictions <- .predictswitch(prep$model, res$predictions, groups)
    }
  }

  return(list(
    cvl = res$cvl,
    cvls = res$cvls,
    predictions = res$predictions, 
    fold = groups, 
    fullfit = .makepenfit(res$fit, pu, prep$model, lambda1, fusedl = fusedl, 
      lambda2, prep$orthogonalizer, prep$weights, prep$formula)
  ))
}


######################################
# Finds the curve of the cross-validated likelihood for a given L2-penalty and a range of L1-penalty values
######################################
profL1 <- function(response, penalized, unpenalized, minlambda1, maxlambda1, base1, lambda2 = 0,
  fusedl = FALSE,positive = FALSE, data, model = c("cox", "logistic", "linear", "poisson"), 
  startbeta, startgamma, fold, epsilon = 1e-10, maxiter = Inf, standardize = FALSE, steps = 100, 
  minsteps = steps/4, log = FALSE, save.predictions = FALSE, trace = TRUE, plot = FALSE) {

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # prepare the model
  fit <- .modelswitch(prep$model, prep$response, prep$offset, prep$strata)

  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  chr <- prep$chr
  fusedl = prep$fusedl

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)
  
  #check if fold=1 and save.predictions=TRUE at the same time
  if(fold <= 1 && save.predictions)
  {
    warning("save.predictions is set to FALSE, because the number of folds equals 1")
    save.predictions <- FALSE
  }

  # make vectors of lambda1 and lambda2

  if (missing(base1)) {
    baselambda1 <- prep$baselambda1
    base1 <- 1
  } else
    baselambda1 <- c(numeric(pu), base1) * prep$baselambda1
  if (length(lambda2) == 1)
    uselambda2 <- lambda2 * prep$baselambda2
  else
    uselambda2 <- c(numeric(pu), lambda2) * prep$baselambda2    

   
 
  # benchmark: cvl at infinite penalty
  if (pu > 0) {
    nullfit <- .cvl(prep$X[,1:pu, drop=FALSE], lambda1 = rep(0,pu), chr = chr[1:pu], 
               fusedl = prep$fusedl,lambda2 = rep(0,pu), positive = FALSE, beta = prep$nullgamma, 
               fit=fit, groups=groups, epsilon=epsilon, maxiter=maxiter, trace = FALSE, save.predictions = FALSE)
    nullgammas <- nullfit$betas
    nullcvl <- nullfit$cvl
    nullfit <- nullfit$fit
  } else {
    nullcvl <- fit$cvl(numeric(n), !logical(n))
    nullfit <- list()
    nullfit$fit <- fit$fit(numeric(n))
    nullfit$iterations <- 1
    nullfit$converged <- TRUE
  }
  
  # find the maxlambda1 and minlambda1
  if (missing(maxlambda1)) {
    if (pu > 0) 
      lps <- drop(prep$X[,1:pu,drop=FALSE] %*% nullgammas)
    else 
      lps <- matrix(0, n, fold)
    gradients <- sapply(1:fold, function(ff) {
      drop(crossprod(prep$X[groups!=ff,pu+1:pp,drop=FALSE], fit$fit(lps[groups!=ff,ff], groups==ff)$residuals))
    })
    rel <- gradients / matrix(baselambda1[pu+1:pp], pp, ncol(lps))
    maxlambda1 <- max(apply(rel, 2, function (ff)
      max(ifelse(prep$positive[pu+1:pp], ff, abs(ff)))
    ))
  }  
  if (missing(minlambda1)) {
    if (log) 
      stop("argument \"minlambda1\" is missing. please specify \"minlambda1\" or set log = FALSE", call. = FALSE)
    else
      minlambda1 <- 0
  }  
  
  # find the sequence from maxlambda1 to minlambda1
  if (steps < 2) stop("please set \"steps\" >= 2", call. = FALSE)
  if (log) {
    lambda1s <- exp(seq(log(maxlambda1), log(minlambda1), length.out = steps))
  } else {
    lambda1s <- seq(maxlambda1, minlambda1, length.out = steps)
  }
  if (pp+pu >= n) # HD model
    lambda1s <- lambda1s[lambda1s != 0]

  # the actual repeated cvl-calculation
  betas <- NULL
  beta <- prep$beta
  cvls <- rep(NA,length(lambda1s))
  finished <- FALSE
  iter <- 0
  fits <- vector("list", length = length(lambda1s))
  predictions <- vector("list", length = length(lambda1s))
  while (!finished) {
    iter <- iter + 1
    rellambda <- lambda1s[iter]
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(prep$X, rellambda*baselambda1, uselambda2, chr = chr,
      fusedl=fusedl,positive = prep$positive, 
      beta = beta, fit=fit, groups=groups,epsilon=epsilon, maxiter=maxiter, 
      trace = trace, betas = betas,quit.if.failed=FALSE, save.predictions = save.predictions)
    if (trace) cat("cvl=", out$cvl, "\n")
    beta <- out$fit$beta
    betas <- out$betas
    cvls[iter] <- out$cvl
    fits[[iter]] <- out$fit
    if (save.predictions)
      predictions[[iter]] <- out$predictions
    finished <- ((fold > 1) && ((cvls[[iter]] < min(c(nullcvl, cvls[1:(iter-1)]))) && (iter >= minsteps))) || (iter == length(lambda1s))
  }

  # remove the tail of the output
  if (fold > 1) {
    lambda1s <- lambda1s[!is.na(cvls)]
    fits <- fits[!is.na(cvls)]                                           
    predictions <- predictions[!is.na(cvls)]
    cvls <- cvls[!is.na(cvls)]
  }

  # merge the cross-validated predictions
  if (save.predictions)
    predictions <- lapply(predictions, function(preds) .predictswitch(prep$model, preds, groups))
  
  # create all the penfit objects
  makethisfit <- function(iter) {
    .makepenfit(fits[[iter]], pu, fusedl = prep$fusedl, prep$model, lambda1s[[iter]]*base1, lambda2,
    prep$orthogonalizer, prep$weights, prep$formula)
  }

  if (plot && !fusedl)
    plot(lambda1s, cvls, type="l", log="x", ylab="cvl", xlab="lambda")
   if (plot && fusedl) 
     stop("not for fused lasso. please set \"fusedl\" = FALSE", call. = FALSE)
  return(list(
    lambda = lambda1s, 
    fold = groups, 
    cvl = cvls, 
    predictions = if (save.predictions) predictions else "not saved", 
    fullfit = lapply(1:length(cvls), makethisfit)
  ))
} 

######################################
# Finds the curve of the cross-validated likelihood for a given L2-penalty and a range of L1-penalty values
######################################
profL2 <- function(response, penalized, unpenalized, lambda1 = 0, minlambda2, maxlambda2, base2, 
 fusedl = FALSE,positive = FALSE, data, model = c("cox", "logistic", "linear", "poisson"), 
  startbeta, startgamma, fold, epsilon = 1e-10, maxiter, standardize = FALSE, steps = 100, minsteps = steps/2, 
  log = TRUE, save.predictions = FALSE, trace = TRUE, plot = FALSE, approximate = FALSE) {

  # Maximum number of iterations depends on the input
  if (missing(maxiter)) maxiter <- if (all(lambda1 == 0) && !any(positive)) 25 else Inf

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # prepare the model
  fit <- .modelswitch(prep$model, prep$response, prep$offset, prep$strata)
  
  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  chr <- prep$chr
  fusedl = prep$fusedl

  # make vectors of lambda1 and lambda2
  if (length(lambda1) == 1)
    uselambda1 <- lambda1 * prep$baselambda1
  else
    uselambda1 <- c(numeric(pu), lambda1) * prep$baselambda1
  if (missing(base2)) {
    base2 <- 1
    baselambda2 <- prep$baselambda2
  } else
    baselambda2 <- c(numeric(pu), base2) * prep$baselambda2

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)
  
  # choose the fastest way of computing in case of a linear model
  if(!fusedl){
  if(prep$model == "linear" && fold>1 && all(lambda1==0) && !any(positive)) approximate = TRUE }   
  
  #check if fold=1 and save.predictions=TRUE at the same time
  if(fold <= 1 && save.predictions)
  {
    warning("save.predictions is set to FALSE, because the number of folds equals 1")
    save.predictions <- FALSE
  }

  # Find the sequence from maxlambda2 to minlambda2
  if (missing(minlambda2)) 
    if (!log) minlambda2 <- 0 else stop("Agrument \"minlambda2\" is missing with no default")
  if (missing(minlambda2)) stop("Agrument \"maxlambda2\" is missing with no default")
  if (steps < 2) stop("please set \"steps\" >= 2", call. = FALSE)
  if (log) 
    lambda2s <- exp(seq(log(maxlambda2), log(minlambda2), length.out = steps))
  else
    lambda2s <- seq(maxlambda2, minlambda2, length.out = steps)

  # benchmark: cvl at infinite penalty
  if (pu > 0) {
    nullfit <- .cvl(prep$X[,1:pu, drop=FALSE], lambda1 = rep(0,pu), lambda2 = rep(0,pu),  
      chr = chr[1:pu], fusedl =prep$fusedl, positive = FALSE,
      beta = prep$nullgamma, fit=fit, groups=groups, epsilon=epsilon, 
      maxiter=maxiter, trace = FALSE, save.predictions = FALSE)
    nullcvl <- nullfit$cvl
    nullfit <- nullfit$fit
  } else {
    nullcvl <- fit$cvl(numeric(n), !logical(n))
    nullfit <- list()
    nullfit$fit <- fit$fit(numeric(n))
    nullfit$iterations <- 1
    nullfit$converged <- TRUE
  }

  # the actual repeated cvl-calculation
  betas <- NULL
  beta <- prep$beta
  cvls <- rep(NA,length(lambda2s))
  finished <- FALSE
  iter <- 0
  fits <- vector("list", length = length(lambda2s))
  predictions <- vector("list", length = length(lambda2s))
  while (!finished) {
    iter <- iter + 1
    rellambda <- lambda2s[iter]
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    if(approximate && !fusedl)
    {
      # check if lambda1 = 0
      if (!all(lambda1 == 0))
      stop("Approximation method only works for ridge penalty, so lambda1 can not differ from 0", call.=FALSE)
      # check if fold > 1
      if(fold <= 1)
      stop("Approximation method is not implemented for 1-fold cross-validation", call.=FALSE)
      # check if positive = FALSE everywhere
      if(any(positive))
      stop("Approximation method can not be used in combination with positivity constraints", call.=FALSE)            
    
      out <- .cvlapprox(prep$X, uselambda1, rellambda*baselambda2,
        positive = prep$positive, beta = beta, fit=fit, groups=groups, 
        epsilon=epsilon, maxiter=maxiter, trace = trace, betas = betas, 
        quit.if.failed=FALSE, save.predictions = save.predictions)
    }
    else
    {
      out <- .cvl(prep$X, uselambda1, rellambda*baselambda2, chr = chr,
        fusedl =prep$fusedl,positive = prep$positive, 
        beta = beta, fit=fit, groups=groups, epsilon=epsilon, maxiter=maxiter,
        trace = trace, betas = betas, quit.if.failed=FALSE, save.predictions = save.predictions)
    }    
    if (trace) if (fold > 1) cat("cvl=", out$cvl, "\n") else cat("\n")
    beta <- out$fit$beta
    betas <- out$betas
    cvls[iter] <- out$cvl
    fits[[iter]] <- out$fit
    if (save.predictions) predictions[[iter]] <- out$predictions
    finished <- ((fold > 1) && (cvls[[iter]] < min(c(nullcvl, cvls[1:(iter-1)]))) &&
      (iter >= minsteps)) || (iter == length(lambda2s))
  }

  # remove the tail of the output
  if (fold > 1) {
    lambda2s <- lambda2s[!is.na(cvls)]
    fits <- fits[!is.na(cvls)]
    predictions <- predictions[!is.na(cvls)]
    cvls <- cvls[!is.na(cvls)]
  }
         
  # merge the cross-validated predictions 
  if (save.predictions)               
    predictions <- lapply(predictions, function(preds) .predictswitch(prep$model, preds, groups))
    
  # create the penfit objects                          
  makethisfit <- function(iter) {
    .makepenfit(fits[[iter]], pu, prep$model, lambda1, lambda2s[[iter]]*base2, fusedl = prep$fusedl,
      prep$orthogonalizer, prep$weights, prep$formula)
  }

  if (plot && !fusedl)
    plot(lambda2s, cvls, type="l", log="x", ylab="cvl", xlab="lambda")
    
   if (plot && fusedl) 
     stop("not for fused lasso. please set \"fusedl\" = FALSE", call. = FALSE)
  return(list(
    lambda = lambda2s, 
    fold = groups, 
    cvl = cvls, 
    predictions = if (save.predictions) predictions else "not saved", 
    fullfit = lapply(1:length(cvls), makethisfit)
  ))
} 

######################################
# Finds the optimal cross-validated L1-penalty for a given L2-penalty
######################################
optL1 <- function(response, penalized, unpenalized, minlambda1, maxlambda1, base1, lambda2 = 0, 
  fusedl = FALSE, positive = FALSE, data, model = c("cox", "logistic", "linear", "poisson"), 
  startbeta, startgamma, fold, epsilon = 1e-10, maxiter = Inf, standardize = FALSE, tol = .Machine$double.eps^0.25, 
  trace = TRUE) {

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # prepare the model
  fit <- .modelswitch(prep$model, prep$response, prep$offset, prep$strata)

  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  chr <- prep$chr
  fusedl = prep$fusedl
  # make vectors of lambda1 and lambda2
  if (missing(base1)) {
    baselambda1 <- prep$baselambda1
    base1 <- 1
  } else
    baselambda1 <- c(numeric(pu), base1) * prep$baselambda1
  if (length(lambda2) == 1)
    uselambda2 <- lambda2 * prep$baselambda2
  else
    uselambda2 <- c(numeric(pu), lambda2) * prep$baselambda2

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)
  
  # check if fold > 1
  if(fold <= 1)
  stop("Method is not implemented for 1-fold cross-validation", call.=FALSE)  

  # find the maxlambda1 and minlambda1
  if (missing(maxlambda1)) {
    if (pu > 0) {
      nullfit <- .cvl(prep$X[,1:pu, drop=FALSE], lambda1 = rep(0,pu), chr = chr[1:pu],
        fusedl =prep$fusedl,lambda2 = rep(0,pu), positive = FALSE, 
        beta = prep$nullgamma, fit=fit, groups=groups, epsilon=epsilon, maxiter=maxiter, trace = FALSE)
      nullgammas <- nullfit$betas
      lps <- drop(prep$X[,1:pu,drop=FALSE] %*% nullgammas)
    } else {
      lps <- matrix(0, n, fold)
    }
    gradients <- sapply(1:fold, function(ff) {
      drop(crossprod(prep$X[groups!=ff,pu+1:pp,drop=FALSE], fit$fit(lps[groups!=ff,ff], groups==ff)$residuals))
    })
    rel <- gradients / matrix(baselambda1[pu+1:pp], pp, ncol(lps))
    maxlambda1 <- max(apply(rel, 2, function (ff)
      max(ifelse(prep$positive[pu+1:pp], ff, abs(ff)))
    ))
  }  
  if (missing(minlambda1)) minlambda1 <- 0

  # The function to be optimized
  # Note: in passing it keeps track of the fit with the best cvl so far
  # It does this to be able to return more than just the optimal cvl-value
  betas <- NULL
  beta <- prep$beta
  maxcvl <- -Inf
  best <- NULL
  thiscvl <- function(rellambda) {
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    out <- .cvl(prep$X, rellambda*baselambda1, uselambda2, chr = chr, 
      fusedl =prep$fusedl,positive = prep$positive, 
      beta = beta, fit=fit, groups=groups, epsilon=epsilon, 
      maxiter=maxiter, trace = trace, betas = betas)
    if (trace) cat("cvl=", out$cvl, "\n")
    beta <<- out$fit$beta
    if (is.null(betas)) {
      between <- mean(abs(beta - out$fit$beta))
      within <- mean(abs(out$betas - matrix(out$fit$beta, pp+pu, fold)))
      if (between < within) {
        betas <<- out$betas
      } 
    } else {
      betas <<- out$betas
    }
    if (out$cvl > maxcvl) {
      maxcvl <<- out$cvl
      best <<- out
    }
    out$cvl
  }
  
  #optimize it
  opt <- opt.brent(thiscvl, c(minlambda1, maxlambda1), maximum = TRUE, tol = tol)
  
  # merge the cross-validated predictions of the optimal model
  best$predictions <- .predictswitch(prep$model, best$predictions, groups)

  return(list(
    lambda = opt$argmax, 
    cvl = opt$max, 
    predictions = best$predictions, 
    fold = groups, 
    fullfit = .makepenfit(best$fit, pu, prep$model, opt$argmax*base1, lambda2, fusedl =prep$fusedl,
      prep$orthogonalizer, prep$weights, prep$formula)
  ))
}

######################################
# Finds the optimal cross-validated L2-penalty for a given L1-penalty
######################################                       
optL2 <- function(response, penalized, unpenalized, lambda1 = 0, minlambda2, maxlambda2, base2, fusedl = FALSE ,positive = FALSE, data, 
  model = c("cox", "logistic", "linear", "poisson"), startbeta, startgamma, 
  fold, epsilon = 1e-10, maxiter, standardize = FALSE, tol = .Machine$double.eps^0.25, 
  trace = TRUE, approximate = FALSE) {
                                                                  
  # maximum number of iterations depends on the input
  if (missing(maxiter)) maxiter <- if (all(lambda1 == 0) && !any(positive)) 25 else Inf

  # call the general input checking function
  prep <- .checkinput(match.call(), parent.frame())

  # prepare the model
  fit <- .modelswitch(prep$model, prep$response, prep$offset, prep$strata)

  # retrieve the dimensions for convenience
  pu <- length(prep$nullgamma)
  pp <- ncol(prep$X) - pu
  n <- nrow(prep$X)
  chr <- prep$chr
  fusedl = prep$fusedl
  # make vectors of lambda1 and lambda2
  if (length(lambda1) == 1)
    uselambda1 <- lambda1 * prep$baselambda1
  else
    uselambda1 <- c(numeric(pu), lambda1) * prep$baselambda1
  if (missing(base2)) {
    base2 <- 1
    baselambda2 <- prep$baselambda2
  } else
    baselambda2 <- c(numeric(pu), base2) * prep$baselambda2

  # divide the samples into folds for cross-validation
  if (missing(fold)) fold <- n
  groups <- .getFolds(fold, n)
  names(groups) <- rownames(prep$X)
  fold <- max(groups)
  
  # check if fold > 1
  if(fold <= 1)
  stop("Method is not implemented for 1-fold cross-validation", call.=FALSE)
  
  # choose the fastest way of computing in case of a linear model
  if(prep$model == "linear" && fold>1 && all(lambda1==0) && !any(positive)) approximate = TRUE   

  # benchmark: cvl at infinite penalty
  thiscvlinf <- function() {
    if (pu > 0) {
      null <- .cvl(prep$X[,1:pu, drop=FALSE], lambda1 = rep(0,pu), lambda2 = rep(0,pu), 
        chr = chr[1:pu], fusedl =prep$fusedl, positive = FALSE, 
        beta = prep$nullgamma, fit=fit, groups=groups, epsilon=epsilon, 
        maxiter=maxiter, trace = FALSE)
      null$fit$beta <- c(null$fit$beta, numeric(pp))
    } else {
      null <- list()
      null.lp <- numeric(n)
      names(null.lp) <- rownames(prep$X)
      null$cvl <- do.call(sum, lapply(1:fold, function(ff) fit$cvl(numeric(n), groups == ff)))
      null$fit <- list()
      null$fit$beta <- c(numeric(pu+pp))
      null$fit$penalty <- c(L1 = 0, L2 = 0)
      names(null$fit$beta) <- colnames(prep$X)
      null$fit$fit <- fit$fit(null.lp)
      null$predictions <- lapply(as.list(1:n), function(i) fit$prediction(null.lp[i], nuisance= null$fit$fit$nuisance, which=i))
      null$fit$iterations <- 1
      null$fit$converged <- TRUE
    }
    null
  }
                 
  # The function to be optimized
  # Note: in passing it keeps track of the fit with the best cvl so far
  # It does this to be able to return more than just the optimal cvl-value
  betas <- NULL
  beta <- prep$beta
  best <- list(cvl=-Inf)
  thiscvl <- function(rellambda) {
    if (trace) {
      cat("lambda=", rellambda, "\t")
      flush.console()
    }
    if (rellambda == Inf) 
      out <- thiscvlinf()
    else {
       if(approximate && !fusedl)
       {
          # check if lambda1 = 0            
          if (!all(lambda1 == 0))
          stop("Approximation method only works for ridge penalty, so lambda1 can not differ from 0", call.=FALSE)
          # check if positive = FALSE everywhere
          if(any(positive))
          stop("Approximation method can not be used in combination with positivity constraints", call.=FALSE)               

            out <- .cvlapprox(prep$X, uselambda1, rellambda*baselambda2,
              positive = prep$positive, beta = beta, fit=fit, groups=groups, epsilon=epsilon, 
              maxiter=maxiter, trace = trace, betas = betas)                                                      
        }
        else
        {
            out <- .cvl(prep$X, uselambda1, rellambda*baselambda2, chr = chr, 
            fusedl =prep$fusedl,positive = prep$positive, 
            beta = beta, fit=fit, groups=groups, epsilon=epsilon, 
            maxiter=maxiter, trace = trace, betas = betas)
        }

      if (out$cvl > - Inf) {
        beta <<- out$fit$beta
        if (is.null(betas) && !is.null(out$fit$betas)) {
          between <- mean(abs(beta - out$fit$beta))
          within <- mean(abs(out$betas - matrix(out$fit$beta, pp+pu, fold)))
          if (between < within) {
            betas <<- out$betas
          } 
        } else {
          betas <<- out$betas
        }
      }
    }
    if (trace) cat("cvl=", out$cvl, "\n")
    if (out$cvl > best$cvl) {
      best <<- out
    }
    out$cvl
  }
  
  # phase 1: find the order of magnitude of lambda if not given by the user
  if (missing(maxlambda2))
    nullcvl <- thiscvl(Inf)
  else
    nullcvl <- thiscvl(maxlambda2)
  if (missing(minlambda2) || missing(maxlambda2)) {
    if (missing(minlambda2) && missing(maxlambda2)) {
      left <- 1
    } else if (missing(maxlambda2)) {
      left <- minlambda2
    } else {
      left <- maxlambda2 / 10
    }
    leftcvl <- thiscvl(left)
    right <- 10*left
    if (missing(maxlambda2))
      rightcvl <- thiscvl(right)
    else
      rightcvl <- nullcvl
    if (leftcvl < rightcvl || rightcvl == -Inf) {
      high <- right
      highcvl <- rightcvl
      low <- left 
      lowcvl <- leftcvl
      fac <- 10
    } else {
      high <- left
      highcvl <- leftcvl
      low <- right
      lowcvl <- rightcvl
      fac <- 0.1
    }
    ready <- (!missing(maxlambda2) && fac > 1) || (!missing(minlambda2) && fac < 1)
    # infmax: the maximum is (numerically) at infinite penalty
    # infmin: the maximum is (numerically) at zero penalty
    infmax <- FALSE
    infmin <- FALSE          
    if (!ready) {
      while (!ready && !infmax) {
        nxt <- high*fac
        nxtcvl <- thiscvl(nxt)
        ready <- nxtcvl < highcvl
        if (!ready) {
          low <- high
          lowcvl <- highcvl
          high <- nxt
          highcvl <- nxtcvl
        }
        infmax <- ((abs(lowcvl - nullcvl) / abs(nullcvl + 0.1) < epsilon) && 
          (abs(highcvl - nullcvl) / abs(nullcvl + 0.1) < epsilon))
        infmin <- (fac < 1) && (abs(lowcvl - nxtcvl) / abs(nxtcvl + 0.1) < epsilon) 
      }
      minlambda2 <- min(low, nxt)
      maxlambda2 <- max(low, nxt)
    } else {
      minlambda2 <- min(low, high)
      maxlambda2 <- max(low, high)
    }
  } else {
    infmax <- infmin <- FALSE
  }
                                
  # phase 2: optimize lambda within the order of magnitude found
  if (!infmax && !infmin) {
    opt <- opt.brent(thiscvl, sort(c(minlambda2,maxlambda2)), maximum = TRUE, 
      tol = tol)
    if (opt$max <= nullcvl) {
      opt <- list(argmax = maxlambda2, max = nullcvl)
    }
  } else {
    if (infmax) {
      opt <- list(argmax = Inf, max = nullcvl)
      names(best$fit$beta) <- colnames(prep$X)
    } else {
      best$cvl <- -Inf    # clears bestfit
      best$cvl <- thiscvl(0)
      opt <- list(argmax = 0, max = best$cvl) 
    } 
  }
                                
  # merge the cross-validated predictions at the optimal fit
  best$predictions <- .predictswitch(prep$model, best$predictions, groups)
                                                  
  return(list(
    lambda = opt$argmax, 
    cvl = opt$max, 
    predictions = best$predictions, 
    fold = groups, 
    fullfit = .makepenfit(best$fit, pu, prep$model, lambda1, opt$argmax*base2, fusedl =prep$fusedl,
      prep$orthogonalizer, prep$weights, prep$formula)
  ))
}
      
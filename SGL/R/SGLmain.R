SGL <- function(data, index, type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.1, nlam = 20, gamma = 0.8, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL){

  X.transform <- NULL

  if(standardize == TRUE){
    X <- data$x
    means <- apply(X,2,mean)
    X <- t(t(X) - means)
    var <- apply(X,2,function(x)(sqrt(sum(x^2))))
    X <- t(t(X) / var)
    data$x <- X
    X.transform <- list(X.scale = var, X.means = means)
  }

  if(type == "linear"){
    if(standardize == TRUE){
      intercept <- mean(data$y)
      data$y <- data$y - intercept
    }   
    Sol <- oneDim(data, index, thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, reset = reset, alpha = alpha)
    if(standardize == TRUE){
      Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, intercept = intercept, X.transform = X.transform)
    }
     if(standardize == FALSE){
      Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, X.transform = X.transform)
    }
  }

  if(type == "logit"){
    Sol <- oneDimLogit(data, index, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset)

    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, intercept = Sol$intercept, X.transform = X.transform)

  }

  if(type == "cox"){
    Sol <- oneDimCox(data, index, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset) 

  Sol = list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, X.transform = X.transform)
  }
  class(Sol) = "SGL"
    return(Sol)
}


cvSGL <- function(data, index = rep(1,ncol(data$x)), type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.05, nlam = 20, gamma = 0.8, nfold = 10, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL){

  if(standardize == TRUE){
    X <- data$x
    means <- apply(X,2,mean)
    X <- t(t(X) - means)
    var <- apply(X,2,function(x)(sqrt(sum(x^2))))
    X <- t(t(X) / var)
    data$x <- X
  }

  if(type == "linear"){
   if(standardize == TRUE){
      intercept <- mean(data$y)
      data$y <- data$y - intercept
    }   
    Sol <- linCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, reset = reset, alpha = alpha)

    if(standardize == TRUE){
       Sol$fit = list(beta = Sol$fit$beta, lambdas = Sol$fit$lambdas, intercept = intercept, step = step)
    }
  }
  if(type == "logit"){
    Sol <- logitCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset)

  }
  if(type == "cox"){
    Sol <- coxCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset)

  }

  Sol = list(fit = Sol$fit, lldiff = Sol$lldiff, lambdas = Sol$lambdas, type = type, llSD = Sol$llSD)

class(Sol) = "cv.SGL"

    return(Sol)
}

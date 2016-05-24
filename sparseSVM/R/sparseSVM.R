sparseSVM <- function (X, y, alpha = 1, gamma = 0.1, nlambda=100, lambda.min = ifelse(nrow(X)>ncol(X), 0.01, 0.05), lambda, 
                       preprocess = c("standardize", "rescale", "none"),  screen = c("ASR", "SR", "none"), max.iter = 1000, 
                       eps = 1e-5, dfmax = ncol(X)+1, penalty.factor=rep(1, ncol(X)), message = FALSE) {
  
  # Error checking
  preprocess <- match.arg(preprocess)
  screen <- match.arg(screen)
  if (alpha < 0 || alpha > 1) stop("alpha should be between 0 and 1")
  if (gamma < 0 || gamma > 1) stop("gamma should be between 0 and 1")
  if (missing(lambda) && nlambda < 2) stop("nlambda should be at least 2")
  if (length(penalty.factor)!=ncol(X)) stop("the length of penalty.factor should equal the number of columns of X")

  if (is.factor(y)) {
    levels <- levels(y)
  } else {
    levels <- unique(y)    
  }
  if (length(levels) != 2) stop("currently the function only supports binary classification")
  
  call <- match.call()
  # convert response to +1/-1 coding
  n <- length(y)
  yy <- double(n)
  yy[y == levels[1]] <- 1
  yy[y == levels[2]] <- -1
  # Include a column for intercept
  XX <- cbind(rep(1,n), X)
  penalty.factor <- c(0, penalty.factor) # no penalty for intercept term
  p <- ncol(XX)
  
  if(missing(lambda)) {
    lambda <- double(nlambda)
    user <- 0
  } else {
    nlambda <- length(lambda)
    user <- 1
  }
  
  # Flag for preprocessing and screening
  ppflag = switch(preprocess, standardize = 1, rescale = 2, none = 0)
  scrflag = switch(screen, ASR = 1, SR = 2, none = 0)
  # Fitting
  fit <- .C("sparse_svm", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1),
            as.double(XX), as.double(yy), as.double(penalty.factor), as.double(gamma), as.double(alpha), 
            as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), 
            as.integer(ppflag), as.integer(scrflag), as.integer(dfmax), as.integer(max.iter), 
            as.integer(user), as.integer(message))
  weights <- matrix(fit[[1]],nrow = p)
  iter <- fit[[2]]
  lambda <- fit[[3]]
  saturated <- fit[[4]]
  # Eliminate saturated lambda values
  ind <- !is.na(iter)
  weights <- weights[, ind]
  iter <- iter[ind]
  lambda <- lambda[ind]
  
  # Names
  vnames <- colnames(X)
  if (is.null(vnames)) vnames=paste0("V",seq(p-1))
  vnames <- c("(Intercept)", vnames)
  dimnames(weights) <- list(vnames, round(lambda, 4))

  # Output
  structure(list(call = call,
                 weights = weights,
                 iter = iter,
                 saturated = saturated,
                 lambda = lambda,
                 alpha = alpha,
                 gamma = gamma,
                 penalty.factor = penalty.factor[-1],
                 levels = levels),
            class = "sparseSVM")
}

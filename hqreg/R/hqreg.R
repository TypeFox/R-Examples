hqreg <- function (X, y, method = c("huber", "quantile", "ls"), gamma, tau = 0.5, alpha=1, nlambda=100, lambda.min = ifelse(nrow(X)>ncol(X), 0.001, 0.05), lambda, 
		  preprocess = c("standardize", "rescale", "none"),  screen = c("ASR", "SR", "none"), max.iter = 10000, eps = 1e-7, dfmax = ncol(X)+1, penalty.factor=rep(1, ncol(X)), message = FALSE) {
  
  # Error checking
  method <- match.arg(method)
  preprocess <- match.arg(preprocess)
  screen <- match.arg(screen)
  if (missing(lambda) && nlambda < 2) stop("nlambda should be at least 2")
  if (alpha < 0 || alpha > 1) stop("alpha should be between 0 and 1")
  if (method == "huber" && !missing(gamma) && gamma <= 0) stop("gamma should be positive for Huber loss")
  if (method == "quantile" && (tau < 0 || tau > 1)) stop("tau should be between 0 and 1 for quantile loss")
  if (length(penalty.factor)!=ncol(X)) stop("the length of penalty.factor should equal the number of columns of X")

  call <- match.call()
  # Include a column for intercept
  n <- nrow(X)
  XX <- cbind(rep(1,n), X)
  penalty.factor <- c(0, penalty.factor) # no penalty for intercept term
  p <- ncol(XX)

  shift <- 0
  if (method == "huber") {
    shift <- mean(y)
    yy <- y-shift
    if(missing(gamma)) gamma <- quantile(abs(yy), 0.05)
  } else if (method == "ls") {
    shift <- mean(y)
    yy <- y-shift
  } else if (method == "quantile") {
    shift <- quantile(y, tau)
    yy <- y-shift
    # initialize gamma to determine lambda sequence
    gamma <- quantile(abs(yy), 0.05)
    if (gamma < 0.00001) gamma = 0.00001
  }

  # Constant for quantile loss
  c <- 2*tau-1

  # For Huber loss, scale gamma accordingly when the response is standardized


  # Setup vector d for generating the lambda sequence
  d <- 0
  user <- 0
  if (missing(lambda)) {
    lambda <- double(nlambda)
    d <- initDeriv(yy, X, method, gamma, c, penalty.factor)
  } else {
    nlambda <- length(lambda)
    user <- 1
  }
  
  # Flag for preprocessing and screening
  ppflag = switch(preprocess, standardize = 1, rescale = 2, none = 0)
  scrflag = switch(screen, ASR = 1, SR = 2, none = 0)
  # Fitting
  if (alpha > 0) {
    if (method == "huber") {
      fit <- .C("huber", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy), as.double(d), as.double(penalty.factor), 
                as.double(gamma), as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), as.integer(ppflag),
                as.integer(scrflag), as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(message))
    } else if (method == "quantile") {
      fit <- .C("quant", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy), as.double(d), as.double(penalty.factor), 
                as.double(gamma), as.double(tau), as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), 
                as.integer(ppflag), as.integer(scrflag), as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(message))
    } else {
      fit <- .C("squared", double(p*nlambda), integer(nlambda), as.double(lambda), integer(1), integer(1), as.double(XX), as.double(yy), as.double(d), as.double(penalty.factor), 
                as.double(alpha), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), as.integer(ppflag), as.integer(scrflag),
                as.integer(dfmax), as.integer(max.iter), as.integer(user), as.integer(message))
    }
    beta <- matrix(fit[[1]],nrow = p)
    iter <- fit[[2]]
    lambda <- fit[[3]]
    saturated <- fit[[4]]
    nv <- fit[[5]]
    # Eliminate saturated lambda values
    ind <- !is.na(iter)
    beta <- beta[, ind]
    iter <- iter[ind]
    lambda <- lambda[ind]
  } else {
    if (method == "huber") {
      fit <- .C("huber_l2", double(p*nlambda), integer(nlambda), as.double(lambda), as.double(XX), as.double(yy), as.double(d), as.double(penalty.factor), 
                as.double(gamma), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), as.integer(ppflag),
                as.integer(max.iter), as.integer(user), as.integer(message))
    } else if (method == "quantile") {
      fit <- .C("quantile_l2", double(p*nlambda), integer(nlambda), as.double(lambda), as.double(XX), as.double(yy), as.double(d), as.double(penalty.factor), 
                as.double(gamma), as.double(tau), as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), as.integer(ppflag),
                as.integer(max.iter), as.integer(user), as.integer(message))      
    } else {
      fit <- .C("squared_l2", double(p*nlambda), integer(nlambda), as.double(lambda), as.double(XX), as.double(yy), as.double(d), as.double(penalty.factor), 
                as.double(eps), as.double(lambda.min), as.integer(nlambda), as.integer(n), as.integer(p), as.integer(ppflag),
                as.integer(max.iter), as.integer(user), as.integer(message))      
    }
    beta <- matrix(fit[[1]],nrow = p)
    iter <- fit[[2]]
    lambda <- fit[[3]]
    saturated <- 0
    nv <- 0
  }
 
  # Intercept
  beta[1,] <- beta[1,] + shift
  
  # Names
  vnames <- colnames(X)
  if (is.null(vnames)) vnames=paste0("V",seq(p-1))
  vnames <- c("(Intercept)", vnames)
  dimnames(beta) <- list(vnames, round(lambda, 4))

  # Output
  structure(list(call = call,
                 beta = beta,
                 iter = iter,
                 saturated = saturated,
                 lambda = lambda,
                 alpha = alpha,
                 gamma = if (method == "huber") gamma else NULL,
                 tau = if (method == "quantile") tau else NULL,
                 penalty.factor = penalty.factor[-1],
                 method = method,
                 nv = nv),
            class = "hqreg")
}

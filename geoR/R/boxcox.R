##
## Box-Cox transformation in the package geoR
## -------------------------------------------------------------
##

"boxcox.geodata" <- function(object, trend = "cte", ...)
{
  xmat <- unclass(trend.spatial(trend = trend, geodata = object))
  if (nrow(xmat) != length(object$data)) 
    stop("coords and trend have incompatible sizes")
  MASS::boxcox(object$data ~ xmat + 0, ...)
}

"boxcoxfit" <-
  function(object, xmat, lambda, lambda2 = NULL, add.to.data = 0,...)
{
  call.fc <- match.call()
  data <- object + add.to.data
  if(is.null(lambda2) && any(data <= 0))
    stop("Transformation requires positive data")
  ##
  data <- as.vector(data)
  n <- length(data)
  if(missing(xmat)) xmat <- rep(1, n)
  xmat <- as.matrix(xmat)
  if(any(xmat[,1] != 1)) xmat <- cbind(1, xmat)
  ## do not reverse order of the next two lines:
  xmat <- xmat[!is.na(data),,drop=FALSE]
  data <- data[!is.na(data)]
  n <- length(data)
  ##
  beta.size <- ncol(xmat)
  if(nrow(xmat) != length(data))
    stop("xmat and data have incompatible lengths")
  ##  lik.method <- match.arg(lik.method, choices = c("ML", "RML"))
  lik.method <- "ML"
  ##
  if(all(data > 0)) absmin <- 0
  else absmin <- abs(min(data)) + 0.00001 * diff(range(data))
  if(!is.null(lambda2)){
    if(missing(lambda)) lambda.ini <- seq(-2, 2, by=0.2)
    else lambda.ini <- lambda
    lambda2.ini <- 0
    if(isTRUE(lambda2)) lambda2.ini <- absmin
    else if(mode(lambda2) == "numeric") lambda2.ini <- lambda2
    lambdas.ini <- as.matrix(expand.grid(lambda.ini, lambda2.ini))
    ##
    if(length(as.matrix(lambdas.ini)) > 2){
      lamlik <- apply(lambdas.ini, 1, .negloglik.boxcox, data=data + absmin,
                      xmat=xmat, lik.method=lik.method)
      lambdas.ini <- lambdas.ini[which(lamlik == min(lamlik)),]
    }
    lambdas.ini <- unname(drop(lambdas.ini))
    lik.lambda <- optim(par=lambdas.ini, fn = .negloglik.boxcox,
                        method="L-BFGS-B",
                        #hessian = TRUE, 
                        lower = c(-Inf, absmin), 
                        data = data, xmat = xmat, lik.method = lik.method)
  }
  else{
    lik.lambda <- optimize(.negloglik.boxcox, interval = c(-5, 5), data = data,
                           xmat = xmat, lik.method = lik.method)
    lik.lambda <- list(par = lik.lambda$minimum, value = lik.lambda$objective,
                       convergence = 0, message = "function optimize used")
  }
  ##
  ##  hess <- sqrt(diag(solve(as.matrix(lik.lambda$hessian))))
  lambda.fit <- lik.lambda$par
  if(length(lambda.fit) == 1) lambda.fit <- c(lambda.fit, 0)
  data <- data + lambda.fit[2]
  ##
#  if(abs(lambda.fit[1]) < 0.0001) yt <- log(data)
  if(isTRUE(all.equal(unname(lambda.fit[1]),0))) yt <- log(data)
  else yt <- ((data^lambda.fit[1]) - 1)/lambda.fit[1]
  beta <- solve(crossprod(xmat), crossprod(xmat, yt))
  mu <- drop(xmat %*% beta)
  sigmasq <- sum((yt - mu)^2)/n
  if(lik.method == "ML")
    loglik <- drop((-(n/2) * (log(2*pi) + log(sigmasq) + 1)) + (lambda.fit[1]-1) * sum(log(data)))
  ##  if(lik.method == "RML")
  ##    loglik <- drop(-lik.lambda$value - (n/2)*log(2*pi) - (n-beta.size)*(log(n) - 1))
  ##
  temp <- 1 + lambda.fit[1] * mu
  fitted.y <- ((temp^((1/lambda.fit[1]) - 2)) *
               (temp^2 + ((1-lambda.fit[1])/2) * sigmasq))
  variance.y <-  (temp^((2/lambda.fit[1]) - 2)) * sigmasq
  if(beta.size == 1){
    fitted.y <- unique(fitted.y)
    variance.y <- unique(fitted.y)
  }
  ##
  beta <- drop(beta)
  if(length(beta) > 1)
    names(beta) <- paste("beta", 0:(beta.size-1), sep="")
  if(length(lik.lambda$par) == 1) lambda.fit <- lambda.fit[1]
  if(length(lik.lambda$par) == 2) names(lambda.fit) <- c("lambda", "lambda2")
  res <- list(lambda = lambda.fit, beta.normal = drop(beta),
              sigmasq.normal = sigmasq, 
              loglik = loglik, optim.results = lik.lambda)
  ## res$hessian <- c(lambda = hess) 
  res$call <- call.fc
  oldClass(res) <- "boxcoxfit"
  return(res)
}

".negloglik.boxcox" <-
  function(lambda.val, data, xmat, lik.method = "ML")
{
  if(length(lambda.val) == 2){
    data <- data + lambda.val[2]
    lambda <- lambda.val[1]
  }
  else lambda <- lambda.val
  lambda <- unname(lambda)
  n <- length(data)
  beta.size <- ncol(xmat)
  if(isTRUE(all.equal(unname(lambda), 0))) yt <- log(data)
  else yt <- ((data^lambda) - 1)/lambda
  beta <- solve(crossprod(xmat), crossprod(xmat, yt))
  ss <- sum((drop(yt) - drop(xmat %*% beta))^2)
  if(lik.method == "ML")
    neglik <- (n/2) * log(ss) - ((lambda - 1) * sum(log(data)))
  if(lik.method == "RML"){
    xx <- crossprod(xmat)
    if(length(as.vector(xx)) == 1)
      choldet <- 0.5 * log(xx)
    else
      choldet <- sum(log(diag(chol(xx))))
    neglik <- ((n-beta.size)/2) * log(ss) + choldet -
      ((lambda - 1) * sum(log(data)))
  }
  if(mode(neglik) != "numeric") neglik <- Inf
  return(drop(neglik))
}

"print.boxcoxfit" <-
  function(x, ...)
{
  if(length(x$lambda) == 1) names(x$lambda) <- "lambda"
  if(length(x$beta.normal) == 1) names(x$beta.normal) <- "beta"
  res <- c(x$lambda, x$beta.normal, sigmasq = x$sigmasq.normal)
  cat("Fitted parameters:\n")
  print(res)
  cat("\nConvergence code returned by optim: ")
  cat(x$optim.results$convergence)
  cat("\n")
  return(invisible())
}

"plot.boxcoxfit" <-
  function(x, hist = TRUE, data = eval(x$call$object), ...)
{
  if(is.null(data)) stop("data object not provided or not found")
  if(!is.null(x$call$xmat) && ncol(eval(x$call$xmat)) > 1)
    stop("plot.boxcoxfit not valid when covariates are included")
  if(!is.null(x$call$add.to.data))
    data <- data + eval(x$call$add.to.data)
  y <- data
  rd <- range(y)
  obj <- x
  x <- NULL
  f <- function(x, res = obj){
    dboxcox(x, lambda=res$lambda[1],
            lambda2 = ifelse(res$lambda[2], res$lambda[2], 0),
            mean=res$beta, sd=sqrt(res$sigmasq))
  } 
  ldots <- list()
  if(hist){
    if(is.null(ldots$ylim)){
      lim.hist <- max(hist(y, prob=TRUE, ...)$dens)
      lim.dens <- max(f(seq(rd[1], rd[2], l=200)))
      hist(y, prob=TRUE, ylim=c(0, max(lim.hist, lim.dens)), ...)
    }
    else
      hist(y, prob=TRUE, ...)
    curve(f(x), from=rd[1], to=rd[2], add=TRUE)
  }
  else{
    if(is.null(ldots$from)) ini <- rd[1]
    else ini <- ldots$from
    if(is.null(ldots$to)) fim <- rd[2]
    else fim <- ldots$to
    curve(f(x), from=ini, to=fim, ...)
  }
  return(invisible())
}

"lines.boxcoxfit" <-
  function(x, data = eval(x$call$object), ...)
{
  if(is.null(data)) stop("data object not provided or not found")
  if(!is.null(x$call$xmat) && ncol(eval(x$call$xmat)) > 1)
    stop("lines.boxcoxfit not valid when covariates are included")
  y <- data
  rd <- range(y)
  obj <- x
  x <- NULL
  rd <- range(data)
  ldots <- list()
  if(is.null(ldots$from)) ini <- rd[1]
  else ini <- ldots$from
  if(is.null(ldots$to)) fim <- rd[2]
  else fim <- ldots$to
  f <- function(x, res = obj){
    dboxcox(x, lambda=res$lambda,
            lambda2 = ifelse(res$lambda[2], res$lambda[2], 0),
            mean=res$beta, sd=sqrt(res$sigmasq))
  }
  curve(f(x), from=rd[1], to=rd[2], add=TRUE, ...)
  return(invisible())
}

"rboxcox" <-
  function(n, lambda, lambda2 = NULL, mean = 0, sd = 1)
{
  if(is.null(lambda2)) lambda2 <- 0
  if(is.na(lambda2)) lambda2 <- 0
  xn <- rnorm(n = n, mean = mean, sd = sd)
  if(isTRUE(all.equal(unname(lambda), 0))) xbc <- exp(xn)
  else{
    xbc <- rep(NA, n)
    ind <- xn < -1/lambda
    sum.ind <- sum(ind)
    if(sum.ind > 0)
      cat(paste("rboxcox: WARNING ", sum.ind, "values truncated to 0\n")) 
    xn[ind] <- -1/lambda
    xbc <- ((xn * lambda) + 1)^(1/lambda)
  }
  return(xbc - lambda2)
}

"dboxcox" <-
  function(x, lambda, lambda2 = NULL, mean = 0, sd = 1)
{
  if(is.null(lambda2)) lambda2 <- 0
  if(is.na(lambda2)) lambda2 <- 0
  x <- x + lambda2
  lx <- length(x)
  dval <- rep(0, lx)
  for(i in 1:lx){
    if(x[i] <=0) dval[i] <- 0
    else{
      if(isTRUE(all.equal(unname(lambda), 0))) xt <- log(x[i])
      else xt <- ((x[i]^lambda) - 1)/lambda 
      dval[i] <- ((1/sqrt(2*pi)) * (1/sd) * x[i]^(lambda-1) *
                  exp(-((xt-mean)^2)/(2*sd^2)))
    }
  }
  return(dval)
}

"BCtransform" <-
  function(x, lambda, add.to.data = 0,
           inverse = FALSE, log.jacobian = FALSE)
{
  x <- x + add.to.data
  if(inverse){
    if(log.jacobian)
      stop("options log.jacobian not allowed with inverse = TRUE")
    if(isTRUE(all.equal(unname(lambda), 0)))
      x <- exp(x)
    else{
      if(lambda > 0)
        x[x < (-1/lambda)] <- -1/lambda
      if(lambda < 0)
        x[x > (-1/lambda)] <- -1/lambda
      x <- ((x * lambda) + 1)^(1/lambda)
    }
    return(list(data = x))
  }
  else{
    if(!isTRUE(all.equal(unname(lambda), 1))){
      if(any(x <= 0))
        stop("Transformation requires positive data")
      if(log.jacobian){
        Jdata <- x^(lambda - 1)
        if(any(Jdata <= 0))
          temp.list$log.jacobian <- log(prod(Jdata))
        else temp.list$log.jacobian <- sum(log(Jdata))
        Jdata <- NULL
      }
      if(isTRUE(all.equal(unname(lambda), 0)))
        x <- log(x)
      else x <- ((x^lambda) - 1)/lambda
      if(any(c(is.na(x), is.nan(x))))
        stop("transformation has generated NA or NaN values")
      if(any(!is.finite(x)))
        stop("transformation has generated Inf values")
    }
    else
      if(log.jacobian) log.jacobian <- 0
    if(log.jacobian)
      return(list(data = x, log.jacobian = log.jacobian))
    else
      return(list(data = x))
  }
}

"backtransform.moments" <-
  function(lambda, mean, variance, distribution,
           simul.back = FALSE, n.simul = 1000) 
{
  ##
  ## This function is called by krige.bayes and krige.conv functions 
  ## in the package geoR to backtransform predictions when the original variable
  ## was transformed by using the Box-Cox transformation.
  ##
  ## WARNING: The Box-Cox transformation internal to the functions are:
  ##    Z = ((Y^lambda) -1)/lambda  for Y != 0 and Y != 1
  ##    Z = log(Y) for Y = 0
  ##    NO TRANSFORMATION (i.e. Z = Y) is performed when lambda = 1
  ##
  ## The transformations can be done:
  ##    - by analytical approximation (setting simul.back = FALSE)
  ##  or  
  ##    - by simulation (setting simul.back = TRUE)
  res <- list(mean = numeric(), variance = numeric(), distribution = character()) 
  ni <- length(mean)
  mean <- as.vector(mean)
  variance <- as.vector(variance)
  if (ni != length(variance)) stop("mean and variances must have same length")
#  if(abs(lambda-1) > 0.001){
  if(!isTRUE(all.equal(unname(lambda), 1))){
    if(simul.back){
      res$distribution <- "back-transformed (Box-Cox) from Gaussian by simulation"
#      ap.warn <- options()$warn
#      options(warn = -1)
      temp.data <- suppressWarnings(matrix(rnorm(n = ni * n.simul,
                                mean = mean,
                                sd = sqrt(variance)),
                          nrow = ni))
#      options(warn = ap.warn)
      ind.zero <- (variance < 1e-12)
      temp.data[ind.zero,  ] <- mean[ind.zero]
      remove(ind.zero)
      temp.data <- BCtransform(x = temp.data, lambda = lambda, 
                               inverse = TRUE)$data
      if(lambda < 0) {
        res$mean  <-  "resulting distribution has no mean for negative lambda. Medians returned"
        res$variance  <-  "resulting distribution has no variance for negative lambda"
      }
      else{
        res$mean <- as.vector(rowMeans(temp.data))
        res$variance <- as.vector(apply(temp.data, 1, var))
        quants <- apply(temp.data, 1, quantile, prob=c(.025,.5, .975))
        res$median <- drop(quants[2,])
        res$uncertainty <- drop((quants[3,] - quants[1,])/4)
      }
    }
    else{
      res$distribution <- "back-transformed (Box-Cox) from Gaussian"
      if(isTRUE(all.equal(unname(lambda), 0))) {
        temp <- mean
        res$mean <- exp(mean + 0.5 * (variance))
        res$variance <- (exp(2 * temp + variance)) * expm1(variance)
        res$median <- exp(mean)
#        res$mode <- exp(mean - variance)
      }
      else{
        temp <- 1 + (lambda * mean)
        res$mean <- (temp^((1/lambda)-2) *
                     ((temp^2) + ((1-lambda)/2) * variance)) 
        if(isTRUE(all.equal(unname(lambda), 0.5)))
          res$variance <- variance * ((variance/8) + temp^2)
        else
          res$variance <- temp^((2/lambda)-2) *  variance
        res$median <- temp^(1/lambda)
#        res$mode <- ((temp+sqrt(1+4*lambda*variance*(lambda-1)+(temp-1)*(1+temp)))/2)^(1/lambda)
      }
      res$uncertainty <- (qnorm(.975, mean = mean, sd = sqrt(variance)) -
                          qnorm(.025, mean = mean, sd = sqrt(variance)))/4
    }
  }
  else{
    res$distribution <- ifelse(missing(distribution), NULL, distribution)
    res$mean <- mean
    res$variance <- variance
  }
  return(res) 
}


profile.maxstab <- function(fitted, param, range, n = 10,
                            plot = TRUE, conf = 0.90,
                            method = "RJ", square = "chol", ...){

  if (!(method %in% c("RJ", "CB", "none")))
    stop("'method' must be one of 'RJ', 'CB' or 'none'")

  if (!(square %in% c("chol", "svd")))
    stop("'square' must be one of 'chol' or 'svd'")

  if (is.null(fitted$std.err)){
    warning("Profile confidence intervals cannot be computed whithout standard errors.")
    method <- "none"
  }

  param.names <- names(fitted$fitted.values)
  n.param <- length(param.names)
  idx.prof <- which(param.names == param)
  fixed.param <- fitted$fixed
  fixed.names <- names(fixed.param)

  ihessian <- fitted$ihessian
  var.cov <- fitted$var.cov

  if (method == "CB"){
    ##Chandler and Bate approach:
    ivar.cov <- solve(var.cov)
    hessian <- solve(ihessian)

    if (square == "svd"){
      svd.hessian <- svd(hessian)
      svd.ivar.cov <- svd(ivar.cov)
      M <- svd.hessian$u %*% diag(sqrt(svd.hessian$d)) %*%
        t(svd.hessian$u)
      Madj <- svd.ivar.cov$u %*% diag(sqrt(svd.ivar.cov$d)) %*%
        t(svd.ivar.cov$u)
    }

    else{
      M <- chol(hessian)
      Madj <- chol(ivar.cov)
    }

    C <- solve(M) %*% Madj
    colnames(C) <- rownames(C) <- colnames(ihessian)
  }

  if (method == "RJ"){
    ##Rotnitzky and Jewell approach:
    hessian <- solve(ihessian[param, param])
    var.cov <- var.cov[param, param]

    Q <- var.cov %*% hessian
    eigen.val <- eigen(Q)$values
  }

  start <- as.list(fitted$fitted.values)
  start <- start[-idx.prof]

  fixed.values <- seq(range[1], range[2], length = n)
  optfun <- fitted$lik.fun

  if (method == "CB"){
    par <- NULL
    for (fixed.val in range){
      if (length(fixed.names) > 0)
        body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-1),"]", collapse = ","),
                                ",", paste(fixed.names, "=", as.numeric(fixed.param),
                                           collapse = ","), ",", param, "=", fixed.val, ")"))

      else
        body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-1),"]", collapse = ","),
                                ",", param, "=", fixed.val, ")"))

      par <- rbind(par, optim(start, optfun)$par)
    }

    par <- cbind(range, par)
    colnames(par) <- c(param, param.names[-idx.prof])
    par <- par[,names(fitted$fitted)]
    par.adj <- t(fitted(fitted) + C %*% (t(par) - fitted(fitted)))
    range.adj <- par.adj[,param]
    fixed.values <- seq(range.adj[1], range.adj[2], length = n)

  }

  llik <- rep(NA, n)
  par <- matrix(NA, ncol = n.param - 1, nrow = n)

  for (i in 1:n){
    fixed.val <- fixed.values[i]

    if (length(fixed.names) > 0)
      body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-1),"]", collapse = ","),
                              ",", paste(fixed.names, "=", as.numeric(fixed.param),
                                         collapse = ","), ",", param, "=", fixed.val, ")"))

    else
      body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-1),"]", collapse = ","),
                              ",", param, "=", fixed.val, ")"))

    if (optfun(unlist(start)) >= 1e15)
      reltol <- 1e-10

    else
      reltol <- 1e-6

    opt <- optim(start, optfun, control = list(reltol = reltol, maxit = 10000))
    llik[i] <- -opt$value
    par[i,] <- opt$par
  }

  llik[llik <= -1e15] <- NA
  ans <- cbind(fixed.values, llik, par)
  colnames(ans) <- c(param, "llik", param.names[-idx.prof])
  ans <- ans[,c("llik", names(fitted$fitted))]

  if (plot){
    llikmax <- fitted$logLik

    if (method == "none")
      plot(ans[,param], ans[,"llik"], xlab = param, ylab = "log-likelihood", ...)

    if (method == "RJ"){
      b.conf <- llikmax - .5 * eigen.val * qchisq(conf, 1)
      plot(ans[,param], ans[,"llik"], xlab = param, ylab = "log-likelihood", ...)
      abline(h = llikmax)
      abline(h = b.conf)
    }

    if (method == "CB"){
      ans[,-1] <- t(fitted(fitted) + solve(C) %*% (t(ans[,-1]) - fitted(fitted)))
      b.conf <- llikmax - .5 * qchisq(conf, 1)
      plot(ans[,param], ans[,"llik"], xlab = param, ylab = "log-likelihood", ...)
      abline(h = llikmax)
      abline(h = b.conf)
    }
  }

  return(ans)

}

profile2d.maxstab <- function(fitted, params, ranges, n = 10,
                              plot = TRUE, ...){
  param.names <- names(fitted$fitted.values)
  n.param <- length(param.names)
  idx.prof <- which(param.names %in% params)
  fixed.param <- fitted$fixed
  fixed.names <- names(fixed.param)

  start <- as.list(fitted$fitted.values)
  start <- start[-idx.prof]

  fixed.values1 <- seq(ranges[1,1], ranges[1,2], length = n)
  fixed.values2 <- seq(ranges[2,1], ranges[2,2], length = n)

  llik <- matrix(NA, n, n)

  optfun <- fitted$lik.fun

  for (i in 1:n){
    for (j in 1:n){
      fixed.val <- c(fixed.values1[i], fixed.values2[j])
      ##We need to modify the body of optfun for each fixed value
      if (length(fixed.names) > 0)
        body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-2),"]", collapse = ","),
                                ",", paste(fixed.names, "=", as.numeric(fixed.param),
                                           collapse = ","), ",",
                                paste(params, "=", fixed.val, collapse = ","), ")"))

      else
        body(optfun) <- parse(text = paste("nplk(", paste("p[",1:(n.param-2),"]", collapse = ","),
                                ",", paste(params, "=", fixed.val, collapse = ","), ")"))


      llik[i,j] <- -optim(start, optfun)$value
    }
  }

  ans <- list(coord = cbind(fixed.values1, fixed.values2), llik = llik)
  colnames(ans$coord) <- c(params)

  if (plot)
    filled.contour(fixed.values1, fixed.values2, llik,
                   xlab = params[1], ylab = params[2], ...)

  return(ans)

}

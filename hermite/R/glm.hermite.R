glm.hermite <- function (formula, data, link = "log", start = NULL, m = NULL) 
{
  Call <- match.call()
  m.fixed   <- FALSE
  intercept <- FALSE
  if (!is.null(m)) 
  {
    m.fixed <- TRUE  
    if (floor(m) != m ) {
      stop("improper m parameter specification")
    }
  }
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call()
  mm <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, mm)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf)
  Y <- model.response(mf, "numeric")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  Xnames <- dimnames(X)[[2L]]
  kx <- NCOL(X)
  if (all(X[, 1] == 1)) 
    intercept <- TRUE
  maxY <- max(Y)
  n <- length(Y)
  linkstr <- link
  linkobj <- make.link(linkstr)
  linkinv <- linkobj$linkinv
  prob <- function(y, mu, d) {
    ny <- length(mu)
    maxy <- max(m, max(y), 2)
    p <- matrix(NA, ny, maxy + 1)
    p[, 1] <- exp(mu * (-1 + (d - 1)/m))
    p[, 2:m] <- vapply(1:(m - 1), function(k) p[, 1] * ((mu^k)/(factorial(k))) * 
                           (((m - d)/(m - 1))^k), rep(1, ny))
    if(m <= maxy)
    {
      for (i in m:maxy) {
           p[, i + 1] <- mu * (p[, i - m + 1] * (d - 1) + p[, i] * (m - d))/(i * (m - 1))
      }
    }  
    IJ <- as.matrix(cbind(row = 1:ny, col = y + 1))
    p[IJ] <- ifelse(p[IJ] < 0, 0, p[IJ])
    return(p[IJ])
  }
  
  loglik <- function(parms) {
    mu <- as.vector(linkinv(X %*% parms[1:kx]))
    d <- parms[kx + 1]
    if (d < 1) d <- 1.1
    loglikh <- sum(log(prob(Y, mu, d))[is.finite(log(prob(Y, mu, d)))])
    loglikh
  }
  mloglik <- function(parms) {
    -1*loglik(parms)
  }
    
  m.naive <- FALSE
  if (is.null(m)) {
    d <- var(Y)/mean(Y)
    p0 <- length(Y[Y == 0])/length(Y)
    m.n <- ifelse(p0 != 0 & dim(X)[2]==1, round((d - 1)/(1 + log(p0)/mean(Y))), 
                NA)
    if (!is.na(m.n) & m.n < 0) m.n <- 1
    m <- m.n
    if (!is.na(m) & !m.fixed) m.naive <- TRUE
  }
  
  if (is.na(m)) m <- 1
  
  if (is.null(start) & link == "identity") start <- c(rep(1, NCOL(X)), 1.1)
  if (is.null(start) & link == "log") 
    start <- c(as.numeric(glm.fit(X, Y, family = poisson(link = link))$coefficients), 
               1.1)

  ex.solMLE <- function(m)
    {
      if (dim(X)[2]>1) return(TRUE)
      res <- TRUE
      fm <- 0
      for (i in 1:length(Y)) {
      mult <- 1
      for (j in 0:(m - 1)) {
        mult <- mult * (Y[i] - j)
      }
      fm <- fm + mult
      }
      fm <- fm/length(Y)
      if (fm <= (mean(Y))^m) {
        res <- FALSE
      }
      return(res)
    }
    if (m.fixed & m > 1) 
    {
      if (ex.solMLE(m)==FALSE)
      {
        coefs <- c(mean(Y), 1)
        llik <- loglik(c(coefs[1], coefs[2]))
        hess <- length(Y)/mean(Y)
        w <- 2 * (loglik(c(coefs[1], coefs[2])) - loglik(c(coefs[1], 1)))
        pval <- pchisq(w, 1, lower.tail = F)
        warning("MLE equations have no solution")
        output <- list()
        output$coefs <- coefs
        output$loglik <- llik
        output$hess <- hess
        output$w <- w
        output$pval <- pval
        class(output) <- "glm.hermite"
        attr(output, "Call") <- Call 
        attr(output, "x") <- X
        return(output)
      }
    }

  if (link == "log") {
    A <- matrix(c(rep(0, kx), 1, c(rep(0,kx), -1)), 2, kx + 1, byrow=TRUE)
    B <- c(1,m)
  }
  else {
    if (intercept) {
      A <- rbind(cbind(X,0), c(rep(0,kx), 1), c(1,rep(0,kx)), c(rep(0,kx), -1)) 
      B <- c(rep(0,n),1,0,m)
    }
    
    if (!intercept)
    {
      if (m.fixed)
      {
        if (m > 1) fit2 <- optim(start, mloglik, method="L-BFGS-B", lower=c(rep(-Inf, NCOL(X)), 1), upper=c(rep(Inf,NCOL(X)),m), hessian=TRUE)
        if (m == 1)
        {
          if(class(try(glm(formula, family = poisson(link = link), data=data, start=start[-length(start)]),silent = TRUE))!="try-error")
          {
            fit2 <- suppressWarnings(glm(formula, family = poisson(link = link), data=data, start=start[-length(start)]))
            fit2$maximum <- logLik(fit2)
          }else{
            stop("Error on glm Poisson fit")
          }
        }
      }
      if (!m.fixed)
      {
        if(class(try(glm(formula, family = poisson(link = link), data=data, start=start[-length(start)]),silent = TRUE))!="try-error")
        {
          fit2 <- suppressWarnings(glm(formula, family = poisson(link = link), data=data, start=start[-length(start)]))
          fit2$maximum <- logLik(fit2)
        }else{
          m <- 2
          if (ex.solMLE(m))
          {
            fit2 <- optim(start, mloglik, method="L-BFGS-B", lower=c(rep(-Inf, NCOL(X)), 1), upper=c(rep(Inf,NCOL(X)),m), hessian=TRUE)
            fit2$maximum <- -fit2$value
          }else{
            warning("MLE equations have no solution for m=", m)
            fit2$maximum <- -Inf
          }
        }
        m.f  <- m
        j    <- 3
        while((j <= m.f+1 & j <= min(max(Y), 10)) | (m.f == 1 & j <= 5))
        {
          m <- j
          if (ex.solMLE(m))
          {
            fit <- optim(start, mloglik, method="L-BFGS-B", lower=c(rep(-Inf, NCOL(X)), 1), upper=c(rep(Inf,NCOL(X)),m), hessian=TRUE)
            fit$maximum <- -fit$value
            if (!is.null(fit$maximum) & fit$convergence==0)
            {
              if (fit$maximum > fit2$maximum)
              {
                fit2 <- fit
                m.f  <- m
              }
            }
          }else{
            warning("MLE equations have no solution for m=", m)
            fit2$maximum <- -Inf
          }
          m <- m.f
          j <- j + 1
        }
      }  
    
      coefs <- c(fit2$par, m)
      l1 <- length(fit2$par)
      l2 <- l1 - 1
      names(coefs) <- c(Xnames, "dispersion.index", "order")
      mu <- as.vector(linkinv(X %*% coefs[1:kx]))
      hess <- fit2$hessian
      if (is.null(hess)) hess <- hessian(fit2)
      if (m == 1) w <- ifelse(l2 > 1, 2 * (fit2$maximum - loglik(c(coefs[1], coefs[2:l2], 1))), 2 * (fit2$maximum - loglik(c(coefs[1], 1))))
      if (m > 1) w <- ifelse(l2 > 1, 2 * (-fit2$value - loglik(c(coefs[1], coefs[2:l2], 1))), 2 * (-fit2$value - loglik(c(coefs[1], 1))))
      pval <- pchisq(w, 1, lower.tail = F)/2
      output <- list()
      output$coefs <- coefs
      output$loglik <- -fit2$value
      output$vcov <- solve(hess)
      if(class(fit2)!="list") output$vcov <- vcov(fit2)
      output$hess <- hess
      output$fitted.values <- mu
      output$w <- w
      output$pval <- pval
      class(output) <- "glm.hermite"
      attr(output, "Call") <- Call 
      attr(output, "x") <- X
      return(output)
    }
  }    
  constraints <- list(ineqA = A, ineqB = B)
  if (any(A %*% start + B < 0)) stop("initial value not feasible")
  
  if (m.fixed & m>1) fit2 <- maxLik(logLik = loglik, start = start, constraints = constraints, iterlim = 1000)
  if (m.fixed & m==1) 
  {
    fit2 <- suppressWarnings(glm(formula, family = poisson(link = link), data=data, start=start[-length(start)]))
    fit2$maximum <- logLik(fit2)
  }
  if (!m.fixed)
  {
    fit2 <- suppressWarnings(glm(formula, family = poisson(link = link), data=data, start=start[-length(start)]))
    fit2$maximum <- logLik(fit2)
    m.f  <- m
    j    <- 2
    while((j <= m.f+1 & j <= min(max(Y), 10)) | (m.f == 1 & j <= 5))
    {
      m <- j
      if (ex.solMLE(m))
      {
        if (link=="log") B <- c(1, m)
        if (link=="identity") B <- B <- c(rep(0,n),1,0,m)
        constraints <- list(ineqA = A, ineqB = B)
        fit <- maxLik(logLik = loglik, start = start, constraints = constraints, iterlim = 1000)
        if (!is.null(fit$maximum))
        {
          if (fit$maximum > fit2$maximum)
          {
            fit2 <- fit
            m.f  <- m
          }
        }
      }else{
        warning("MLE equations have no solution for m=", m)
        fit$maximum <- -Inf
      }
      m <- m.f
      j <- j + 1
    }
  }

  if (m==1)
  {
    fit2 <- suppressWarnings(glm(formula, family = poisson(link = link), data=data, start=start[-length(start)]))
    coefs <- c(fit2$coefficients, 1, m)
    names(coefs) <- c(Xnames, "dispersion.index", "order")
    mu <- as.vector(linkinv(X %*% coefs[1:kx]))
    output <- list()
    output$coefs <- coefs
    output$loglik <- as.numeric(logLik(fit2))
    output$vcov <- vcov(fit2) 
    colnames(output$vcov) <- NULL
    rownames(output$vcov) <- NULL
    output$hess <- solve(output$vcov)
    output$fitted.values <- mu
    output$w <- NA
    output$pval <- NA
    class(output) <- "glm.hermite"
    attr(output, "Call") <- Call 
    attr(output, "x") <- X
    return(output)
  }
  
  if (m.naive==TRUE)
  {
    if(m.n > j) 
    { 
      if (link=="log") B <- c(1, m.n)
      if (link=="identity") B <- B <- c(rep(0,n),1,0,m.n)
      constraints <- list(ineqA = A, ineqB = B)
      fit3 <- maxLik(logLik = loglik, start = start, constraints = constraints, iterlim = 1000)
      if (!is.null(fit3$maximum))
      {
        if (fit3$maximum > fit2$maximum)
        {
          fit2 <- fit3
          m  <- m.n
        }
      }
    }
  }
  
  coefs <- c(fit2$estimate, m)
  l1 <- length(fit2$estimate)
  l2 <- l1 - 1
  names(coefs) <- c(Xnames, "dispersion.index", "order")
  mu <- as.vector(linkinv(X %*% coefs[1:kx]))
  hess <- hessian(fit2)
  w <- ifelse(l2 > 1, 2 * (fit2$maximum - loglik(c(coefs[1], 
                                                  coefs[2:l2], 1))), 2 * (fit2$maximum - loglik(c(coefs[1], 
                                                                                                 1))))
  pval <- pchisq(w, 1, lower.tail = F)/2
  output <- list()
  output$coefs <- coefs
  output$loglik <- fit2$maximum
  output$vcov <- solve(-hess)
  output$hess <- hess
  output$fitted.values <- mu
  output$w <- w
  output$pval <- pval
  class(output) <- "glm.hermite"
  attr(output, "Call") <- Call
  attr(output, "x") <- X
  return(output)
}
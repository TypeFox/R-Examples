##=================================================================
## Author: Y. Deville
## Rider's moments estimator for a mixture of two exponentials
## Might fail to exist and are then turned to NAs
##=================================================================

mom.mixexp2 <- function(x) {
  
  n <- length(x)

  if (any(x <= 0)) stop("x must contain only positive values")
  
  xbar <- mean(x)

  ## scale 
  ## xmod <- x / xbar
  
  Mp1 <- mean(x)
  Mp2 <- mean(x^2)
  Mp3 <- mean(x^3)

  a <- 6 * (2*Mp1^2- Mp2)
  b <- 2 * (Mp3 - 3*Mp1*Mp2)
  c <- 3*Mp2^2 - 2*Mp1*Mp3

  Delta <- b^2 - 4*a*c

  if (Delta < 0) {
    resest <- rep(NA, 3)
    names(resest) <- c("alpha", "lambda1", "lambda2")
    res <- list(estimate = resest,
                method = "moments")
    return(res)
  }

  theta1 <-  (-b - sqrt(Delta)) / 2 / a
  theta2 <- (-b + sqrt(Delta)) / 2 / a

  p <- (Mp1 - theta2) / (theta1 - theta2)

  ## check moments?
  if (FALSE) {
    test <- cbind(c(xbar, Mp2/2, Mp3/6),
                  c(p*theta1+(1-p)*theta2,
                    p*theta1^2+(1-p)*theta2^2,
                    p*theta1^3+(1-p)*theta2^3))
    colnames(test) <- c("Mp", "check")
    print(test)
    
  }

  if ( (theta1 < 0) || (theta2 < 0) || (p < 0) || (p > 1) ) {
    theta1 <- NA
    theta2 <- NA
    resest <- rep(NA, 3)
  } else if (theta1 > theta2) {
    resest <- c(p, 1/theta1, 1/theta2)
  } else {
    resest <- c(1-p, 1/theta2, 1/theta1)
  }
  
  names(resest) <- c("prob1", "rate1", "rate2")

  rest <- list(estimate = resest,
               method = "moments")

}

##==============================================================================
## Author: Yves Deville
## Find initial estimates for the three parameters of the mixexp2 dist.
## Firts try a moment estimation. If this fails, we use an empirical method
## based upon the cumulative hazard rate function H(y) which can be regressed
## on the order statistics.
##==============================================================================
  
ini.mixexp2 <- function(x, plot = FALSE) {

  res <- mom.mixexp2(x = x)

  needsmore <- any(is.na(res$estimate))

  if (needsmore) {
    
    pct.U <- c(90, 95, 85, 80, 75)
    pct.L <- c(10, 15, 20, 25, 30)
    
    if (any(x < 0)) stop("x must contain only values >= 0")
    
    n <- length(x)
    ## pct <- c(10, 90)
    xbar <- mean(x)
    xbarinv <- 1/xbar
    
    xs <- sort(x)
    H <- -log( (n+1-(1:n))/(n+1) )
    
    i1 <- 1
    i2 <- 1
    OK1 <- FALSE
    OK2 <- FALSE
    
    while( (i1 <= length(pct.U)) && (!OK1 || !OK2)) {
      
      ind.U <- (1:n) >=  round(n*pct.U[i1]/100)
      x1 <- xs[ind.U]
      y1 <- H[ind.U]
      coef1 <- lm(y1 ~ x1)$coef
      lambda1.ini <- coef1[2]
      OK1 <- (lambda1.ini < xbarinv)
      
      if (OK1) {
        
        pct.u <- pct.U[i1]
        den <- 1 - lambda1.ini*xbar
        i2 <- 1
        OK2 <- FALSE
        
        while( (i2 <= length(pct.L)) && !OK2) { 
          ind.L <- (1:n) <=  round(n*pct.L[i2]/100)
          x1 <- xs[ind.L]
          y1 <- H[ind.L] 
          lambdatilde <- lm(y1 ~ x1 -1)$coef
          
          if (lambdatilde > lambda1.ini) {
            lambda2.ini <- (lambdatilde - lambda1.ini)  / den
            OK2 <-(lambda2.ini > lambdatilde)
          }
          
          i2 <- i2+1
        }
      }
      
      i1 <- i1+1
      
    }
    
    
    if (OK1 && OK2) {
      lambda2.ini <- (lambdatilde - lambda1.ini)  / den
      alpha1.ini <-   (lambda2.ini - lambdatilde )   / (lambda2.ini - lambda1.ini)
      resest <- c(alpha1.ini, lambda1.ini, lambda2.ini)
      names(resest) <- c("prob1", "rate1", "rate2")
      resmeth <- "Hreg"
      res <- list(estimate = resest,
                  method = resmeth)
    } else {
      lambda1.ini <- 0.9/xbar
      lambda2.ini <- 1.1/xbar
      alpha1.ini <- 0.5
      resest <- c(alpha1.ini, lambda1.ini, lambda2.ini)
      names(resest) <- c("prob1", "rate1", "rate2")
      resmeth <- "arbitrary"
      res <- list(estimate = resest,
                  method = resmeth)
    } 
    
  }
      
  if (plot) {

    if( !needsmore) {
      n <- length(x)
      ## pct <- c(10, 90)
      xbar <- mean(x)
      xbarinv <- 1/xbar
      xs <- sort(x)
      H <- -log( (n+1-(1:n))/(n+1) )
    }
    
    plot(xs, H, type ="p", pch = 21, col = "black")

    Hs <- Hmixexp2(xs,
                   prob1 = res$estimate["prob1"],
                   rate1 = res$estimate["rate1"],
                   rate2 = res$estimate["rate2"])
  
    lines(xs, Hs, col = "purple", lwd =2)
    
    if (needsmore && OK2) {
      points(xs[ind.L], H[ind.L], pch = 21, bg = "red3")
      abline(b = lambdatilde, a = 0, col = "red3")
    }
    if (needsmore && OK1) {
      points(xs[ind.U], H[ind.U], pch = 21, bg = "SpringGreen3")
      abline(coef = coef1, col = "SpringGreen3")
    }
  }

  res 
  
}

##==========================================================================
## Expectation Maximization
##
##==========================================================================

EM.mixexp <- function(x, m = 2) {

  n <- length(x)

  Tau <- matrix(0, nrow = n, ncol = m)
  dx <- matrix(0, nrow = n, ncol = m)
  xx <- as.matrix(x)[ , rep(1, m)]

  ## Insiitalize groups
  if (m == 2) {

    Ini <- ini.mixexp2(x)$estimate
    alpha <- c(Ini["prob1"], 1 - Ini["prob1"])
    theta <- 1/Ini[c("rate1", "rate2")]

  } else {

    g <- cut(x, breaks = m)
    gi <- as.integer(g)
    for (j in 1:m) Tau[ gi == j , j] <- 1
    
    alpha <- apply(Tau, 2, mean)
    theta <- apply(xx*Tau, 2, sum) / apply(Tau, 2, sum) 
  }
     
  maxiter <- 1000

  ## Store results
    
  Theta <- matrix(NA, nrow = maxiter, ncol = m)
  Alpha <- matrix(NA, nrow = maxiter, ncol = m)
  logL     <- rep(NA, length.out = maxiter)

  iter <- 1
  
  cvg <- FALSE
  vec.new <- c(alpha[1], theta)
  
  while ( (iter <= maxiter) && (!cvg) )  {

    vec.old <- c(alpha[1], theta)
    
    for (j in 1:m) {
      dx[ , j] <- dexp(x, rate = 1/theta[j])*alpha[j]
    }
    
    d <- apply(dx, 1, sum)
    
    Tau <- dx / as.matrix(d)[ , rep(1, m)]

    alpha <- apply(Tau, 2, mean)
    theta <- apply(xx*Tau, 2, sum) / apply(Tau, 2, sum)
    
    Theta[iter, ] <- theta
    Alpha[iter, ] <- alpha
    
    if (m == 2)  logL[iter] <-
      sum(dmixexp2(x, prob1 = alpha[1], rate1 = 1/theta[1], rate2 = 1/theta[2], log = TRUE))
                                   
    vec.old <- vec.new
    vec.new <- c(alpha[1], theta)

    if ( max(abs(vec.new-vec.old) / vec.old) < 0.0005) {
      cvg <- TRUE
      logL <- logL[1:iter]
      Alpha <- Alpha[1:iter, ]
      Theta <- Theta[1:iter, ]
    }
    
    iter <- iter + 1
    
  }

  if (FALSE) matplot(Theta, type = "l", lwd = 2)

  estimate <- c(alpha, 1/theta)
  names(estimate) <- c(paste("prob", 1:m, sep = ""),
                       paste("rate", 1:m, sep = ""))
  
  list(estimate = estimate,
       logL = logL,
       Alpha = Alpha,
       Theta = Theta) 
  
}


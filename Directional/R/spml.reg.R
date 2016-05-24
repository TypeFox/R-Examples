################################
#### Projected multivariate normal for circular or angular regression
#### Tsagris Michail 1/2014
#### mtsagris@yahoo.gr
#### Presnell, Morrison and Littell (1998), JASA
################################


spml.reg <- function(y, x, rads = TRUE, xnew = NULL) {
  ## y is the angular dependent variable
  ## x contains the independent variable(s)
  ## xnew is some new data or the current ones
  ## pred is either TRUE (xnew is new data) or
  ## FALSE (xnew is the same as x)
  ## if the data are in degrees we transform them into radians
  if (rads == F)   y <- y/180 * pi
  u <- cbind( cos(y), sin(y) )  ## bring the data onto the circle
  n <- nrow(u)
  mat <- model.matrix(y~ ., as.data.frame(x) )
  x <- as.matrix(mat[1:n, ])
  XX <- solve( crossprod(x), t(x) )
  p <- ncol(x)
  funa <- function(beta) {
    mu <- x %*% beta
    tau <- diag( u %*% t(mu) )
    ell <-  -0.5 * sum( diag( tcrossprod( mu ) ) ) +
      sum( log( 1 + tau * pnorm(tau) / dnorm(tau) ) ) - n * log(2 * pi)
    ell
  }
  para <- as.vector( coef( lm(u ~ x[, -1]) ) )  ## starting values
  ### E-M algorithm is implemented below
  lik <- NULL
  B <- matrix(para, ncol = 2)
  lik[1] <- funa(B)
  mu <- x %*% B
  tau <- diag( u %*% t(mu) )
  psit <- tau + pnorm(tau) / ( dnorm(tau) + tau * pnorm(tau) )
  M <- diag(psit)
  B <- ( XX %*% M ) %*% u
  lik[2] <- funa(B)
  i <- 2
  while ( lik[i] - lik[i - 1] > 1e-06 ) {
    i <- i + 1
    mu <- x %*% B
    tau <- diag( u %*% t(mu) )
    psit <- tau + pnorm(tau) / ( dnorm(tau) + tau * pnorm(tau) )
    M <- diag(psit)
    B <- ( XX %*% M ) %*% u
    lik[i] <- funa(B)
  }
  loglik <- lik[i]
  mu <- x %*% B
  tau <- diag( u %*% t(mu) )
  psit <- tau + pnorm(tau)/( dnorm(tau) + tau * pnorm(tau) )
  psit2 <- diag( 2 - tau * pnorm(tau) / ( dnorm(tau) + tau * pnorm(tau) )  -
                   ( pnorm(tau)/( dnorm(tau) + tau * pnorm(tau) ) )^2  )
  C <- u[, 1]   ;   S <- u[, 2]
  A1 <-  - crossprod(x)
  A2 <-  t(x) %*% psit2
  s11 <-  A1 + A2 %*% ( tcrossprod(C) %*% x )
  s12 <-  ( A2 %*% C ) %*% ( t(S) %*% x )
  s21 <- t(s12)
  s22 <-  A1 + A2 %*% ( tcrossprod(S) %*% x )
  se1 <- cbind(s11, s12)
  se2 <- cbind(s21, s22)
  se <-  - rbind(se1, se2)  ## negative Hessian of the log-likelihood
  se <- solve(se)
  se <- sqrt( diag(se) )  ## standard errors of the coefficients
  seb <- matrix(se, ncol = 2)
  colnames(B) <- colnames(seb) <- c("Cosinus of y", "Sinus of y")
  if ( is.null( colnames(x) ) )  {
    rownames(B) <- c( "Intercept", paste("X", 1:c(p - 1), sep = "") )
    rownames(seb) <- c( "Intercept", paste("X", 1:c(p - 1), sep = "") )
  } else rownames(B) <- rownames(seb) <- colnames(x)
  ## rho is the correlation between the fitted and the observed values
  ## the fitted values are in radians
  if ( !is.null(xnew) ) {  ## predict new values?
    xnew <- cbind(1, xnew)
    xnew <- as.matrix(xnew)
    est <- xnew %*% B
    est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
  } else {
    est <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)
  }
  est <- as.vector(est)
  if (rads == FALSE)  est = est * 180 /pi
  list(beta = B, seb = seb, loglik = loglik, est = est)
}

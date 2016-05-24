rbpspline <- function(y, x, knots, degree, penalty = "gcv", ...){

  if ((degree %% 2) == 0)
    stop("''degree'' must be an odd number")

  if (!is.numeric(penalty) && (penalty != "cv") && (penalty != "gcv"))
    stop("'penalty' must be either a numeric value, either 'cv' or 'gcv'")

  if (is.numeric(penalty))
    cv.hat <- NA

  if (penalty == "cv"){
    cv.fit <- cv(y, x, knots, degree, plot = FALSE, ...)
    penalty <- cv.fit$penalty
    cv.hat <- cv.fit$cv
  }

  if (penalty == "gcv"){
    gcv.fit <- gcv(y, x, knots, degree, plot = FALSE, ...)
    penalty <- gcv.fit$penalty
    cv.hat <- gcv.fit$gcv
  }

  if (is.null(dim(x))){
    idx <- order(x)
    x <- x[idx]
    y <- y[idx]
  }

  dsgn.mat <- rb(x, degree = degree, knots = knots,
                 penalty = penalty)
  pen.mat <- dsgn.mat$pen.mat
  dsgn.mat <- dsgn.mat$dsgn.mat

  ##To compute the beta vector, we use the Demmler-Reinch
  ##orthogonalization
  chol.mat <- chol(crossprod(dsgn.mat))
  chol.mat.inv <- solve(chol.mat)
  svd.mat <- svd(t(chol.mat.inv) %*% pen.mat %*% chol.mat.inv)

  A <- dsgn.mat %*% chol.mat.inv %*% svd.mat$u
  b <- t(A) %*% y

  fitted.values <- A %*% (b / (1 + penalty^degree * svd.mat$d))
  smooth.mat <- A %*% diag(1 / (1 + penalty^degree * svd.mat$d)) %*% t(A)
  beta <- chol.mat.inv %*% svd.mat$u %*% (b / (1 + penalty^degree * svd.mat$d))
  
  df <- sum(diag(smooth.mat))
  res.df <- length(y) - 2 * sum(diag(smooth.mat)) +
    sum(diag(smooth.mat %*% t(smooth.mat)))
  rank <- sum(eigen(smooth.mat, only.values = TRUE)$values != 0)
  
  
  fitted <- list(x = x, obs = y, penalty = penalty, knots = knots,
                 degree = degree, y = fitted.values,
                 dsgn.mat = dsgn.mat, pen.mat = pen.mat, cv = cv.hat,
                 smooth.mat = smooth.mat, df = df, res.df = res.df,
                 rank = rank, beta = beta, call = match.call())

  class(fitted) <- "pspline"
  return(fitted)
}

cv <- function(y, x, knots, degree, plot = TRUE, n.points = 150, ...){

  dsgn.mat <- rb(x, degree = degree, knots = knots, penalty = NULL)
  pen.mat <- dsgn.mat$pen.mat
  dsgn.mat <- dsgn.mat$dsgn.mat

  ##To compute the beta vector, we use the Demmler-Reinch
  ##orthogonalization
  chol.mat <- chol(crossprod(dsgn.mat))
  chol.mat.inv <- solve(chol.mat)
  svd.mat <- svd(t(chol.mat.inv) %*% pen.mat %*% chol.mat.inv)

  A <- dsgn.mat %*% chol.mat.inv %*% svd.mat$u
  b <- t(A) %*% y
  
  obj.fun <- function(penalty){

    if (penalty < 0)
      return(1e6)

    fitted.values <- A %*% (b / (1 + penalty^degree * svd.mat$d))
    smooth.mat <- A %*% diag(1 / (1 + penalty^degree * svd.mat$d)) %*% t(A)
    
    cv <- sum(((y - fitted.values) / (1 - diag(smooth.mat)))^2)
    return(cv)
  }

  opt <- nlm(obj.fun, sd(y), ...)

  if (opt$code >= 3)
    warning("Optimization may have not succeeded")

  if (plot){
    cv.val <- NULL
    lambda.vals <- seq(0, 2 * opt$estimate, length = n.points)
    for (lambda in lambda.vals)
      cv.val <- c(cv.val, obj.fun(lambda))

    cv.val[cv.val == 1e6] <- NA
         
    plot(lambda.vals, cv.val, xlab = expression(lambda), ylab = "CV",
         type = "l")
  }
  
  return(list(cv = opt$minimum, penalty = opt$estimate, nlm.code = opt$code))    

}

gcv <- function(y, x, knots, degree, plot = TRUE, n.points = 150, ...){

  dsgn.mat <- rb(x, degree = degree, knots = knots, penalty = NULL)
  pen.mat <- dsgn.mat$pen.mat
  dsgn.mat <- dsgn.mat$dsgn.mat

  ##To compute the beta vector, we use the Demmler-Reinch
  ##orthogonalization
  chol.mat <- chol(crossprod(dsgn.mat))
  chol.mat.inv <- solve(chol.mat)
  svd.mat <- svd(t(chol.mat.inv) %*% pen.mat %*% chol.mat.inv)

  A <- dsgn.mat %*% chol.mat.inv %*% svd.mat$u
  b <- t(A) %*% y
  n <- length(y)
  
  obj.fun <- function(penalty){

    if (penalty < 0)
      return(1e6)

    fitted.values <- A %*% (b / (1 + penalty^degree * svd.mat$d))
    smooth.mat <- A %*% diag(1 / (1 + penalty^degree * svd.mat$d)) %*% t(A)
    
    gcv <- n^2 * sum((y - fitted.values)^2) / (n - sum(diag(smooth.mat)))^2
    return(gcv)
  }

  opt <- nlm(obj.fun, sd(y), ...)

  if (opt$code >= 3)
    warning("Optimization may have not succeeded")

  if (plot){
    gcv.val <- NULL
    lambda.vals <- seq(0, 2 * opt$estimate, length = n.points)
    for (lambda in lambda.vals)
      gcv.val <- c(gcv.val, obj.fun(lambda))

    gcv.val[gcv.val == 1e6] <- NA
         
    plot(lambda.vals, gcv.val, xlab = expression(lambda), ylab = "GCV",
         type = "l")
  }
  
  return(list(gcv = opt$minimum, penalty = opt$estimate, nlm.code = opt$code))    

}

rb <- function(..., knots, degree, penalty){

  data <- cbind(...)
    
  if ((degree %% 2) == 0)
    stop("''degree'' must be an odd number")
  
  if (missing(penalty))
    penalty <- NULL

  X0 <- NULL
  for (i in 1:((degree - 1) / 2))
    X0 <- cbind(X0, data^i)
  
  X0 <- cbind(1, X0)
  
  ##Define the number of ``purely parametric'' parameters to be
  ##estimated i.e. weights related to radial basis function are not
  ##taken into account
  n.ppar <- ncol(X0)
  
  X1 <- NULL

  if (is.null(dim(knots))){
    n.knots <- length(knots)

    for (k in 1:n.knots)
      X1 <- cbind(X1, as.matrix(dist(c(knots[k], data)))[-1,1])

  }

  else{
    n.knots <- nrow(knots)
    for (k in 1:n.knots)
      X1 <- cbind(X1, as.matrix(dist(rbind(knots[k,], data)))[-1,1])
  }   

  X1 <- X1^degree
  dsgn.mat <- cbind(X0, X1)

  pen.mat <- matrix(0, nrow = dim(dsgn.mat)[2], ncol = dim(dsgn.mat)[2])
  K <- as.matrix(dist(knots, upper = TRUE, diag = TRUE))^degree
  pen.mat[-(1:n.ppar),-(1:n.ppar)] <- t(sqrt(K)) %*% sqrt(K)

  dsgn.mat <- as.matrix(dsgn.mat)
  pen.mat <- as.matrix(pen.mat)

  ans <- list(dsgn.mat = dsgn.mat, pen.mat = pen.mat, degree = degree,
              penalty = penalty, knots = knots, data = data, n.ppar = n.ppar,
              call = match.call())

  return(ans)
}

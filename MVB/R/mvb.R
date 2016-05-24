# fit generalized linear regression
unifit <- function(formula, data = list(),
                   family = c("gaussian", "binomial"),
                   output = 0) {
  # prepare the generic arguments
  this.call <- match.call()
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])

  result <- .Call("unifit", x, y, family, as.integer(output))
  result$family <- family[1]
  result$call <- this.call
  class(result) <- c("mvbfit", class(result))
  result
}

# fit generalized linear regression with lasso penalty
unilps <- function(formula, data = list(),
                   family = c("gaussian", "binomial"),
                   lambda = NULL, nlambda = 100,
                   lambda.min.ratio = ifelse(nobs<nvars, .01, .0001),
                   output = 0, tune = c("AIC", "BIC", "GACV", "BGACV")) {
  # prepare the generic arguments
  this.call <- match.call()
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  nlambda <- as.integer(nlambda)

  # NOTES: reset lambda.max here
  lambda.max <- max(abs(t(x) %*% y / nobs) * 1.1)
  if (is.null(lambda)) {
    if (is.null(lambda.min.ratio)) 
      lambda.min.ratio = ifelse(nobs < nvars, .05, 1e-3)
    if (lambda.min.ratio >= 1)
      stop("lambda.min.ratio should be less than 1")
    wlambda <- exp( seq(log(as.double(lambda.max)),
                        log(as.double(lambda.min.ratio)),
                        log(as.double(lambda.min.ratio / lambda.max))
                        /nlambda) )
    wlambda <- wlambda[1:nlambda]
  } else {
    if (any(lambda < 0)) stop("lambda should be non-negative")
    wlambda <- as.double(rev(sort(lambda)))
    nlambda <- length(wlambda)
  }
  
  lambda_grid <- matrix(wlambda, nrow = 1)
  result <- .Call("unilps", x, y, lambda_grid, family, as.integer(output), tune)
  result$family <- family[1]
  result$call <- this.call
  result$lambda <- wlambda
  result$beta <- as.matrix(result$beta)
  # dim(result$beta) <- c(nvars, dim(result$beta)[1] / nvars)
  result$optbeta <- as.vector(result$beta[,which.min(result$score)])
  class(result) <- c("mvbfit", "lps", class(result))
  result
}

# fit multivariate bernoulli model
mvbfit <- function(x, y, maxOrder = 2,
                   output = 0, printIter = 100) {
  this.call <- match.call()
  x <- as.matrix(x)
  y <- as.matrix(y)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])

  result <- .Call("mvbfit", x, y, as.integer(maxOrder), as.integer(output),
                  as.integer(printIter))
  result$family <- "mvbernoulli"
  result$maxOrder <- maxOrder
  result$call <- this.call
  result$beta <- as.matrix(result$beta)
  dim(result$beta) <- c(nvars, dim(result$beta)[1] / nvars)
  class(result) <- c("mvbfit", class(result))
  result
}

# lasso pattern search for multivariate Bernoulli
mvblps <- function(x, y, maxOrder = 2, lambda = NULL, nlambda = 100,
                   lambda.min.ratio = ifelse(nobs<nvars, .01, .0001),
                   output = 0, printIter = 100, search = c('nm', 'grid'),
                   tune = c("AIC", "BIC", "GACV", "BGACV")) {
  this.call <- match.call()
  # a constant column is added to x
  x <- as.matrix(x)
  y <- as.matrix(y)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2]) + 1 # include constant
  numCol <- 0
  for (order in 1:maxOrder)
    numCol <- numCol + choose(ncol(y), order)
  nvars <- nvars * numCol
  nlambda <- as.integer(max(nlambda^(1/numCol), 2))
  x <- as.matrix(cbind(rep(1, nobs), x))

  # NOTES: reset lambda.max here
  lambda.max <- max(abs(t(x) %*% y / nobs) * 1.1)
  if (search[1] == "grid") {
    if (is.null(lambda)) {
      if (is.null(lambda.min.ratio)) 
        lambda.min.ratio = ifelse(nobs < nvars, .05, 1e-3)
      if (lambda.min.ratio >= 1)
        stop("lambda.min.ratio should be less than 1")
      wlambda <- exp( seq(log(as.double(lambda.max)),
                          log(as.double(lambda.min.ratio)),
                          log(as.double(lambda.min.ratio / lambda.max))
                          /nlambda) )
      wlambda <- wlambda[1:nlambda]
    } else {
      if (any(lambda < 0)) stop("lambda should be non-negative")
      wlambda <- as.double(rev(sort(lambda)))
      nlambda <- length(wlambda)
    }
    lambda_grid <- t(as.matrix(expand.grid(rep(list(wlambda), numCol)))[,numCol:1])
  } else {
    lambda_grid <- matrix(rep(lambda.max, numCol * 2), numCol, 2)
  }

  result <- .Call("mvblps", x, y, as.integer(maxOrder),
                  lambda_grid, as.integer(output),
                  as.integer(printIter), tune, search)
  result$family <- "mvbernoulli"
  result$maxOrder <- maxOrder
  result$call <- this.call
  result$lambda <- lambda_grid
  result$beta <- as.matrix(result$beta)
  if (search[1] == 'grid') {
    dim(result$beta) <- c(nvars, dim(result$beta)[1] / nvars)
    result$optbeta <- as.vector(result$beta[,which.min(result$score)])
  } else {
    result$optbeta = result$beta
  }
  dim(result$optbeta) <- c(ncol(x), length(result$optbeta) / ncol(x))
  result$tune <- tune[1]
  class(result) <- c("mvbfit", "lps", class(result))
  result
}

# step-wise fit multivariate Bernoulli model
# NOTES: implement backward elimination only (not consider order yet)
stepfit <- function(x, y, maxOrder = 2,
                    output = 0,
                    direction = c("backward", "forward"),
                    tune = c("AIC", "BIC", "GACV", "BGACV"),
                    start = NULL) {
  this.call <- match.call()
  x <- as.matrix(x)
  y <- as.matrix(y)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])

  if (is.null(start)) {
    obj = rep(0, 1)
  } else {
    maxOrder = start$maxOrder
    if ("lps" %in% class(start)) {
      obj = as.vector(start$optbeta)
      x <- as.matrix(cbind(rep(1, nobs), x))
      nvars <- nvars + 1
    } else {
      obj = as.vector(start$beta)
    }
  }
  result <- .Call("stepfit", x, y, as.integer(maxOrder),
                  direction, tune, as.integer(output), obj)
  result$family <- "mvbernoulli"
  result$call <- this.call
  result$beta <- as.matrix(result$beta)
  dim(result$beta) <- c(nvars, dim(result$beta)[1] / nvars)
  class(result) <- c("mvbfit", class(result))
  result
}

# convert decimal to other
integer.base.b <- function(x, b = 2) {
  xi <- as.integer(x)
  if(any(is.na(xi) | ((x-xi)!=0)))
    print(list(ERROR="x not integer", x=x))
  N <- length(x)
  xMax <- max(x)	
  ndigits <- (floor(logb(xMax, base=2))+1)
  Base.b <- array(NA, dim=c(N, ndigits))
  for(i in 1:ndigits){#i <- 1
    Base.b[, ndigits-i+1] <- (x %% b)
    x <- (x %/% b)
  }
  if(N ==1) Base.b[1, ] else Base.b
}

# sort vectors with bubble algorithm
isLess <- function(vec1, vec2) {
  if (sum(vec1) < sum(vec2)) 
    return(1)
  if (sum(vec1) > sum(vec2))
    return(0)
  for (k in 1:length(vec1)) {
    if (vec1[k] > vec2[k]) 
      return(1)
    if (vec1[k] < vec2[k]) 
      return(0)
  }
}
sort.vec <- function(mat) {
  for (i in 1:(nrow(mat) - 1)) {
    current <- i
    for (j in (i + 1):nrow(mat)) {
      # compare the rows
      if (isLess(mat[j,], mat[current,])) {
        current <- j
      }
    }
    tmp <- mat[i,]
    mat[i,] <- mat[current,]
    mat[current,] <- tmp
  }
  mat
}
       
       
# generate simulated data from bernoulli (for order 2 only)
mvb.simu <- function(coefficients, x, K = 2, offset = as.double(0)) {
  # assume design matrix x comes without constant term (unit column not included)
  p <- dim(coefficients)[1]
  givenCol <- dim(coefficients)[2]
  n <- dim(x)[1]
  beta <- matrix(0, p, 2^K - 1)
  beta[,1:givenCol] = coefficients
  if (length(offset) < ncol(beta))
    offset <- c(offset, rep(0, ncol(beta) - length(offset)) )
  fx <- cbind(rep(1, n), x) %*% rbind(offset[1:ncol(beta)], beta)
  eS <- .Call("get_eS", fx, as.integer(K))
  eb <- apply(eS, 1, sum)
  probs <- cbind(rep(1, n), eS) / eb
  prob <- cbind(rep(0, n), t(apply(probs, 1, cumsum)))
  unifrand <- runif(n, 0, 1)
  response <- apply(unifrand - prob, 1, function(vec){min(which(vec < 0) - 1)})
  category <- rep(0, K)
  for (num in 1:(2^K - 1) ) {
    binary <- integer.base.b(num)
    category <- rbind(category, c(rev(binary), rep(0, K - length(binary))) )
  }
  category <- sort.vec(category)
  res <- category[response,]
  rownames(res) <- NULL
  list(response = res, beta = beta)
}


# calculate negative log-likelihood for the distribution
loglike <- function(x, y, input,
                    family = c("gaussian", "bernoulli", "mvbernoulli")) {
  this.call <- match.call()
  x <- as.matrix(x)
  y <- as.matrix(y)
  input <- as.vector(input)
  if ( length(input) %% ncol(x) != 0)
    x <- cbind(rep(1, nrow(x)), x)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])

  .Call("loglike", x, y, input, family)
}

# fit mixed effects model
mvbme <- function(x, y, z, maxOrder = 2,
                  output = 0, printIter = 100) {
  this.call <- match.call()
  x <- as.matrix(x)
  y <- as.matrix(y)
  zz <- as.factor(z) # only one mixed effect allowed
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  numLevel <- length(levels(zz))
  z <- as.numeric(zz)

  result <- .Call("mvbme", x, y, z, as.integer(maxOrder), as.integer(output),
                  as.integer(printIter))
  result$family <- "mvbernoulli"
  result$maxOrder <- maxOrder
  result$call <- this.call
  result$beta <- as.matrix(result$beta)
  result$b <- as.matrix(result$b)
  dim(result$beta) <- c(nvars, dim(result$beta)[1] / nvars)
  dim(result$b) <- c(dim(result$b)[1] / dim(result$beta)[2], dim(result$beta)[2])
  result$numLevel <- numLevel
  class(result) <- c("mvbfit", "mvbme", class(result))
  result
}

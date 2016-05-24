wqs.est <-
function(y.train, x.train, z.train = NULL, y.valid = y.train, x.valid = x.train, z.valid = z.train, n.quantiles = 4, B = 100, b1.pos = TRUE){

  # Check training data
  check_train <- check_xyz(x.train, y.train, z.train)
  if(check_train[1])stop("x.train must be a matrix")
  if(check_train[2])stop("check dimensions of training data")
  if(check_train[3])stop("check dimensions of training data")

  # Check validation data
  check_valid <- check_xyz(x.valid, y.valid, z.valid)
  if(check_valid[1])stop("x.valid must be a matrix")
  if(check_valid[2])stop("check dimensions of validation data")
  if(check_valid[3])stop("check dimensions of validation data")

  # Check other inputs
  if(B < 2)stop("value of B must be at least 2")
  if(class(b1.pos)!= "logical")stop("b1.pos must be logical value of TRUE or FALSE")
  if(n.quantiles < 2 | n.quantiles > 10)stop("n.quantiles must be at least 2 and no greater than 10")

  c <- dim(x.train)[2] # number of components

  # calculate quantiles 
  q.train <- quantile.fn(x.train, n.quantiles)
  q.valid <- quantile.fn(x.valid, n.quantiles)

  # specify lower and upper bounds for b1 and the weights
  bounds <- specify.bounds(b1.pos, c)
  ineqLB <- bounds$ineqLB
  ineqUB <- bounds$ineqUB 

  # specify initial values
  init <- specify.init(z.train, y.train, b1.pos, c)

  # Estimate weights across bootstrap samples for WQS Regression
  result <- wqs_b.est(y.train, q.train, z.train, B, pars = init, fun = objfn.cont, eqfun = lincon, eqB = 1, ineqfun = ineq, ineqLB, ineqUB, LB = NULL, UB = NULL)

  wts.matrix <- result$wts.matrix

  # Calculate final weights for each subset using relative test statistic 
  weights <- teststat.fn(wts.matrix, result$test_stat)

  # Apply weights to validation data Set
  final <- wqs.fit(q.valid, z.valid, y.valid, weights)

  out <- list(q.train, q.valid, wts.matrix, weights, final$WQS, final$fit)
  names(out) <- c("q.train", "q.valid", "wts.matrix", "weights", "WQS", "fit")
  return(out)
}

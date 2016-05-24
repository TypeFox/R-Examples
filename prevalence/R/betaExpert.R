betaExpert <-
function(best, lower, upper, p = 0.95, method = "mode"){
  ## check presence
  if (missing(best))
    stop("'best' is missing")
  if (missing(lower) & missing(upper))
    stop("at least 'lower' or 'upper' must be specified")

  ## check input values: range(0,1)
  checkInput(best, "best", range = c(0, 1))
  checkInput(p, "p", range = c(0, 1))
  if (!missing(lower))
    checkInput(lower, "lower", range = c(0, 1), minEq = 0)
  if (!missing(upper))
    checkInput(upper, "upper", range = c(0, 1), maxEq = 1)

  ## check input values: order
  if (!missing(lower))
    if (lower > best) stop("'lower' cannot be greater than 'best'")
  if (!missing(upper))
    if (upper < best) stop("'upper' cannot be smaller than 'best'")
  if (!missing(lower) & !missing(upper)) # useless??
    if (lower > upper) stop("'lower' cannot be greater than 'upper'")

  ## functions to optimize ~ mode
  f_mode <-
  function(x, mode, p, target){
    return(
      sum(
        (qbeta(p = p,
               shape1 = x,
               shape2 = (x * (1 - mode) + 2 * mode - 1) / mode) -
         target) ^ 2
    ))
  }

  f_mode_zero <-
  function(x, p, target){
    return((qbeta(p = p, shape1 = 1, shape2 = x) - target) ^ 2)
  }

  f_mode_one <-
  function(x, p, target){
    return((qbeta(p = p, shape1 = x, shape2 = 1) - target) ^ 2)
  }

  ## functions to optimize ~ mean
  f_mean <-
  function(x, mean, p, target){
    return(
      sum(
        (qbeta(p = p,
               shape1 = x,
               shape2 = (x * (1 - mean)) / mean) -
         target) ^ 2
    ))
  }

  ## define 'target' and 'p'
  if (!missing(lower) & missing(upper)){
    target <- lower
    p <- 1 - p
  } else if (!missing(upper) & missing(lower)){
    target <- upper
  } else if (!missing(upper) & !missing(lower)){
    target <- c(lower, upper)
    p <- c(0, p) + (1 - p) / 2
  }

  ## derive a and b (=shape1 and shape2)
  if (method == "mode"){
    if (best == 0){
      a <- 1
      b <- optimize(f_mode_zero, c(0, 1000), p = p, target = target)$minimum
    } else if (best == 1) {
      a <- optimize(f_mode_one, c(0, 1000), p = p, target = target)$minimum
      b <- 1
    } else {
      a <- optimize(f_mode, c(0, 1000),
                    mode = best, p = p, target = target)$minimum
      b <- (a * (1 - best) + 2 * best - 1) / best
    }
  } else if (method == "mean"){
      a <- optimize(f_mean, c(0, 1000),
                    mean = best, p = p, target = target)$minimum
      b <- (a * (1 - best)) / best
  }

  ## create 'out' dataframe
  out <- list(alpha = a, beta = b)
  class(out) <- "betaExpert"

  ## return 'out'
  return(out)
}
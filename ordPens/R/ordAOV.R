ordAOV <- function(x, y, type = c("RLRT", "LRT"), nsim = 10000,
null.sample = NULL, ...){

  type <- match.arg(type)
  type <- switch(type, RLRT="RLRT", LRT="LRT")

  ## Check x and y
  if(!(is.matrix(x) | is.numeric(x)))
    stop("x has to be a matrix or numeric vector")

  tol <- .Machine$double.eps^0.5
  if(any(abs(x - round(x)) > tol) | any(x < 1))
      stop("x has to contain positive integers")

  if(length(ncol(y)) > 0 | !is.numeric(y))
    stop("y has to be a numeric vector")

  if(any(is.na(x)))
    stop("Missing values in x are not allowed")

  if(any(is.na(y)))
    stop("Missing values in y are not allowed")

  # one-factorial anova
  if (ncol(cbind(x)) == 1)
    ordAOV1(x = x, y = y, type = type, nsim = nsim, null.sample = null.sample,
    ...)
  
  # multi-factorial anova (main effects only)
  else
    ordAOV2(x = x, y = y, type = type, nsim = nsim, null.sample = null.sample,
    ...)
}

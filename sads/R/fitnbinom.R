fitnbinom <- function(x, trunc=0, start.value, ...){
  dots <- list(...)
  if (!is.null(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){ 
    muhat <- length(x)/(length(x) + mean(x))
    sizehat <- muhat*mean(x) 
  }
  else{
    sizehat <- start.value[[1]]
    muhat <- start.value[[2]]
  }
  if (is.null(trunc)){
    LL <- function(size, mu) -sum(dnbinom(x, size=size, mu=mu, log = TRUE))
  } else{
    LL <- function(size, mu) -sum(dtrunc("nbinom", x = x, coef = list(size=size, mu=mu), trunc = trunc, log = TRUE))
  }
  result <- do.call("mle2", c(list(LL, start = list(size = sizehat, mu = muhat), data = list(x = x)), dots))
  new("fitsad", result, sad="nbinom", distr = distr.depr, trunc = ifelse(is.null(trunc), NaN, trunc))
}

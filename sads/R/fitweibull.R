fitweibull <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)) {
	  f <- function(shape, x){
		  n <- length(x)
		  (1 / shape) + (1 / n) * sum(log(x)) - (sum((x ^ shape) * log(x)) / sum(x ^ shape))
	  }
	  ka <- uniroot(f, interval = c(0.0000001, 10), x=x)$root
	  theta <- ((1 / length(x)) * sum(x ^ ka)) ^ (1 / ka)
  } else{
    ka <- start.value[1]
    theta <-start.value[2]
  }
  if (missing(trunc)){
    LL <- function(shape, scale) -sum(dweibull(x, shape, scale, log = TRUE))
  } else {
    LL <- function(shape, scale) -sum(dtrunc("weibull", x = x, coef = list(shape, scale), trunc = trunc, log = TRUE))
  }  
  result <- do.call("mle2", c(list(LL, start = list(shape = ka, scale = theta), data = list(x = x)), dots))
  new("fitsad", result, sad="weibull", distr = distr.depr, trunc = ifelse(missing(trunc), NaN, trunc)) 
}

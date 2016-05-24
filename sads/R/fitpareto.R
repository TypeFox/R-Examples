fitpareto <- function(x, trunc, start.value, upper = 20, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    alpha <- length(x)/sum(log(x)-log(min(x)))
  } else{
    alpha <- start.value[1]
  }
  if (missing(trunc)){
    LL <- function(shape,scale) -sum(dpareto(x, shape, scale, log = TRUE))
  } else {
    LL <- function(shape, scale) -sum(dtrunc("pareto", x = x,
                                             coef = list(shape = shape, scale = scale), trunc = trunc, log = TRUE))
  }  
  result <- do.call("mle2", c(list(LL, start = list(shape = alpha, scale=min(x)), fixed=list(scale=min(x)),
                 data = list(x = x), method = "Brent", lower = 0, upper = upper), dots))
  if(abs(as.numeric(result@coef) - upper) < 0.001) 
    warning("mle equal to upper bound provided. \n Try value for the 'upper' argument")
  new("fitsad", result, sad="pareto", distr = distr.depr, trunc = ifelse(missing(trunc), NaN, trunc)) 
}

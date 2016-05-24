fitlnorm <- function(x, trunc, start.value, ...){
  dots <- list(...)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if(missing(start.value)){
    meanlog <- mean(log(x))
    sdlog <- sd(log(x))
  } else{
    meanlog <- start.value[1]
    sdlog <-start.value[2]
  }
  if (missing(trunc)){
    LL <- function(meanlog, sdlog) -sum(dlnorm(x, meanlog, sdlog, log = TRUE))
  } else {
    LL <- function(meanlog, sdlog) -sum(dtrunc("lnorm", x, coef = list(meanlog = meanlog, sdlog = sdlog), trunc = trunc, log = TRUE))
  }  
  result <- do.call("mle2", c(list(LL, start = list(meanlog = meanlog, sdlog = sdlog), data = list(x = x)), dots))
  new("fitsad", result, sad="lnorm", distr = distr.depr, trunc = ifelse(missing(trunc), NaN, trunc)) 
}

fitzipf <- function(x, N, trunc, start.value, upper = 20, ...){
  if(class(x)!="rad") rad.tab <- rad(x)
  else rad.tab <- x
  y <- rep(rad.tab$rank, rad.tab$abund)
  dots <- list(...)
  if (!missing(trunc)){
    if (min(y)<=trunc) stop("truncation point should be lower than the lowest rank")
  }
  if(missing(N)){
    N <- max(rad.tab$rank)
  }
  if(missing(start.value)){
    p <- rad.tab$abund/sum(rad.tab$abund)
    lzipf <- function(s, N) -s*log(1:N) - log(sum(1/(1:N)^s))
    opt.f <- function(s) sum((log(p) - lzipf(s, length(p)))^2)
    opt <- optimize(opt.f, c(0.5, length(p)))
    sss <- opt$minimum
  }
  else{
    sss <- start.value
  }
  if(missing(trunc)){
    LL <- function(N, s) -sum(dzipf(y, N=N, s=s, log = TRUE))
  }
  else{
    LL <- function(N, s) -sum(dtrunc("zipf", x = y, coef = list(N = N, s = s), trunc = trunc, log = TRUE))
  }
  result <- do.call("mle2", c(list(LL, start = list(s = sss), data = list(x = y), fixed=list(N=N), method = "Brent", lower = 0, upper = upper), dots))
  if(abs(as.numeric(result@coef) - upper) < 0.001)
    warning("mle equal to upper bound provided. \n Try increase value for the 'upper' argument")
  new("fitrad", result, rad="zipf", distr = distr.depr, trunc = ifelse(missing(trunc), NaN, trunc), rad.tab=rad.tab)
}

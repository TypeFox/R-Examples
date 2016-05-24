fitrbs <- function(x, trunc, ...){
	dots <- list(...)
  if(class(x)!="rad") rad.tab <- rad(x)
  else rad.tab <- x
  y <- rep(rad.tab$rank, rad.tab$abund)
  N <- sum(rad.tab$abund)
  S <- length(rad.tab$abund)
  if (!missing(trunc)){
    if (min(y)<=trunc) stop("truncation point should be lower than the lowest rank")
  }
  if(missing(N)){
    N <- max(rad.tab$rank)
  }
  if(missing(trunc)){
    LL <- function(N, S) -sum(drbs(y, N, S, log = TRUE))
  }
  else{
    LL <- function(N, S) -sum(dtrunc("rbs", x = y, coef = list(N = N, S = S), trunc = trunc, log = TRUE))
  }
  result <- do.call("mle2", c(list(LL, start = list(N=N, S = S), data = list(x = y), fixed=list(N=N, S=S), eval.only=T), dots))
  #BUGFIX: the show method on mle2 class always expects the "convergence" slot to be set.
  result@details$convergence = 0 
  new("fitrad", result, rad="rbs", distr = distr.depr, trunc = ifelse(missing(trunc), NaN, trunc), rad.tab=rad.tab)
}

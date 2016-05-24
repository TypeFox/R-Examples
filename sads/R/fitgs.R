fitgs <- function(x, trunc, ...){
	dots <-list(...)
  if(class(x)!="rad") rad.tab <- rad(x)
  else rad.tab <- x
  y <- rep(rad.tab$rank, rad.tab$abund)
  S <- length(rad.tab$abund)
  if (!missing(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
    else{
      LL <- function(S, k) -sum(dtrunc("gs", x = y, coef = list( k = k, S = S),
                                          trunc = trunc, log = TRUE))
      }
  }
  if (missing(trunc)){
    LL <- function(S, k) -sum(dgs(y, k, S, log = TRUE))
  }
  result <- do.call("mle2", c(list(LL, start = list(k=0.01), data = list(x = y), fixed=list(S=S), method = "Brent", lower = 1e-16, upper = 1-1e-16), dots))
  new("fitrad", result, rad = "gs", distr = distr.depr, trunc = ifelse(missing(trunc), NaN, trunc), rad.tab=rad.tab)
}

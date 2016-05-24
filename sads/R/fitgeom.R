fitgeom <- function(x, trunc = 0, start.value, ...){
	dots <- list(...)
  if (!is.null(trunc)){
    if (min(x)<=trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    phat <- 1/(mean(x))
  } else{
    phat <- start.value
  }
  if (is.null(trunc)){
    LL <- function(prob) -sum(dgeom(x, prob, log = TRUE))
  } else{
    LL <- function(prob) -sum(dtrunc("geom", x = x, coef = prob, trunc = trunc, log = TRUE))
  }
        ##result <- do.call("mle2", c(list(LL, start = list(prob = phat), data = list(x = x), method = "Brent", lower = 0, upper = 1), dots)) ## Brent method does not converge sometimes.
        result <- do.call("mle2", c(list(LL, start = list(prob = phat), data = list(x = x)), dots))  
  new("fitsad", result, sad = "geom", distr = distr.depr, trunc = ifelse(is.null(trunc), NaN, trunc))
}

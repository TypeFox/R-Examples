

rvci <- function(obj, interval=0.95, one.sided=FALSE, left=TRUE)
{
  # NAME
  #  rvci - Uncertainty (Credible) Interval
  # 
  if (length(interval)>1) {
    ci <- rvquantile(obj, probs=range(interval))
  } else if (one.sided) {
    q <- if (left) interval else (1-interval)
    ci <- rvquantile(obj, q)
  } else {
    lower <- (1-interval)/2
    upper <- (lower + interval)
    ci <- rvquantile(obj, c(lower,upper))
  }
  return(ci)
}


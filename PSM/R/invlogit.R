`invlogit` <-
function(y,xmin=0,xmax=1) {
  res <- exp(y)/(1+exp(y))*(xmax-xmin)+xmin
  idx <- (y >  500) & !is.nan(y) & !is.na(y)
  res[idx] <- xmax[idx]
  idx <- (y < -500) & !is.nan(y) & !is.na(y)
  res[idx] <- xmin[idx]
  res
}


Detrend <- function(x, demean=TRUE) {
  # detrend using the row means
  if (class(x) == "matrix" | class(x) == "data.frame") {
    xx <- rowMeans(x, na.rm=TRUE)
  } else {
    xx <- x
  }
  N <- length(xx)
  if (N == 1) {
    # if only one time step is given, detrending amounts to substracting the mean
    trnd <- x
  } else {
    # otherwise estimate a linear function of time and take the fitted values
    # as "the trend" 
    t <- 1:N
    lmod <- lm(xx~t)
    trnd <- drop(cbind(1, t) %*% coef(lmod))
  }
  # if demean is false, add the grand mean back to x minus trend
  m <- ifelse(demean, 0, mean(unlist(x), na.rm=TRUE))
  return(x - trnd + m)
}


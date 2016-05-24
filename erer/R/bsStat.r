bsStat <- function(y, two = NULL, digits = c(2, 2), 
  use = 'complete.obs', na.rm = TRUE, ...)
{
  digits <- rep(digits, 2)[1:2]    
  ff <- data.frame(matrix(data = 0, nrow = ncol(y), ncol = ncol(y) + 6))
  colnames(ff) <- c("name", "mean", "stde", "mini", "maxi", "obno", colnames(y))
  ff$name <- colnames(y)
  for (i in 1:ncol(y)){
    ff[i, "mean"] <- round(mean(y[, i], na.rm = na.rm), digits=digits[2])
    ff[i, "stde"] <- round(  sd(y[, i], na.rm = na.rm), digits=digits[2])    
    ff[i, "mini"] <- round( min(y[, i], na.rm = na.rm), digits=digits[2])    
    ff[i, "maxi"] <- round( max(y[, i], na.rm = na.rm), digits=digits[2])    
    ff[i, "obno"] <- length(na.omit(y[, i]))
    for (j in 1:ncol(y)) {
      ff[i, j+6] <- round(cor(y[,i], y[,j], use = use), digits=digits[1])
    }
  }
  fstat <- ff[, 1:6]
  corr <- ff[, -c(2:6)]

  if (is.null(two)) {
      if (ncol(y) < 11) { two <- FALSE } else { two <- TRUE }
  }
  if (two) {
      result <- list(fstat = fstat, corr = corr)
  } else {
      result <- ff
  }
  return(result)
}
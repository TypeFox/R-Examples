randomPermutation <- function(x) {
  if (!(inherits(x, "jags") && inherits(x$model, "BMMmodel"))) 
    stop("Use only with 'jags' objects with model of class 'BMMmodel'.")
  k <- x$model$data$k
  n <- dim(x$results)
  permutedIndex <- as.vector(t(apply(matrix(seq_len(n[1]*k), ncol = k), 1, sample, size = k)))
  variables <- x$variables
  dropIndex <- NULL
  for (i in seq_along(variables)) {
    name <- variables[i]
    ii <- grep(name, colnames(x$results))
    if (length(ii) == k) {
      dummy <- x$results[,ii]
      dummy <- dummy[permutedIndex]
      x$results[,ii] <- dummy
    }
    else if (length(ii) > 1) {
      x$results <- x$results[,-ii]
      dropIndex <- c(dropIndex, i)
    }
  }
  if (length(dropIndex)) {
    x$variables <- x$variables[-dropIndex]
    warning("Variables have been dropped.")
  }
  x
}

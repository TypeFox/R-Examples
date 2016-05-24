dummy <- function (v){ #dummy coding with NA
  t <- table(v)
  lt <- length(t)
  n.obs <- length(v)
  new <- matrix(NA, nrow = n.obs, ncol = lt)
  vv <- as.numeric(as.factor(v))
  for (i in 1:n.obs) {
    new[i, vv[i]] <- 1
    new[i, -vv[i]] <- 0
  }
  colnames(new) <- names(t)
  new
  return(new)
}
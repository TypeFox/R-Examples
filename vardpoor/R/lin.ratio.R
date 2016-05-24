
lin.ratio <- function(Y, Z, weight, Dom = NULL, percentratio = 1) {
  
  if (length(percentratio) != 1 | !any(is.integer(percentratio) | percentratio > 0)) stop("'percentratio' must be the positive integer value")
 
  # Y
  Y <- data.table(Y, check.names = TRUE) 
  if (!all(sapply(Y, is.numeric))) stop("'Y' must be numerical")

  if (any(is.na(Y))) stop("'Y' has unknown values")
  if (!is.null(Dom)) Yd <- domain(Y, Dom) else Yd <- Y

  # Z
  Z <- data.table(Z, check.names = TRUE)
  if (nrow(Z)!= nrow(Y)) stop("'Y' and 'Z' have different row count")
  if (ncol(Z)!= ncol(Y)) stop("'Y' and 'Z' have different column count") 
  if (!all(sapply(Z, is.numeric))) stop("'Z' must be numerical")
  if (any(is.na(Z))) stop("'Z' has unknown values")
  if (!is.null(Dom)) Zd <- domain(Z, Dom) else Zd <- Z

  # weight
  weight <- data.frame(weight)
  if (nrow(weight) != nrow(Y)) stop("'weight' length is not equal with 'Y' row count")
  if (ncol(weight) != 1) stop("'weight' must be a vector or 1 column data frame, data matrix, data table")
  weight <- weight[,1]
  if (!is.numeric(weight)) stop("'weight' must be numerical")
  if (any(is.na(weight))) stop("'weight' has unknown values")
 
  Y_est <- colSums(Yd * weight)
  Z_est <- colSums(Zd * weight)
  R_est <- Y_est / Z_est

  percentratio <- as.integer(percentratio)
  U <- percentratio * t((1 / Z_est) * t(Yd - t(R_est * t(Zd))))
  return(data.table(U))
}
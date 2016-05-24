#' Creating stratified cross validation folds
#'
#' This function creates stratified cross validation folds, in case of unbalanced case-control groups.
#' @param y binary vector with group assignment
#' @param k Default k=10. Number of folds
#' @param stratified Defaults to TRUE. Stratification TRUE/FALSE
#' @keywords internal
#' @export

crossvalFolds <- function(y, k=10, stratified=T) {
  if (stratified) {
    t <- table(y)
    f <- matrix(nrow=length(t), ncol=k)
    for (i in 1:length(t)) {
      f[i,] <- ipred::kfoldcv(k, t[i])
    }
    rownames(f) <- names(t)
    idx <- rep(NA, length(y))
    data <- data.frame(y=y, x=1:length(y))
    for (i in 1:(k-1)) {
      s <- sampling::strata(data, "y", f[match(unique(data$y), rownames(f)), i], "srswor") #sample from y without replacement, taking as many observations as indicated by f
      sy <- table(data[s$ID, "y"])
      stopifnot(sy[rownames(f)] == f[,i]) # stop if the sampling was not according to f[,i]
      idx[data[s$ID, "x"]] <- i
      stopifnot(sum(!is.na(idx)) == sum(f[, 1:i]))
      data <- data[-s$ID,] # sampling without replacement
    }
    idx[data[, "x"]] <- k
    stopifnot(all(!is.na(idx)))
    return (idx)     # idx indicates which observations have been chosen for which outer cross-validation loop
  } else {
    folds.size <- ipred::kfoldcv(k, length(y))  #if no stratification, divide observations in k groups, without consideration of 0/1
    idx <- c()
    for (i in 1:k) {
      idx <- c(idx, rep(i, folds.size[i]))
    }
    idx <- sample(idx)
    return (idx)
  }
}

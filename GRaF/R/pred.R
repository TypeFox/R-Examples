pred <-
function(predx, fit, mn, std = TRUE, maxn = 250, same = FALSE) {
  predx <- as.matrix(predx)
  n <- nrow(predx)  
  if (n > maxn & !same) {
    inds <- split(1:n, ceiling((1:n) / maxn))
    fun <- function(ind, X, fit, std, maxn) {
      pred(X[ind, , drop = FALSE], fit, mn[ind], std, maxn, same)
    }    
    prediction <- lapply(inds, fun, predx, fit, std, maxn)
    prediction <- do.call('rbind', prediction)
  } else {
    if(same) {
      # if predicting back to input data re-use covariance matrix
      Kx <- fit$K
      prediction <- fit$MAP
    } else {
      Kx <- cov.SE(fit$x, predx, e1 = fit$e, l = fit$ls)
      mpred <- qnorm(mn)
      prediction <- crossprod(Kx, fit$a) + mpred
    }

    if (std) {
      v <- backsolve(fit$L, sqrt(as.vector(fit$W)) * Kx, transpose = T)
      # using correlation matrix, so diag(kxx) is all 1s, no need to compute kxx
      predvar <- 1 - crossprod(v)
      prediction <- cbind(prediction, sqrt(diag(predvar)))
      colnames(prediction) <- c("MAP", "std")
    }
  }
  prediction
}
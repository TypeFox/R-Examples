diag.disc <-function(x, centroids, prior, weight) {
### Computes the class discriminant functions assuming scaled x and centroids
  if(! missing(weight)) {
    posid <- (weight > 0)
    if(any(posid)) {
      weight <- sqrt(weight[posid])
      centroids <- centroids[posid,  , drop = FALSE] * weight
      x <- x[posid,  , drop = FALSE] * weight
    }
    else {
      mat <- outer(rep(1, ncol(x)), log(prior), "*")
      dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
      return(mat)
    }
  }
  dd <- t(x) %*% centroids
  dd0 <- drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  names(dd0) <- NULL
  scale(dd, dd0, FALSE)
}

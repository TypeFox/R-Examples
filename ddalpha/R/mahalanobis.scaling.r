
mah.moment <- function(x){
  mu <- colMeans(x)
  
  scale.eig <- eigen(cov(x))
  B <- scale.eig$vectors %*% diag(sqrt(scale.eig$values)) 
  B_inv <- solve(B)
  return (list(mu = as.numeric(mu), b = B_inv, s = cov(x)))
}

mah.mcd <- function(x, alpha = 1/2){
  #library(robustbase)
  estimate <- covMcd(x, alpha = alpha)
  mu <- estimate$center
  
  scale.eig <- eigen(estimate$cov)
  B <- scale.eig$vectors %*% diag(sqrt(scale.eig$values)) 
  B_inv <- solve(B)
  return (list(mu = as.numeric(mu), b = B_inv, s = estimate$cov))
}

mah.transform <- function (x, mu, B_inv, inv = F)  {
  if (inv)
    return (t(solve(B_inv) %*% (t(x))+mu))
  return (t(B_inv %*% (t(x)-mu)))
}

mah.transform.back <- function (x, mu, B_inv)  {
  return (t(solve(B_inv) %*% (t(x))+mu))
}

MahMomentTransformer <- function(mu, b){
  f <- function(points, inv = F){
    return(mah.transform(points, mu, b, inv))
  }
  environment(f) <- new.env()
  environment(f)$mu = mu
  environment(f)$b = b
  return (f)
}

# mahalanobisRegionsX <- function(x, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), col = "black"){
#   c <- c(0,0)
#   for (i in 1:nrow(x)){
#     c <- c + x[i,]
#   }
#   mu <- c / nrow(x)
#   
#   mahalanobisRegions(mu, cov(x), col)
# }
# 
# mahalanobisRegions <- function(mu, s, depths = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), col = "black"){
#   sigma.inv <- solve(s)
#   for (i in 1:length(depths)){
#     ellipsem(mu = mu, amat = sigma.inv, c2 = 1/depths[i] - 1, showcentre = F, col = col)
#   }
# }

depth.Mahalanobis <- function(x, data, mah.estimate = "moment", mah.parMcd = 0.75){
  if (mah.estimate == "moment") s = solve(cov(data))
  else if (mah.estimate == "MCD") s = solve(covMcd(data, alpha = mah.parMcd)$cov)
  else stop ("Wrong parameter 'estimate'")
  
  depths <- .Mahalanobis_depth(x, center = colMeans(data), sigma = s)
  return (depths)
}

depth.space.Mahalanobis <- function(data, cardinalities, mah.estimate = "moment", mah.parMcd = 0.75){
  if (is.data.frame(data))
    data <- data.matrix(data)
  if (!is.numeric(data)
      || !is.matrix(data) 
      || ncol(data) < 2){
    stop("Argument \"data\" should be a numeric matrix of at least 2-dimensional data")
  }
  if (!is.vector(cardinalities, mode = "numeric") 
      || is.na(min(cardinalities)) 
      || sum(.is.wholenumber(cardinalities)) != length(cardinalities) 
      || min(cardinalities) <= 0 
      || sum(cardinalities) != nrow(data)){
    stop("Argument \"cardinalities\" should be a vector of cardinalities of the classes in \"data\" ")
  }
  if (sum(cardinalities < ncol(data) + 1) != 0){
    stop("Not in all classes sufficiently enough objetcs")
  }
  
  depth.space <- NULL
  for (i in 1:length(cardinalities)){
    pattern <- data[(1 + sum(cardinalities[0:(i - 1)])):sum(cardinalities[1:i]),]
    pattern.depths <- depth.Mahalanobis (data, pattern, mah.estimate, mah.parMcd)
    depth.space <- cbind(depth.space, pattern.depths, deparse.level = 0)
  }
  
  return (depth.space)
}

gdf <- function(y) {
  a <- matrix(y,nrow=1)
  n <- length(y)
  A <- MW(n)
  A%*%t(a)
}
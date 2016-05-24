#' To determine the Optimum value of Cluster size (K) based on Average silhouette method
#'
#' Takes data and returns value of optimum K
#' @param data_in as Input data, in any format (data matrix data frame or vector). All variables should be numeric and NA values will get removed while execution.
#' @return K value, an integer as optimum number clusters.
#' @export
#' @importFrom cluster silhouette
#' @importFrom  stats dist kmeans na.omit
#' @examples
#' # iris dataset as input
#' x <- optimum_k(iris[2])

optimum_k <- function(data_in)
{
  #library(cluster)
  data_in <- na.omit(data_in)
  #data_in <- data.frame(data_in)
  #--- to avoid the Error "more cluster centers than distinct data points."
  # in "kmeans", k.max is decides as follows
  c <- table(data_in)

  k.max1 <- length(c)
  k.max <- min(10,k.max1)
  avg <- rep(0, k.max)
  for(i in 2:k.max){
    km.res <- kmeans(data_in, centers = i, nstart = 25)
    ss <- silhouette(km.res$cluster, dist(data_in))
    avg[i] <- mean(ss[, 3])
  }
  options(warn=-1)
  return(which.max(avg))
}

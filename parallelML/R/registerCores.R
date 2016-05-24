#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @export
registerCores <- function(numberCores){
  cluster <- makeCluster(numberCores)
  registerDoParallel(cluster)
}
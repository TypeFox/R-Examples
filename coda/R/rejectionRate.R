rejectionRate.mcmc <- function (x) {
  x <- as.matrix(x)
  apply(x[-nrow(x),,drop=FALSE] == x[-1,, drop=FALSE],2,mean)
}

rejectionRate.mcmc.list <- function (x) {
  apply(sapply(x,rejectionRate.mcmc),1,mean)
}

rejectionRate <- function(x) {
  UseMethod("rejectionRate")
}



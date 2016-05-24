marginal <- function(x, grid, d){
  fit <- x
  addfunction <- fit$func
  addfunction.var <- fit$index
  theta <- fit$theta
  rho <- fit$rho
  bigfm <- function(x, d){
    temp <- 0
    for (m in seq(1, length(theta))) {
      if (addfunction.var[[m]] == d) {
        temp <- temp + rho * theta[[m]] * predict(addfunction[[m]], x)$y
      }
    }
    return(temp)
  }
  return(bigfm(grid, d))
}



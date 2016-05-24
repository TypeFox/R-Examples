# ------------------------------------------------------------------------------
# Internal function 'localenv.get'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
localenv.get <- function(sprel, data, power, useExp, maxdist, tol) {
  
  grps <- colnames(data); rows <- rownames(data)
  
  if (class(sprel) == "nb") {
    xmat <- spdep::nb2mat(sprel, style = "W") # needs the package 'spdep'
    if (nrow(xmat) != nrow(data)) 
      stop("'data' must have the same number of rows as 'sprel'", call. = FALSE)
    env <- xmat %*% data
  }
  
  else if (class(sprel) == "dist") {
    sprel <- as.matrix(sprel)
    if (nrow(sprel) != nrow(data)) 
      stop("'data' must have the same number of rows as 'sprel'", call. = FALSE)
    env <- matrix(nrow = nrow(data), ncol = ncol(data))
    for (i in 1:nrow(data)) {
      if (useExp)
        weight <- exp(power * sprel[i,] * -1)
      else
        weight <- 1/(sprel[i,] + tol)^power
      if (maxdist >= 0)
        weight[which(sprel[i,] > maxdist)] <- 0
      env[i,] <- apply(data, 2, function(z) sum(z * weight)/sum(weight))
    }
  } 
  
  else {
    if (nrow(sprel) != nrow(data))
      stop("'data' must have the same number of rows as 'sprel'", call. = FALSE)
    xval <- sprel[,1]; yval <- sprel[,2]
    dim <- ncol(data); data <- as.vector(data)
    env <- .Call("envconstruct", xval, yval, data, as.integer(dim), power, 
                 as.integer(useExp), maxdist, tol)
    # --------------------------------------------------------------------------
    # R version 'envconstruct()'
    # --------------------------------------------------------------------------
    # envconstruct <- function(x, y, v, d) {
    #   n <- length(v)
    #   env <- rep(0, n)
    #   
    #   for (i in 1:n) {
    #     m <- 0
    #     for (j in 1:n) {
    #       dx <- x[i] - x[j]
    #       dy <- y[i] - y[j]
    #       if (dx <= d && dy <= d) {
    #         if (dx^2 + dy^2 <= d^2) {
    #           env[i] <- env[i] + v[j]
    #           m <- m + 1
    #         }
    #       }
    #     }
    #     if (m > 1)
    #       env[i] <- env[i] / m
    #   }
    #   return(env)
    # }
    # --------------------------------------------------------------------------
  }
  
  colnames(env) <- grps; rownames(env) <- rows
  env
}

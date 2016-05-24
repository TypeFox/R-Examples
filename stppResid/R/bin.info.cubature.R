bin.info.cubature <- function(X, cifunction, theta = NULL, gf, type = c(1,2), tol, maxEval, absError)
{
  if(is.null(theta)) {
    if(type == 1) {
      fun2 <- function(z) 
      {
        Y <- rbind(data, z)
        Y <- Y[which(Y$t <= z[3]),]
        Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
        tail(cifunction(Y),1)
      }
    } else {
      fun2 <- function(z) 
      {
        Y <- rbind(data, z)
        Y <- Y[which(Y$t <= z[3]),]
        Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
        sqrt(tail(cifunction(Y),1))
      }
    }
    compute.fun <- function(node) 
    {
      lower.limit <- c(node[c(1,3)], X$tcoord[1])
      upper.limit <- c(node[c(2,4)], X$tcoord[2])
      integral <- adaptIntegrate(fun2, lower.limit, upper.limit, tol = tol, maxEval = maxEval, absError = absError)
    }
  }
  if(!is.null(theta)) {
    if(type == 1) {
      fun2 <- function(z) 
      {
        Y <- rbind(data, z)
        Y <- Y[which(Y$t <= z[3]),]
        Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
        tail(cifunction(Y, theta),1)
      }
    } else {
      fun2 <- function(z) 
      {
        Y <- rbind(data, z)
        Y <- Y[which(Y$t <= z[3]),]
        Y <- stpp(Y$x, Y$y, Y$t, stwin(X$xcoord, X$ycoord, X$tcoord))
        sqrt(tail(cifunction(Y, theta),1))
      }
    }
    compute.fun <- function(node) 
    {
      lower.limit <- c(node[c(1,3)], X$tcoord[1])
      upper.limit <- c(node[c(2,4)], X$tcoord[2])
      integral <- adaptIntegrate(fun2, lower.limit, upper.limit, tol = tol, maxEval = maxEval, absError = absError)
    }
  }
  data <- data.frame(cbind(X$x, X$y, X$t))
  names(data) <- c("x", "y", "t")
  temp.grid <- gf$grid.full
  info <- apply(temp.grid, 1, compute.fun)
  info <- matrix(unlist(info), ncol = 4, byrow = TRUE)
  int.approx <- info[,1]
  lamb.sd.final <- info[,2]
  n.final <- info[,3]
  bins <- list(n = n.final, integral = int.approx, sd.lambda = lamb.sd.final)
  class(bins) <- "bin.info"
  return(bins)
}

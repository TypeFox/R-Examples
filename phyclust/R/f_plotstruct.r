### This file contains the function for the structure plot.

plotstruct <- function(Z, X.class = NULL, sort.inside.class = TRUE,
    direction = "h", main = "Structure Plot", xlab = "Observations",
    ylab = "Posterior Probabilities", ...){
  if(direction != "h" && direction != "v"){
    stop("The direction is not found.")
  }

  N <- nrow(Z)
  K <- ncol(Z)
  color <- .Color[0:(K - 1) %% length(.Color) + 1]

  if(is.null(X.class)){
    class.id <- apply(Z, 1, which.max)
  } else{
    class.id <- X.class
  }

  id.X <- order(class.id)
  class.id <- class.id[id.X]
  Z <- Z[id.X,]

  if(sort.inside.class && is.null(X.class)){
    for(k in 1:K){
      tmp.Z <- Z[class.id == k, ]
      id.X <- order(tmp.Z[, k], decreasing = TRUE)
      Z[class.id == k,] <- tmp.Z[id.X,]
    }
  }

  Z.cum <- apply(Z, 1, cumsum)
  n.class <- as.vector(table(class.id))
  X.nc.cum <- cumsum(n.class) + 1

  if(direction == "h"){
    xlim <- c(1, N + 1)
    ylim <- c(0, 1)
    plot(NULL, NULL, type = "n", xlim = xlim, ylim = ylim,
         main = main, xlab = xlab, ylab = ylab, ...)
    for(n in 1:N){
      x.left <- n
      y.bottom <- 0
      for(k in 1:K){
        rect(x.left, y.bottom, x.left + 1, Z.cum[k, n],
             col = color[k], border = NA)
        y.bottom <- Z.cum[k, n]
      }
    }
    abline(v = X.nc.cum[-K], lty = 3, lwd = 0.5)
  } else if(direction == "v"){
    xlim <- c(0, 1)
    ylim <- c(N + 1, 1)
    plot(NULL, NULL, type = "n", xlim = xlim, ylim = ylim,
         main = main, xlab = ylab, ylab = xlab, ...)
    for(n in 1:N){
      y.bottom <- n
      x.left <- 0
      for(k in 1:K){
        rect(x.left, y.bottom, Z.cum[k, n], y.bottom + 1,
             col = color[k], border = NA)
        x.left <- Z.cum[k, n]
      }
    }
    abline(h = X.nc.cum[-K], lty = 3, lwd = 0.5)
  }
} # End of plotstruct().


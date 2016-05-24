laymanmetrics <- function(x,y){

  out <- list()

  # --------------------------------------
  # Layman metric # 1 - dN range
  out$dN_range <- max(y) - min(y)

  # --------------------------------------
  # Layman metric # 2 - dC range
  out$dC_range <- max(x) - min(x)

  # --------------------------------------
  # Layman metric #3 - area of convex hull
  # some convex hull stuff
  out$hull <- convexhull(x,y)

  # --------------------------------------
  # Layman metric # 4 - mean distance to centroid CD
  mean_y <- mean(y)
  mean_x <- mean(x)

  out$CD <- mean( ( (mean_x - x)^2 + (mean_y - y)^2 ) ^ 0.5 )

  # --------------------------------------
  # Layman metric # 5 - mean nearest neighbour distance NND
  NNDs <- numeric(length(x))
  for (j in 1:length(x)){
    tmp <- ( (x[j] - x)^2 + (y[j] - y)^2 ) ^ 0.5
    tmp[j] <- max(tmp)
    NNDs[j] <-   min(tmp)
  }

  out$MNND <- mean(NNDs)

  # --------------------------------------
  # Layman metric # 6 - standard deviation of nearest neighbour distance SDNND
  out$SDNND <- sd(NNDs)

  out

}
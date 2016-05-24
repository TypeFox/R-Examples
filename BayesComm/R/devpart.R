devpart <-
function (null, environment, community, full) {
  getdev <- function (y, z) {
    p <- pnorm(z, log.p = TRUE)
    for(j in 1:ncol(y)) {
      ind <- which(!y[, j])
      p[, ind, j] <- log(-expm1(p[, ind, j]))
    }
    -2 * apply(p, 1, colSums)
  }
  nuldev <- getdev(null$call$Y, null$trace$z)
  envdev <- getdev(environment$call$Y, environment$trace$z)
  comdev <- getdev(community$call$Y, community$trace$z)
  fuldev <- getdev(full$call$Y, full$trace$z)
  
  nul <- rbind(rowMeans(nuldev), apply(nuldev, 1, quantile, c(0.25, 0.975)))
  env <- rbind(rowMeans(envdev), apply(envdev, 1, quantile, c(0.25, 0.975)))
  com <- rbind(rowMeans(comdev), apply(comdev, 1, quantile, c(0.25, 0.975)))
  ful <- rbind(rowMeans(fuldev), apply(fuldev, 1, quantile, c(0.25, 0.975)))
  rownames(nul)[1] <- rownames(env)[1] <- rownames(com)[1] <- rownames(ful)[1] <- "Mean"

  propR2 <- t(rbind(1 - env[1, ] / nul[1, ], 1 - com[1, ] / nul[1, ], 1 - ful[1, ] / nul[1, ], rep(1, length(env[1, ]))))
  colnames(propR2) <- c("env", "com", "full", "total")
  rownames(propR2) <- colnames(full$call$Y)
  list(devpart = propR2, null = nul, environment = env, community = com, full = ful)
}
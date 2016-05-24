conout <- function(y, trtment, trt_ini, subgrp, sbp_lev, nsbp) {
  
  betaout <- matrix(0, nsbp, 2)
  n <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))
  mean <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))
  std <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))

  for(i in 1 : nsbp) {
    
    x <- sbp_lev[i]
    y1 <- y[which((subgrp == x) & (trtment == trt_ini[1]) & !is.na(y))]
    n$ref[i] <- length(y1)
    mean$ref[i] <- mean(y1)
    std$ref[i] <- sd(y1)
    y2 <- y[which((subgrp == x) & (trtment == trt_ini[2]) & !is.na(y))]
    n$trt[i] <- length(y2)
    mean$trt[i] <- mean(y2)
    std$trt[i] <- sd(y2)
    betaout[i, ] <- c(mean$trt[i] - mean$ref[i], 
                      sqrt(std$ref[i] ^ 2 / n$ref[i] + std$trt[i] ^ 2 / n$trt[i]))
    
  }
  
  out <- list(n = n, trtment = trt_ini, trtname = trt_ini[1], 
              nsbp = nsbp, subgroup = sbp_lev, beta = betaout)
  out
  
}
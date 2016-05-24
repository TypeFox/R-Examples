binout <- function(y, trtment, trt_ini, subgrp, sbp_lev, nsbp, scale) {
  
  betaout <- matrix(0, nsbp, 2)
  n <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))
  risk <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))
  event <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))
  
  for(i in 1 : nsbp) {
    
    x <- sbp_lev[i]
    y1 <- y[which((subgrp == x) & (trtment == trt_ini[1]) & !is.na(y))]
    event$ref[i] <- sum(y1 == 1)
    n$ref[i] <- length(y1)
    risk$ref[i] <- event$ref[i] / n$ref[i]
    y2 <- y[which((subgrp == x) & (trtment == trt_ini[2]) & !is.na(y))]
    event$trt[i] <- sum(y2 == 1)
    n$trt[i] <- length(y2)
    risk$trt[i] <- event$trt[i] / n$trt[i]
    
    betaout[i, ] <- switch(scale, 
                           RD = c(risk$trt[i] - risk$ref[i], 
                                  sqrt(risk$trt[i] * (1 - risk$trt[i]) / n$trt[i] + 
                                         risk$ref[i] * (1 - risk$ref[i]) / n$ref[i])),
                           RR = c(log((event$trt[i] + 0.5) / (n$trt[i] + 0.5)) - 
                                    log((event$ref[i] + 0.5) / (n$ref[i] + 0.5)), 
                                  sqrt(1 / (event$ref[i] + 0.5) + 1 / (event$trt[i] + 0.5) - 
                                         1 / (n$ref[i] + 0.5) - 1 / (n$trt[i] + 0.5))),
                           OR = c(log((event$trt[i] + 0.5) * (n$ref[i] - event$ref[i] + 0.5) / 
                                        (event$ref[i] + 0.5)/ (n$trt[i] - event$trt[i] + 0.5)), 
                                  sqrt(1 / event$ref[i] + 1 / (n$ref[i] - event$ref[i]) + 
                                         1 / event$trt[i] + 1 / (n$trt[i] - event$trt[i]))))  
    
  }
  
  out <- list(n = n, trtment = trt_ini, trtname = trt_ini[1], 
              nsbp = nsbp, subgroup = sbp_lev, beta = betaout)
  out
  
}


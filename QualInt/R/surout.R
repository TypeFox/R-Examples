#' @importFrom survival coxph

surout <- function(y, trtment, trt_ini, subgrp, sbp_lev, nsbp) {
  
  betaout <- matrix(0, nsbp, 2)
  n <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))
  rate <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))
  event <- list(ref = rep(0, nsbp), trt = rep(0, nsbp))
  
  outmod <- lapply(sbp_lev, function(x) coxph(y ~ trtment, subset = (subgrp == x)))
  
  mod <- outmod[[2]]
  trttar <- substring(tail(names(mod$coef), n = 1), 8)
  trtname <- trt_ini[trt_ini != trttar]
  
  for(i in 1 : nsbp) { 
    
    x <- sbp_lev[i]
    
    y1 <- y[which((subgrp == x) & (trtment == trtname) & !is.na(y))]
    n$ref[i] <- dim(y1)[1]
    y2 <- y[which((subgrp == x) & (trtment == trttar) & !is.na(y))]
    n$trt[i] <- dim(y2)[1]
 
    modd <- outmod[[i]]
    trttarr <- substring(tail(names(modd$coef), n = 1), 8)
    betaout[i, ] <- coef(summary(modd))[, c("coef", "se(coef)")]
    if(trttarr == trtname) betaout[i, 1] <- -betaout[i, 1]
    
  }
  
  out <- list(n = n, trtment = trt_ini, trtname = trtname, 
              nsbp = nsbp, subgroup = sbp_lev, beta = betaout)
  out
  
}
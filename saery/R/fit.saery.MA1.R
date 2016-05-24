fit.saery.MA1 <-
function(X, ydi, D, md, sigma2edi, conf.level){
  
  fit <- REML.saery.MA1(X, ydi, D, md, sigma2edi)
  sigmau1.hat <- fit[[1]][1]
  sigmau2.hat <- fit[[1]][2]
  theta.hat  <- fit[[1]][3]
  stderr <- stderr.saery.MA1(fit)
  alpha <- 1-conf.level
  k <- 1-alpha/2
  z <- qnorm(k)
  
  lower.s <- c(sigmau1.hat-z*stderr[[2]][1], sigmau2.hat-z*stderr[[2]][2], theta.hat-z*stderr[[2]][3])
  upper.s <- c(sigmau1.hat+z*stderr[[2]][1], sigmau2.hat+z*stderr[[2]][2], theta.hat+z*stderr[[2]][3])
  sigmas <- data.frame(Estimate=c(sigmau1.hat, sigmau2.hat, theta.hat), lower=lower.s, upper=upper.s, row.names= c("sigmau1", "sigmau2", "theta"))
  
  beta.u.hat <- BETA.U.saery.MA1(X, ydi, D, md, sigma2edi, sigmau1.hat, sigmau2.hat, theta.hat)
  
  pv <- pvalue(beta.u.hat[[1]], fit)
  lower.b <- beta.u.hat[[1]]-z*stderr[[1]]
  upper.b <- beta.u.hat[[1]]+z*stderr[[1]]
  betas <- data.frame(Estimate=beta.u.hat[[1]], std.err=stderr[[1]], t=beta.u.hat[[1]]/stderr[[1]], 'p-value'=pv, lower=lower.b, upper=upper.b)
  
  return(list("Regression"=betas, "SIGMA"=sigmas, Fsig=fit[[2]], iteration=fit[[3]]))
  
}

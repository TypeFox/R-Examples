eblup.saery.indep <-
function(X, ydi, D, md, sigma2edi, plot){
  
  fit <- REML.saery.indep(X, ydi, D, md, sigma2edi)
  sigmau1.hat <- fit[[1]][1]
  sigmau2.hat <- fit[[1]][2]
  beta.u.hat <- BETA.U.saery.indep(X, ydi, D, md, sigma2edi, sigmau1.hat, sigmau2.hat)
  
  u1d.hat <- beta.u.hat[[2]]
  u1di.hat <- list()
  for(d in 1:D)
    u1di.hat[[d]] <- rep(u1d.hat[d,1],md[d])
  u1di.hat <- unlist(u1di.hat)
  
  mudi.hat <- as.vector(X%*%beta.u.hat[[1]] + u1di.hat + beta.u.hat[[3]])                      
  mse <- mse.saery.indep(X, D, md, sigma2edi, sigmau1.hat, sigmau2.hat, fit[[2]])
  resid <- ydi-mudi.hat
  
  if(plot==TRUE)
    plot.saery(ydi, mudi.hat, mse, sigma2edi)
  
  return(data.frame(Domain=rep(1:D,md), Period=sequence(md), direct=ydi, eblup=mudi.hat, var.direct=sigma2edi, mse.eblup=mse, resid=resid))
  
}

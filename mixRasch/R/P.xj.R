`P.xj` <-
function(delt,th){
  
  ldelt <- length(delt)
  steps <- ldelt - sum(is.na(delt))
  delt <- delt[1:steps]
  delt <- cumsum(delt)    
  stepper <- 1:steps
  all.p <- exp((stepper %o% th) - delt)
  if(! steps == 1) { all.p <- t(t(all.p)/( 1 + colSums(matrix(all.p[stepper,],ncol=length(th))) ))
                 } else all.p <- all.p/(1 + all.p)
  if(steps < ldelt) all.p <- rbind(all.p, array(NA, dim=c(ldelt-steps,ncol(all.p))))
  all.p
}


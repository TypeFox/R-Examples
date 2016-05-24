one.boot.plot <-
function(x, ci, fixed.values, fix.ind, out.ind){

  myboot.sample <- x$functions$boot.sample( x$derived.data$event, 
                                            x$derived.data$trt, 
                                            rho = x$model.fit$cohort.attributes)
  #this makes it work for step function
  if(substr(ci, 1,1)=="h") addind = 0 
  else addind = 1


  rho.b <- myboot.sample[1:7]
  ind   <- myboot.sample[-c(1:7)]

  event.b  <- x$derived.data$event[ind]
  trt.b    <- x$derived.data$trt[ind]

  if(x$model.fit$link == "risks_provided"){
    obsrisk.t0.b <- x$derived.data$fittedrisk.t0[ind]
    obsrisk.t1.b <- x$derived.data$fittedrisk.t1[ind]
    linkinvfun <- NULL
    marker.b <- obsrisk.t0.b - obsrisk.t1.b#
  }else{
    marker.b <- x$derived.data$marker[ind] 
    
    
    coef <- unname(get.coef(event.b, trt.b, marker.b, 
                            x$model.fit$study.design, 
                            rho.b, 
                            link = x$model.fit$link)[,1])
    
    linkinvfun <- binomial(link = x$model.fit$link)$linkinv
    obsrisk.t0.b  <-  get.risk.t0(coef,  marker.b, linkinvfun)
    obsrisk.t1.b  <-  get.risk.t1(coef,  marker.b, linkinvfun)
  }
  
  obsdelta.b <-obsrisk.t0.b - obsrisk.t1.b#

  F.Y <- x$functions$get.F( marker.b,        event.b, trt.b, rho.b)*100#
  F.D <- x$functions$get.F( obsdelta.b, event.b, trt.b, rho.b)*100#
  
  theta.c <- EventRateVec(obsrisk.t0.b, obsrisk.t1.b, F.D, rho.b, event.b, trt.b)

  #all 
  all  <- cbind( F.Y, obsrisk.t0.b, obsrisk.t1.b, F.D, obsdelta.b, theta.c)
  all <- unique(all)
#  browser()
  if(length(fix.ind) > 1){
  myorder <- apply(all[,fix.ind], 2, order)
  
  out <- matrix(0, ncol = length(fixed.values), nrow = length(fix.ind))
  
  for( i in 1:length(fix.ind)){
    #for decreasing pred curves with horizontal bands
    if(is.element(fix.ind[i], c(2,3)) & !(all.equal(order(all[,fix.ind[i]]), order(all[,1]))==TRUE)){ addind = 1} 
    
  tmpind <- sum.I(fixed.values, ">=", all[myorder[,i],fix.ind[i]])+addind
  tmpind[tmpind<=0] <- NA
  tmpall <- all[myorder[,i],out.ind[i]]
  out[i, ] <- tmpall[tmpind]
  
  }

  }else{
  myorder <- order(all[,fix.ind])
  
  out <- numeric( length(fixed.values))
  
  tmpind <- sum.I(fixed.values, ">=", all[myorder,fix.ind]) +addind
  tmpind[tmpind==0] <- NA
  tmpall <- all[myorder,out.ind]
  out <- tmpall[tmpind]
  
  }

  out
}

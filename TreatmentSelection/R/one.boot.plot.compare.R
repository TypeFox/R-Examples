one.boot.plot.compare <-
function(event, trt, marker, marker2, ci,
                                                fixeddeltas.y1, fixeddeltas.y2,
                                                rho = rho, study.design, obp.boot.sample, obp.get.F, fix.ind, out.ind, link, provided_risk = NULL){
  if(substr(ci, 1,1)=="h") addind = 0 
  else addind = 1
  
  myboot.sample <- obp.boot.sample( event, trt, rho)
 
  rho.b <- myboot.sample[1:7]
  ind   <- myboot.sample[-c(1:7)]

  event.b  <- event[ind]
  trt.b  <- trt[ind]
  marker.b  <- marker[ind] 
  marker2.b <- marker2[ind]

  ### marker1


  
  fixed.values <- fixeddeltas.y1
  coef <- unname(get.coef(event.b, trt.b, marker.b, study.design, rho.b, link)[,1])
  if(link == "risks_provided"){
    obsrisk.t0.b <- provided_risk[ind,1]
    obsrisk.t1.b <- provided_risk[ind,2]
  }else{
    linkinvfun <- binomial(link = link)$linkinv
    obsrisk.t0.b  <-  get.risk.t0(coef,  marker.b, linkinvfun)
    obsrisk.t1.b  <-  get.risk.t1(coef,  marker.b, linkinvfun)
    #obsrisk.t0.b  <- c(ifelse(coef[3]>0, 0,1),        get.risk.t0(coef,  marker.b, linkinvfun))
    #obsrisk.t1.b  <- c(ifelse(sum(coef[3:4])>0, 0,1), get.risk.t1(coef,  marker.b, linkinvfun))
    
  }
  
  obsdelta.b <- obsrisk.t0.b- obsrisk.t1.b
  F.D <- obp.get.F( obsdelta.b, event.b, trt.b, rho.b)*100
  
  #obsdelta.b <- c(-1, obsrisk.t0.b[-1] - obsrisk.t1.b[-1])
  #F.D <- c(0, obp.get.F( obsdelta.b[-1], event.b, trt.b, rho.b))*100

  all  <- cbind( F.D, obsdelta.b)
  myorder <- order(all[,fix.ind])
  out <- numeric( length(fixed.values))
  
  tmpind <- sum.I(fixed.values, ">=", all[myorder,fix.ind])
  tmpind[tmpind==0] <- NA
  tmpall <- all[myorder,out.ind]
  out <- tmpall[tmpind]
  
  marker.vals.delta.Y1 <- out

  ##

  ### Y2
fixed.values <- fixeddeltas.y2
coef <- unname(get.coef(event.b, trt.b, marker2.b, study.design, rho.b, link)[,1])
  if(link == "risks_provided"){
    obsrisk.t0.b <- provided_risk[ind,3]
    obsrisk.t1.b <- provided_risk[ind,4]
  }else{
    linkinvfun <- binomial(link = link)$linkinv
    obsrisk.t0.b  <-  get.risk.t0(coef,  marker2.b, linkinvfun)
    obsrisk.t1.b  <-  get.risk.t1(coef,  marker2.b, linkinvfun)
    #obsrisk.t0.b  <- c(ifelse(coef[3]>0, 0,1),        get.risk.t0(coef,  marker2.b, linkinvfun))
    #obsrisk.t1.b  <- c(ifelse(sum(coef[3:4])>0, 0,1), get.risk.t1(coef,  marker2.b, linkinvfun))
    
  }
  obsdelta.b <- obsrisk.t0.b - obsrisk.t1.b
  F.D <- obp.get.F( obsdelta.b, event.b, trt.b, rho.b)*100
  
  #obsdelta.b <- c(-1, obsrisk.t0.b[-1] - obsrisk.t1.b[-1])
  #F.D <- c(0, obp.get.F( obsdelta.b[-1], event.b, trt.b, rho.b))*100

  all  <- cbind( F.D, obsdelta.b)
  myorder <- order(all[,fix.ind])
  out <- numeric( length(fixed.values))
  
  tmpind <- sum.I(fixed.values, ">=", all[myorder,fix.ind]) + addind
  tmpind[tmpind==0] <- NA
  tmpall <- all[myorder,out.ind]
  out <- tmpall[tmpind]

  marker.vals.delta.Y2 <- out

  ##

  return( rbind( marker.vals.delta.Y1, marker.vals.delta.Y2))
}

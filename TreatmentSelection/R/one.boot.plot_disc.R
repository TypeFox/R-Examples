one.boot.plot_disc <-
  function(x){

    myboot.sample <- x$functions$boot.sample( x$derived.data$event, 
                                              x$derived.data$trt, 
                                              x$model.fit$cohort.attributes)
    
    rho.b <- myboot.sample[1:7]
    ind   <- myboot.sample[-c(1:7)]
    
    event.b  <- x$derived.data$event[ind]
    trt.b  <- x$derived.data$trt[ind]
    marker.b  <- x$derived.data$marker[ind] 

    if(x$model.fit$link == "risks_provided") 
    {
      provided_risk <- cbind(x$derived.data$fittedrisk.t0, 
                             x$derived.data$fittedrisk.t1)
    } else provided_risk = NULL
    
    
  
    tmp.trtsel <- trtsel.boot(event.b, 
                              trt.b, 
                              marker.b, 
                              d=x$model.fit$thresh, 
                              rho = rho.b, 
                              x$model.fit$study.design,
                              x$model.fit$link, 
                              x$model.fit$disc.marker.neg, 
                              provided_risk = provided_risk)
    
   unique.fitted.vals <- unique(cbind(marker.b, tmp.trtsel$derived.data[,3:5]))
    
    
    mval <- sort(unique(marker.b))
    out <- with( unique.fitted.vals,  
      c( event.trt0.mkr0 = fittedrisk.t0[marker.b==mval[1]], 
         event.trt1.mkr0 = fittedrisk.t1[marker.b==mval[1]], 
         event.trt0.mkr1 = fittedrisk.t0[marker.b==mval[2]],
         event.trt1.mkr1 = fittedrisk.t1[marker.b==mval[2]], 
         trteff.mkr0 = trt.effect[marker.b==mval[1]], 
         trteff.mkr1 = trt.effect[marker.b==mval[2]]))
    return(out)
  }


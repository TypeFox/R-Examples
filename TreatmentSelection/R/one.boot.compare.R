one.boot.compare <-
function(data1, data2, rho, study.design, obe.boot.sample, obe.get.summary.measures, link, d, disc.marker.neg = NULL){


  myboot.sample <- obe.boot.sample( data1$event, data1$trt, rho)

  rho.b <- myboot.sample[1:7]
  ind   <- myboot.sample[-c(1:7)]

  event.b  <- data1$event[ind]
  trt.b  <- data1$trt[ind]

  if(link == "risks_provided")
  {
     x1.b <- trtsel.boot( event = event.b, trt = trt.b, d = d, study.design = study.design, rho = rho.b, link = link, disc.marker.neg =disc.marker.neg, 
                          provided_risk = cbind(data1$fittedrisk.t0, data1$fittedrisk.t1)[ind,])
     coefs1 <- rep(0,4)
  }else{
     x1.b <- trtsel.boot( event = event.b, trt = trt.b, marker = data1$marker[ind], d = d, study.design = study.design, rho = rho.b, link = link, disc.marker.neg =disc.marker.neg)
     coefs1 <- x1.b$model$coefficients[,1]
  }
  
  if(is.null(data1[["marker.neg"]])){
    x1.b$derived.data$marker.neg <- 1- x1.b$derived.data$marker.neg
    names(x1.b$derived.data)[6] <- "marker.pos"
    
  }
  
  sm1.b <- obe.get.summary.measures(x1.b$derived.data, rho.b)

  if(link == "risks_provided")
  {
    x2.b <- trtsel.boot( event = event.b, trt = trt.b, d = d, study.design = study.design, rho = rho.b, link = link, disc.marker.neg =disc.marker.neg, 
                         provided_risk = cbind(data2$fittedrisk.t0, data2$fittedrisk.t1)[ind,])
    coefs2 <- rep(0, 4)
  }else{
    x2.b <- trtsel.boot( event = event.b, trt = trt.b, marker = data2$marker[ind], d = d, study.design = study.design, rho = rho.b, link = link, disc.marker.neg =disc.marker.neg)
    coefs2 <- x2.b$model$coefficients[,1]
  }
  if(is.null(data2[["marker.neg"]])){
    x2.b$derived.data$marker.neg <- 1- x2.b$derived.data$marker.neg
    names(x2.b$derived.data)[6] <- "marker.pos"
    
  }
  
  sm2.b <- obe.get.summary.measures(x2.b$derived.data, rho.b)

  c(unlist(coefs1), unlist(coefs2), unlist(sm1.b), unlist(sm2.b))



}

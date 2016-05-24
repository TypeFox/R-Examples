one.boot.eval <-
function(data, rho, study.design, obe.boot.sample, obe.get.summary.measures, link, d, disc.marker.neg = NULL, provided_risk = NULL){

  

  sample <- obe.boot.sample( data$event, data$trt, rho)
  rho.b <- sample[1:7]
  ind   <- sample[-c(1:7)]


  x.b <- trtsel.boot( event = data$event[ind], 
                      trt = data$trt[ind], 
                      marker = data$marker[ind], 
                      d = d, 
                      study.design = study.design, 
                      rho = rho.b, 
                      link = link, 
                      disc.marker.neg = disc.marker.neg, 
                      provided_risk = provided_risk[ind,])
 
  if(is.null(data[["marker.neg"]])){
    x.b$derived.data$marker.neg <- 1- x.b$derived.data$marker.neg
    names(x.b$derived.data)[6] <- "marker.pos"
    
  }
  #a3.b <- x.b$model$coefficients[4]
  #a1.b <- x.b$model$coefficients[2]
  coefs <- x.b$model$coefficients[,1]
  #sm = 'summary measures'

  sm.b <- obe.get.summary.measures(x.b$derived.data, rho.b)
  
 
  
#  pdhat  <- sm.b$p.neg
#  neg    <- x.b$derived.data[ind,6] #marker neg or pos
#  marker.b <- data$marker[ind]
# thresh.b <- ifelse(pdhat > 0, max(marker.b[neg==1]), NA)


  #c(a3.b = a3.b, a1.b = a1.b, unlist(sm.b))
  if(is.null(coefs)) coefs <- rep(0, 4)
   c(unlist(coefs), unlist(sm.b))#, thresh.b)

}

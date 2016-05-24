eval.trtsel <-
function(x, bootstraps = 1000, alpha = .05){

  if(!is.trtsel(x)) stop("x must be an object of class 'trtsel' created by using the function 'trtsel' see ?trtsel for more help")
 
  if(alpha<0 | alpha > 1) stop("Error: alpha should be between 0 and 1")
  if(bootstraps ==0 ) print("bootstrap confidence intervals will not be calculated")
  if(bootstraps == 1) warning("Number of bootstraps must be greater than 1, bootstrap confidence intervals will not be computed") 

  data<-x$derived.data
  study.design<-x$model.fit$study.design
  rho<-x$model.fit$cohort.attributes
  boot.sample <- x$functions$boot.sample
  get.summary.measures <- x$functions$get.summary.measures
  marker.bounds <- x$model.fit$marker.bounds
  test.Null.val <- test.Null(x, alpha = alpha)
  
  link <- x$model.fit$link
  
  if(link == "risks_provided") provided_risk <- cbind(x$derived.data$fittedrisk.t0, x$derived.data$fittedrisk.t1)
  else provided_risk = NULL
  

  if(bootstraps > 1){
  #get bootstrap data

  boot.data <- replicate(bootstraps, one.boot.eval(data = data, 
                                                   rho = rho, 
                                                   d = x$model.fit$thresh, 
                                                   study.design = study.design, 
                                                   obe.boot.sample = boot.sample, 
                                                   obe.get.summary.measures = get.summary.measures, 
                                                   link = link, 
                                                   disc.marker.neg = x$model.fit$disc.marker.neg, 
                                                   provided_risk = provided_risk))
 #appease check
  quantile <- NULL
  ## 1. Test the null hypothesis of Theta = 0   


  ## 2. Estimate summary measures
  
  summary.measures <- data.frame(get.summary.measures(data, rho))
  #marker threshold st delta(mthresh) = 0
  if(any(data$marker.neg==0) & any(data$marker.neg==1) &is.null(x$model.fit$disc.marker.neg)& link != "risks_provided"){

    summary.measures$Marker.Thresh <-ifelse( with(data, trt.effect[which.min(marker)]) < 0 , 
                                             max(data$marker[data$marker.neg == 1]), 
                                             min(data$marker[data$marker.neg == 1]))


  }else if(any(data$marker.pos==1) & any(data$marker.pos==0) &is.null(x$model.fit$disc.marker.neg)& link != "risks_provided"){
    
    summary.measures$Marker.Thresh <-ifelse( with(data, trt.effect[which.min(marker)]) < 0 , 
                                             max(data$marker[data$marker.pos == 0]), 
                                             min(data$marker[data$marker.pos == 0]))
    
    
  }else{
  summary.measures$Marker.Thresh <- NA
  }

  conf.intervals <- apply(boot.data[-c(1:4),], 1, quantile, probs = c(alpha/2, 1-alpha/2), type = 1, na.rm = TRUE) 
  row.names(conf.intervals) <- c("lower", "upper")
  #data for more ints, this is to generate ci bounds for our estimates of the thresh hold value 
  # data.mi <- rbind( (-boot.data[2,]/boot.data[4,]), boot.data[15,])
  #more.ints <- apply(data.mi, 1, quantile, probs = c(alpha/2, 1-alpha/2, 1:99/100), na.rm = TRUE) ##########  I am giving off all the quantiles from 1 to 99 here for the simulations, we should take this out for the actual package release. 

  if(!is.null(marker.bounds)){

     #bootstrap distribution of -a1/a3
     a1a3.boot.data <- -boot.data[2,]/boot.data[4,]

     #quantile( a1a3.boot.data , probs = c(alpha/2, 1-alpha/2))
  
     #find the p.value of hat(-a1/a3) falls within the marker bounds
     potential.pvals <- (1:bootstraps)/bootstraps   

     tmp.boot.data <- a1a3.boot.data
     tmp.boot.data <- tmp.boot.data[is.finite(tmp.boot.data)]
      
      if( median(tmp.boot.data, na.rm=TRUE) < marker.bounds[1] |  median(tmp.boot.data, na.rm=TRUE) > marker.bounds[2] ){ 
        a1a3.pval <- 1
      }else if(min(tmp.boot.data) > marker.bounds[1] & max(tmp.boot.data) < marker.bounds[2]){
        a1a3.pval <- 0
      }else{

           reject.all.low <- unname( mapply( cover, 
                                   quantile(tmp.boot.data, potential.pvals/2, , type = 1, na.rm = TRUE),
                                   quantile(tmp.boot.data, 1 - potential.pvals/2, type = 1, na.rm = TRUE), 
                                   rep(marker.bounds[1], bootstraps))  )

           reject.all.high <- unname( mapply( cover, 
                                   quantile(tmp.boot.data, potential.pvals/2, , type = 1, na.rm = TRUE),
                                   quantile(tmp.boot.data, 1 - potential.pvals/2, type = 1, na.rm = TRUE), 
                                   rep(marker.bounds[2], bootstraps))  )

     
          reject.all <- c(!(!reject.all.low*!reject.all.high), FALSE)
          tmp.reject <- which(reject.all==FALSE)[1] 
          a1a3.pval <- potential.pvals[ifelse(tmp.reject==1, 1, tmp.reject - 1)]
      }
      if(a1a3.pval > alpha) test.Null.val$reject <- FALSE
      test.Null.val$a1a3.pval <- a1a3.pval

    }

  result <- list(test.Null          = test.Null.val, 
                 estimates = summary.measures, 
                 conf.intervals   = conf.intervals)#, 
                 #more.ints        = more.ints)
  }else{

  summary.measures <- data.frame(get.summary.measures(data, rho))

  #marker threshold st delta(mthresh) = 0
  if(any(data$marker.neg==0) & any(data$marker.neg==1) &is.null(x$model.fit$disc.marker.neg)& link != "risks_provided"){
    
    summary.measures$Marker.Thresh <-ifelse( with(data, trt.effect[which.min(marker)]) < 0 , 
                                             max(data$marker[data$marker.neg == 1]), 
                                             min(data$marker[data$marker.neg == 1]))
    
    
  }else if(any(data$marker.pos==1) & any(data$marker.pos==0) &is.null(x$model.fit$disc.marker.neg)& link != "risks_provided"){
    
    summary.measures$Marker.Thresh <-ifelse( with(data, trt.effect[which.min(marker)]) < 0 , 
                                             max(data$marker[data$marker.pos == 0]), 
                                             min(data$marker[data$marker.pos == 0]))
    
    
  }else{
    summary.measures$Marker.Thresh <- NA
  }
  
  
  result <- list(test.Null          = test.Null.val, estimates = summary.measures)
  }
  if(!is.null(x$model.fit$disc.marker.neg)) result$discrete.marker = TRUE
  

  class(result) <- "eval.trtsel"

  return(result) 

}

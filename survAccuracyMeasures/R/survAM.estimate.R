survAM.estimate <- function(time, event, marker,
                             data, 
                             predict.time,  
                             marker.cutpoint = 'median', 
                             estimation.method = "IPW", 
                             ci.method = "logit.transformed", 
                             se.method = "bootstrap",
                             bootstraps = 1000, 
                             alpha=0.05){

  # checks
  stopifnot(is.data.frame(data))
  
  time <- eval(substitute(time), data)
  event <- 1*eval(substitute(event), data)
  marker <- eval(substitute(marker), data)

  stopifnot(is.element(estimation.method, c("IPW", "Cox")))
  stopifnot(is.numeric(predict.time))
  if(marker.cutpoint=='median') marker.cutpoint =  median(eval(substitute(marker), data))
  stopifnot(is.numeric(marker.cutpoint))
    stopifnot(is.element(se.method, c("bootstrap", "asymptotic")))
  
  #cant return IPW se estimates 
  if(estimation.method =="IPW" & se.method=="asymptotic") stop("Asymptotic variance calculations are not available for IPW estimates, please use the bootstrap to calculate standard error")
  
  #set some defaults
  measures = c('all')

  if(is.element("all", measures)) measures <- c("AUC","TPR", "FPR", "PPV", "NPV")

  N = nrow(data)

  ## get estimates via getEstimates, also calculate the bootstrap se if necessary
  
  #we handle semiparametric and nonparametric estimates differently
  if(is.element(estimation.method, c("S", "Cox", "Semi-Parametric", "semiparametric"))){
  
    mydata <- prepareDataSP(time, event, marker)
    
    if(se.method == "asymptotic"){
      
      estRawOutput <- getEstimatesSP( data = mydata, 
                                      cutpoint = marker.cutpoint,  
                                      measures = measures,
                                      predict.time = predict.time,
                                      CalVar = TRUE)  
      
    }else if(substr(se.method, 1,4)=="boot"){
      bootstraps = round(bootstraps)
      if(bootstraps <= 1) stop("bootstraps must be larger than 1")
      #estimates
      estRawOutput<-  getEstimatesSP( data = mydata, 
                                                     cutpoint = marker.cutpoint,  
                                                     measures = measures,
                                                     predict.time = predict.time,
                                                     CalVar = FALSE)
      #bootstrap ci's
      bootests <- matrix(ncol = length(estRawOutput$est), nrow = bootstraps)
      for( b in 1:bootstraps){                  
        bootests[b,] <- unlist(getEstimatesSP( data = mydata[sample.int(N, replace = TRUE),], 
                                                        cutpoint = marker.cutpoint,  
                                                        measures = measures,
                                                        predict.time = predict.time,
                                                        CalVar = FALSE)$est) 
      }
      estRawOutput$se <- data.frame(t(apply(bootests, 2, sd, na.rm = TRUE)))
      names(estRawOutput$se) = names(estRawOutput$estimates)
    }
    
  }else if(is.element(estimation.method, c("NP", "IPW"))){
    mydata <- prepareDataNP(time, event, marker)
    
    
    ##no bootstrap
   
      if(se.method == "asymptotic"){
        
        estRawOutput <- getEstimatesNP( data = mydata, 
                                        cutpoint = marker.cutpoint,  
                                        measures = measures,
                                        predict.time = predict.time,
                                        CalVar = TRUE)  
        
      }else if(substr(se.method, 1,4)=="boot"){
        bootstraps = round(bootstraps)
        if(bootstraps <= 1) stop("bootstraps must be larger than 1")
        #estimates
        estRawOutput<-  getEstimatesNP( data = mydata, 
                                                        cutpoint = marker.cutpoint,  
                                                        measures = measures,
                                                        predict.time = predict.time,
                                                        CalVar = FALSE  
                                                       )
        #bootstrap ci's
        bootests <- matrix(ncol = length(estRawOutput$est), nrow = bootstraps)
        for( b in 1:bootstraps){   

          bootests[b,] <-  unlist(getEstimatesNP( data = mydata[sample.int(N, replace = TRUE),], 
                                                          cutpoint = marker.cutpoint,  
                                                          measures = measures,
                                                          predict.time = predict.time,
                                                          CalVar = FALSE)$est) 
        }
 
        estRawOutput$se <- data.frame(t(apply(bootests, 2, sd, na.rm = TRUE)))
        names(estRawOutput$se) = names(estRawOutput$estimates)
      }
    
      
    
    
  }else{
    
    stop("estimation.method not set correctly: it must be one of `Cox` or 'IPW'")
  }
     

     #process the raw estimate data for clean output
     myests <- processRawOutput(estRawOutput, ci.method, alpha)
 
        
  myests$cutpoint = marker.cutpoint; 
  myests$estimation.method = estimation.method; 
  myests$ci.method = ci.method; 
  myests$se.method = se.method;
  myests$predict.time = predict.time; 
  myests$alpha = alpha; 
  
  ## define class
  class(myests) <-  "SurvAM"
  myests

}

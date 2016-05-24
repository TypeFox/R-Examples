print.trtsel <-
function(x, ...){

  cat(paste("Study design:", x$model.fit$study.design, "\n\n"))
  if( substr(x$model.fit$study.design,1,4)  == "nest" ){
  ca = x$model.fit$cohort.attributes
  cat("   Cohort attributes: \n")
  cat(paste("      Cohort sample size =", ca[1], "\n"))
  cat(paste("      Pr(trt == 1) in cohort =", ca[2], "\n"))
  cat(paste("      Pr(event == 1) in cohort =", ca[3], "\n"))
  cat(paste("      fraction of cases sample from cohort =", ca[4], "\n\n"))  

  }
  if( substr(x$model.fit$study.design,1,4)  == "stra" ){
    ca = x$model.fit$cohort.attributes
    cat("   Cohort attributes: \n")
    cat(paste("      Cohort sample size =", ca[1], "\n"))
    cat(paste("      Pr(trt == 0 & event == 0) in cohort =", ca[2], "\n"))
    cat(paste("      Pr(trt == 0 & event == 1) in cohort =", ca[3], "\n"))
    cat(paste("      Pr(trt == 1 & event == 0) in cohort =", ca[4], "\n"))
    cat(paste("      Pr(trt == 1 & event == 1) in cohort =", ca[5], "\n"))
    
    cat(paste("      fraction of cases with trt == 0 sampled from cohort=", ca[6], "\n"))  
    cat(paste("      fraction of cases with trt == 1 sampled from cohort=", ca[7], "\n\n"))
  }
  
  if(ncol(x$derived.data)==6){
    cat("Fitted risks provided. No model fitting implemented.")
    
  }else{
  cat("Model Fit:\n\n")
 
  

    cat(paste(" Link function:", x$model.fit$link, "\n\n"))
  
   
  cat(" Coefficients: \n")
  print(x$model.fit$coefficients)
  }
  cat("\n\n")
  cat("Derived Data: (first ten rows)\n\n")
  print(head(x$derived.data, n=10))

  cat("\n")


}

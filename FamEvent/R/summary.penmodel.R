summary.penmodel <- function(object, ...){

  parms <- round(data.frame(rbind(Estimate=object$parms.est, SE=object$parms.se, RobustSE=object$parms.rse)),4)

  names(parms) <- c("lambda","rho" , "beta.sex","beta.gene")

  pen70 <- round(data.frame(rbind(Estimate=object$pen70.est, SE=object$pen70.se, CI=object$pen70.ci)),4)

  cat("Model Parameters \n")
  print(parms)
  cat("\nPenetrance by age 70 \n")
  print(pen70)
  cat("\n")
  
  out<-list(parameters=parms, pen70=pen70)
  attr(out,"parms") <- parms
  attr(out,"pen70") <- pen70
}

# --------------------------------------------------------
# Description: Function for DTComPair-package
# Author: C. Stock
# Last modified: Feb 15, 2014
# --------------------------------------------------------


# --------------------------------------------------------
# sesp.diff.ci
# --------------------------------------------------------
sesp.diff.ci <- function(tab, ci.method, alpha, cont.corr) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") 
    stop("Table must be of class 'tab.paired'")
  if (missing(alpha)) alpha <- 0.05 
  if (missing(ci.method)) ci.method <- "wald"
  if (ci.method %in% c("wald","agresti-min","bonett-price","tango") == FALSE) 
    stop(paste("The specified ci.method '",ci.method,"' is unknown.", sep="")) 
  if (missing(cont.corr)) cont.corr <- FALSE
  if ((cont.corr==TRUE) & (ci.method!="wald") )
    stop(paste("A continuity correction is only available for ci.method='wald'."))  
  # compute accuracy
  acc <- acc.paired(tab)  
  # sensitivity
  sens.1 <- acc$Test1$sensitivity["est"]; sens.2 <- acc$Test2$sensitivity["est"]
  sens.diff <- (sens.2-sens.1); names(sens.diff) <- NULL  
  # specificity
  spec.1 <- acc$Test1$specificity["est"]; spec.2 <- acc$Test2$specificity["est"]
  spec.diff <- (spec.2-spec.1); names(spec.diff) <- NULL  
  # wald confidence intervals without continuity correction
  if ( (ci.method == "wald") & (cont.corr == FALSE) ) {
    # sensitivity
    b <- tab$diseased[1,2]; c <- tab$diseased[2,1]; n <- tab$diseased[3,3]
    sens.diff.se <- sqrt((b+c) - ((b-c)**2) / n) / n
    sens.diff.cl <- sens.diff + c(-1,1) * qnorm(1-alpha/2) * sens.diff.se
    # specificity
    b <- tab$non.diseased[1,2]; c <- tab$non.diseased[2,1]; n <- tab$non.diseased[3,3]
    spec.diff.se <- sqrt((b+c) - ((b-c)**2) / n) / n
    spec.diff.cl <- spec.diff + c(-1,1) * qnorm(1-alpha/2) * spec.diff.se    
  }
  # wald confidence intervals with continuity correction
  if ( (ci.method == "wald") & (cont.corr == TRUE) ) {
    # sensitivity
    b <- tab$diseased[1,2]; c <- tab$diseased[2,1]; n <- tab$diseased[3,3]
    sens.diff.se <- (sqrt((b+c) - ((b-c)**2) / n) / n) + 1/n
    sens.diff.cl <- sens.diff + c(-1,1) * qnorm(1-alpha/2) * sens.diff.se
    # specificity
    b <- tab$non.diseased[1,2]; c <- tab$non.diseased[2,1]; n <- tab$non.diseased[3,3]
    spec.diff.se <- (sqrt((b+c) - ((b-c)**2) / n) / n) + 1/n
    spec.diff.cl <- spec.diff + c(-1,1) * qnorm(1-alpha/2) * spec.diff.se    
  }
  # agresti-min confidence intervals
  if (ci.method == "agresti-min") {
    k <- 0.5
    # sensitivity    
    b <- tab$diseased[1,2]+k; c <- tab$diseased[2,1]+k; n <- tab$diseased[3,3]+4*k
    sens.diff.se <- (sqrt((b+c) - ((b-c)**2) / n) / n) 
    sens.diff.cl <- sens.diff + c(-1,1) * qnorm(1-alpha/2) * sens.diff.se
    # specificity
    b <- tab$non.diseased[1,2]+k; c <- tab$non.diseased[2,1]+k; n <- tab$non.diseased[3,3]+4*k
    spec.diff.se <- (sqrt((b+c) - ((b-c)**2) / n) / n) 
    spec.diff.cl <- spec.diff + c(-1,1) * qnorm(1-alpha/2) * spec.diff.se    
  }  
  # bonett-price confidence intervals
  if (ci.method == "bonett-price") {
    # sensitivity    
    b <- tab$diseased[1,2]; c <- tab$diseased[2,1]; n <- tab$diseased[3,3]
    p2 <- (b+1) / (n+2); p3 <- (c+1) / (n+2)
    sens.diff.se <- sqrt( (p2 + p3 - (p2-p3)^2) / (n+2) ) 
    sens.diff.cl <- sens.diff + c(-1,1) * qnorm(1-alpha/2) * sens.diff.se
    # specificity
    b <- tab$non.diseased[1,2]; c <- tab$non.diseased[2,1]; n <- tab$non.diseased[3,3]
    p2 <- (b+1) / (n+2); p3 <- (c+1) / (n+2)
    spec.diff.se <- sqrt( (p2 + p3 - (p2-p3)^2) / (n+2) )  
    spec.diff.cl <- spec.diff + c(-1,1) * qnorm(1-alpha/2) * spec.diff.se    
  }  
  # tango confidence intervals
  if (ci.method == "tango") {
    # sensitivity    
    b <- tab$diseased[1,2]; c <- tab$diseased[2,1]; n <- tab$diseased[3,3]
    tango <- scoreci.mp(b, c, n, conf.level=1-alpha)    
    sens.diff.se <- NA    
    sens.diff.cl <- sort(c(tango$conf.int[1], tango$conf.int[2]))
    if ( (tango$conf.int[1] > sens.diff) | (tango$conf.int[2] < sens.diff))
      sens.diff.cl <- sort(-1*sens.diff.cl)
    # specificity
    b <- tab$non.diseased[1,2]; c <- tab$non.diseased[2,1]; n <- tab$non.diseased[3,3]
    tango <- scoreci.mp(b, c, n, conf.level=1-alpha)    
    spec.diff.se <- NA    
    spec.diff.cl <- sort(c(tango$conf.int[1], tango$conf.int[2]))
    if ( (tango$conf.int[1] > spec.diff) | (tango$conf.int[2] < spec.diff))
      spec.diff.cl <- sort(-1*spec.diff.cl)
  } 
  # results
  sensitivity <- c(sens.1, sens.2, sens.diff, sens.diff.se, sens.diff.cl)
  names(sensitivity) <- c("test1","test2","diff","diff.se","diff.lcl","diff.ucl")
  specificity <- c(spec.1, spec.2, spec.diff, spec.diff.se, spec.diff.cl)
  names(specificity) <- c("test1","test2","diff","diff.se","diff.lcl","diff.ucl")
  results <- list(sensitivity, specificity, ci.method, alpha, cont.corr)
  names(results) <- c("sensitivity","specificity","ci.method","alpha","cont.corr")
  return(results)
}

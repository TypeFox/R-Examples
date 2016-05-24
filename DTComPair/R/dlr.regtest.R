# --------------------------------------------------------
# Description: Function for DTComPair-package
# Author: Christian Stock
# Last modified: Jan 29, 2013
# --------------------------------------------------------


# --------------------------------------------------------
# dlr.regtest
# --------------------------------------------------------
dlr.regtest <- function(tab,alpha) {
  # check arguments
  if (missing(tab)) stop("Table is missing.")
  if (class(tab) != "tab.paired") 
    stop("Table must be of class 'tab.paired'")
  if (missing(alpha)) alpha <- 0.05
  acc <- acc.paired(tab)
  # prepare data
  d.wide <- generate.paired(tab)
  d.long <- represent.long(d.wide$d, d.wide$y1, d.wide$y2)
  d.long$y.plus <- 1 - d.long$y
  # pdlr
  pdlr.1 <- acc$Test1$pdlr["est"]; pdlr.2 <- acc$Test2$pdlr["est"]
  names(pdlr.1) <- NULL; names(pdlr.2) <- NULL  
  if (.Platform$OS.type=="windows") sink( tempfile() ) else sink('/dev/null')
  suppressMessages(pdlr <- DLR("~ 1","~  x * y.plus","d", d.long,"id",alpha=alpha)$logDLRModel["x",])
  pdlr[c(1,5,6)] <- exp(pdlr[c(1,5,6)])    
  pdlr <- as.list(c(pdlr.1,pdlr.2,pdlr))
  # ndlr
  ndlr.1 <- acc$Test1$ndlr["est"]; ndlr.2 <- acc$Test2$ndlr["est"]
  names(ndlr.1) <- NULL; names(ndlr.2) <- NULL  
  suppressMessages(ndlr <- DLR("~ 1","~  x * y","d",d.long,"id",alpha=alpha)$logDLRModel["x",])
  ndlr[c(1,5,6)] <- exp(ndlr[c(1,5,6)])
  ndlr <- as.list(c(ndlr.1,ndlr.2,ndlr))
  sink()
  # results
  method <- "DLR regression model (regtest)"
  results <- list(pdlr, ndlr, alpha, method)
  names(results) <- c("pdlr","ndlr","alpha", "method")
  names(results$pdlr) <- c("test1","test2","ratio","se.log","test.statistic","p.value","lcl","ucl")
  names(results$ndlr) <- c("test1","test2","ratio","se.log","test.statistic","p.value","lcl","ucl")
  return(results)  
}



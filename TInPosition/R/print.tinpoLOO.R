print.tinpoLOO <-
function (x,...) {

  res.tinpoLOO <- x
  if (!inherits(res.tinpoLOO, "tinpoLOO")) stop ("no convenient data")
  cat("**TInPosition Leave One Out Cross Validation Data**\n")
  cat("\tFixed Effects (Descriptive) Accuracy: ", res.tinpoLOO$fixed.acc*100,
  "%\n\tRandom Effects (LOO) Accuracy: ", res.tinpoLOO$loo.acc*100, "%\n")  
  cat("*Contains the following objects:\n\n")
  res <- array("", c(6, 2), list(1:6, c("name", "description")))
  
  res[1,] <- c("$loo.assign","The assignment of individuals to groups (after leave one out).")
  res[2,] <- c("$loo.fii","The supplemental factor scores for the items projected after leave one out.")
  res[3,] <- c("$loo.confuse","Confusion matrix for leave one out test.")  
  res[4,] <- c("$loo.acc","Accuracy (% between 0-1) of leave one out estimation.")    
  res[5,] <- c("$fixed.confuse","Confusion matrix for descriptive results (also in $Fixed.Data).")    
  res[6,] <- c("$fixed.acc","Accuracy (% between 0-1) of descriptive results estimation.")      
  
  print(res)

}

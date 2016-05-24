print.mpDOACT.STATIS <-function (x,...) {

  res.mpDOACT.STATIS <- x
  if (!inherits(res.mpDOACT.STATIS, "mpDOACT.STATIS")) stop ("no convenient data")
  cat("**Results for Double STATIS (DO-ACT) **\n")
  cat ("The analysis was performed on ",res.mpDOACT.STATIS$Overview$num.obs,
       "individuals, described by", nrow(res.mpDOACT.STATIS$InnerProduct$C), "tables.\n\n")
  cat("*The results for DO-ACT are available for the following objects:\n\n")
  res <- array("", c(4, 2), list(1:4, c("Name", "Description")))
  
  res[1,] <- c("$Overview","Overview")
  res[2,] <- c("$InnerProduct","Inner Product")
  res[3,] <- c("$Compromise","Compromise")
  res[4,] <- c("$Table","Table")

 
  print(res)
}

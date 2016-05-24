print.mpKPlus1STATIS <-function (x,...) {

  res.mpKPlus1STATIS <- x
  if (!inherits(res.mpKPlus1STATIS, "mpKPlus1STATIS")) stop ("no convenient data")
  cat("**Results for (K+1) STATIS **\n")
  cat ("The analysis was performed on ",res.mpKPlus1STATIS$Overview$num.obs,
       "individuals, described by", nrow(res.mpKPlus1STATIS$InnerProduct$C), "tables.\n\n")
  cat("*The results for (K+1)STATIS are available for the following objects:\n\n")
  res <- array("", c(4, 2), list(1:4, c("Name", "Description")))
  
  res[1,] <- c("$Overview","Overview")
  res[2,] <- c("$InnerProduct","Inner Product")
  res[3,] <- c("$Compromise","Compromise")
  res[4,] <- c("$Table","Table")

 
  print(res)
}

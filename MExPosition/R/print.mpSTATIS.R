print.mpSTATIS <- function (x,...) {

  res.mpSTATIS <- x
  if (!inherits(res.mpSTATIS, "mpSTATIS")) stop ("no convenient data")
  cat("**Results for STATIS **\n")
  cat ("The analysis was performed on ",res.mpSTATIS$Overview$num.obs,
       "individuals, described by", nrow(res.mpSTATIS$InnerProduct$C), "tables.\n\n")
  cat("*The results for STATIS are available for the following objects:\n\n")
  res <- array("", c(4, 2), list(1:4, c("Name", "Description")))
  
  res[1,] <- c("$Overview","Overview")
  res[2,] <- c("$InnerProduct","Inner Product")
  res[3,] <- c("$Compromise","Compromise")
  res[4,] <- c("$Table","Table")

 
  print(res)
}

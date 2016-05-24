print.mpDISTATIS <- function (x,...) {

  res.mpDISTATIS <- x
  if (!inherits(res.mpDISTATIS, "mpDISTATIS")) stop ("no convenient data")
  cat("**Results for DISTATIS **\n")
  cat ("The analysis was performed on ", nrow(res.mpDISTATIS$X),
       "individuals, described by", nrow(res.mpDISTATIS$C), "tables.\n\n")
  cat("*The results for DISTATIS are available for the following objects:\n\n")
  res <- array("", c(4, 2), list(1:4, c("Name", "Description")))
  
  res[1,] <- c("$Overview", "Overview")
  res[2,] <- c("$InnerProduct","Inner Product")
  res[3,] <- c("$Compromise","Compromise")
  res[4,] <- c("$Table","Table")
  
  print(res)
  
}

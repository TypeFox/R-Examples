print.mpANISOSTATIS <- function (x,...) {

  res.mpANISOSTATIS <- x
  if (!inherits(res.mpANISOSTATIS, "mpANISOSTATIS")) stop ("no convenient data")
  cat("**Results for ANISOSTATIS **\n")
  cat ("The analysis was performed on ",res.mpANISOSTATIS$Overview$num.obs,
       "individuals, described by", nrow(res.mpANISOSTATIS$InnerProduct$C), "tables.\n\n")
  cat("*The results for ANISOSTATIS are available for the following objects:\n\n")
  res <- array("", c(5, 2), list(1:5, c("Name", "Description")))
  
  res[1,] <- c("$Overview","Overview")
  res[2,] <- c("$InnerProduct","Inner Product")
  res[3,] <- c("$Compromise","Compromise")
  res[4,] <- c("$Table","Table")
  
  indice = 4
  if(!is.null(res.mpANISOSTATIS$Supp))
  {	res[indice + 1,] <- c("$Supp","Supplementary Projections")
  	indice = indice + 1
  }
 
  print(res[1:indice,])
}

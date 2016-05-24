print.epMDS <-
function (x,...) {

#list(fi=fi,di=di,ci=ci,ri=ri,t=taus,pdq=mds_results)

  res.epMDS <- x
  if (!inherits(res.epMDS, "epMDS")) stop ("no convenient data")
  cat("**Results for Multidimensional Scaling**\n")
  cat ("The analysis was performed on a square matrix of", nrow(res.epMDS$fi),
       "items\n")
  cat("*The results are available in the following objects:\n\n")
  res <- array("", c(10, 2), list(1:10, c("name", "description")))
  
  res[1,] <- c("$fi","Factor scores")
  res[2,] <- c("$di","Squared distances")
  res[3,] <- c("$ci","Contributions")
  res[4,] <- c("$ri", "Cosines")
  res[5,] <- c("$t","Explained Variance")
  res[6,] <- c("$eigs","Eigenvalues")  
  res[7,] <- c("$pdq","SVD data")
  res[8,] <- c("$M","masses")  
  res[9,] <- c("$X","X matrix to decompose")    
  res[10,] <- c("$D","D distance matrix")      
  
  print(res)

}

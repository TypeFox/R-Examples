print.epPCA <-
function (x,...) {

#(list(fi=fi,di=di,ci=ci,ri=ri,fj=fj,cj=cj,rj=rj,dj=dj,t=taus,pdq=pca_results))

  res.epPCA <- x
  if (!inherits(res.epPCA, "epPCA")) stop ("no convenient data")
  cat("**Results for Principal Component Analysis**\n")
  cat ("The analysis was performed on ", nrow(res.epPCA$fi),
       "individuals, described by", nrow(res.epPCA$fj), "variables\n")
  cat("*The results are available in the following objects:\n\n")
  res <- array("", c(16, 2), list(1:16, c("name", "description")))
  
  res[1,] <- c("$fi","Factor scores of the rows")
  res[2,] <- c("$di","Squared distances of the rows")
  res[3,] <- c("$ci","Contributions of the rows")
  res[4,] <- c("$ri", "Cosines of the rows")
  res[5,] <- c("$fj","Factor scores of the columns")
  res[6,] <- c("$dj","square distances of the columns")
  res[7,] <- c("$cj","Contributions for the columns")
  res[8,] <- c("$rj", "Cosines of the columns")
  res[9,] <- c("$t","Explained Variance")
  res[10,] <- c("$eigs","Eigenvalues")  
  res[11,] <- c("$pdq","SVD data")
  res[12,] <- c("$X","X matrix to decompose") 
  res[13,] <- c("$M","Masses - each set to 1")  
  res[14,] <- c("$W","Weights - each set to 1")
  res[15,] <- c("$center","Center of X") 
  res[16,] <- c("$scale","Scale factor of X")   
    
  print(res)

}

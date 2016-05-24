print.tepPLSCA <-
function (x,...) {

#list(fi=fi,di=di,ci=ci,ri=ri,fj=fj,cj=cj,rj=rj,dj=dj,t=taus,M=M,W=W,pdq=pdqFIN)

  res.tepPLSCA <- x
  if (!inherits(res.tepPLSCA, "tepPLSCA")) stop ("no convenient data")
  cat("**Results for Partial Least Squares-Correspondence Analysis**\n")
  cat ("The analysis was performed on ", nrow(res.tepPLSCA$lx),
       "individuals, described by", nrow(res.tepPLSCA$fi), "variables from DATA1 and\n" , nrow(res.tepPLSCA$fj), "variables from DATA2\n")
  cat("*The results are available in the following objects:\n\n")
  res <- array("", c(17, 2), list(1:17, c("name", "description")))
  
  res[1,] <- c("$fi","Factor scores of the groups")
  res[2,] <- c("$di","Squared distances of the groups")
  res[3,] <- c("$ci","Contributions of the groups")
  res[4,] <- c("$ri", "Cosines of the groups")
  res[5,] <- c("$fj","Factor scores of the columns")
  res[6,] <- c("$dj","square distances of the columns")
  res[7,] <- c("$cj","Contributions for the columns")
  res[8,] <- c("$rj", "Cosines of the columns")
  res[9,] <- c("$lx", "Latent variables of X (DATA1)")
  res[10,] <- c("$ly", "Latent variables of Y (DATA2)")    
  res[11,] <- c("$t","Explained Variance")
  res[12,] <- c("$eigs","Eigenvalues")  
  res[13,] <- c("$W1","weights for DATA1")
  res[14,] <- c("$W2","weights for DATA2")
  res[15,] <- c("$pdq","GSVD data")    
  res[16,] <- c("$X","X matrix to decompose")    
  res[17,] <- c("$symmetric","a boolean. TRUE if symmetric scores used for biplot.")  
  
  print(res)

}

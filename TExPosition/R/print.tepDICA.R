print.tepDICA <-
function (x,...) {

#list(fi=fi,di=di,ci=ci,ri=ri,fj=fj,cj=cj,rj=rj,dj=dj,t=taus,M=M,W=W,pdq=pdqFIN)

  res.tepDICA <- x
  if (!inherits(res.tepDICA, "tepDICA")) stop ("no convenient data")
  cat("**Results for Discriminant Correspondence Analysis**\n")
  cat ("The analysis was performed on ", nrow(res.tepDICA$fi),
       "individuals, described by", nrow(res.tepDICA$fj), "variables\n")
  cat("*The results are available in the following objects:\n\n")
  res <- array("", c(23, 2), list(1:23, c("name", "description")))
  
  res[1,] <- c("$fi","Factor scores of the groups")
  res[2,] <- c("$di","Squared distances of the groups")
  res[3,] <- c("$ci","Contributions of the groups")
  res[4,] <- c("$ri", "Cosines of the groups")
  res[5,] <- c("$fj","Factor scores of the columns")
  res[6,] <- c("$dj","square distances of the columns")
  res[7,] <- c("$cj","Contributions for the columns")
  res[8,] <- c("$rj", "Cosines of the columns")
  res[9,] <- c("$lx", "Latent variables of X (DATA)")
  res[10,] <- c("$ly", "Latent variables of Y (DESIGN)")    
  res[11,] <- c("$t","Explained Variance")
  res[12,] <- c("$eigs","Eigenvalues")  
  res[13,] <- c("$M","masses")
  res[14,] <- c("$W","weights")
  res[15,] <- c("$c","center")  
  res[16,] <- c("$pdq","GSVD data")    
  res[17,] <- c("$X","X matrix to decompose")    
  res[18,] <- c("$hellinger","a boolean. TRUE if Hellinger distance was used")        
  res[19,] <- c("$symmetric","a boolean. TRUE if symmetric scores used for biplot.")
  res[20,] <- c("$fii","Factor scores of the individuals")
  res[21,] <- c("$dii","Squared distances of the individuals")
  res[22,] <- c("$rii", "Cosines of the individuals")
  res[23,] <- c("$assign","Information for assignment of individuals to groups")

  
  
  print(res)

}

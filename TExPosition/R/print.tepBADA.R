print.tepBADA <-
function (x,...) {

#list(fi=fi,di=di,ci=ci,ri=ri,fj=fj,cj=cj,rj=rj,dj=dj,t=taus,M=M,W=W,pdq=gpdq_results)

  res.tepBADA <- x
  if (!inherits(res.tepBADA, "tepBADA")) stop ("no convenient data")
  cat("**Results for Barycentric Discriminant Analysis**\n")
  cat ("The analysis was performed on ", nrow(res.tepBADA$fi),
       "individuals, described by", nrow(res.tepBADA$fj), "variables\n")
  cat("*The results are available in the following objects:\n\n")
  res <- array("", c(20, 2), list(1:20, c("name", "description")))
  
  res[1,] <- c("$fi","Factor scores of the rows")
  res[2,] <- c("$di","Squared distances of the rows")
  res[3,] <- c("$ci","Contributions of the rows")
  res[4,] <- c("$ri", "Cosines of the rows")
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
  res[15,] <- c("$pdq","GSVD data")
  res[16,] <- c("$X","X matrix to decompose") 
  res[17,] <- c("$fii","Factor scores of the individuals")
  res[18,] <- c("$dii","Squared distances of the individuals")
  res[19,] <- c("$rii", "Cosines of the individuals")
  res[20,] <- c("$assign","Information for assignment of individuals to groups")  
    
  print(res)

}

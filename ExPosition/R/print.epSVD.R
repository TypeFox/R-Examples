print.epSVD <-
function (x,...) {

#list(fi=fi,di=di,ci=ci,ri=ri,fj=fj,cj=cj,rj=rj,dj=dj,t=taus,M=M,W=W,pdq=pdqFIN)

  res.epSVD <- x
  if (!inherits(res.epSVD, "epSVD")) stop ("no convenient data")
  cat("**Results for SVD**\n")
  cat ("The SVD was performed on ", nrow(res.epSVD$p),
       "individuals, described by", nrow(res.epSVD$q), "variables\n of rank", res.epSVD$rank)
  cat("\n*The results are available in the following objects:\n\n")
  res <- array("", c(7, 2), list(1:7, c("name", "description")))
  
  res[1,] <- c("$p","Left singular vectors.")
  res[2,] <- c("$Dv","Singular values (in a vector).")
  res[3,] <- c("$Dd","Singular values (in diagonal matrix).")
  res[4,] <- c("$q", "Right singular vectors.")   
  res[5,] <- c("$ng", "Number of singular values/vectors.")     
  res[6,] <- c("$rank", "Rank of decomposed matrix. If rank is 1, 0s are padded to singular values and vectors.")       
  res[7,] <- c("$tau", "Explained variance per component.")    
  
  print(res)

}

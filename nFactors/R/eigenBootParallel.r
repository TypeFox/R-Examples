"eigenBootParallel" <-
function(x, quantile=0.95, nboot=30, option="permutation", cor=TRUE, model="components", ...)
 {
 
 if (eigenFrom(x) != "data") stop("Only data from a data.frame must be used as input")

 x          <- data.frame(x)
 res        <- data.frame(matrix(NA, ncol=dim(x)[2], nrow=nboot))
 if (model == "components") { names(res) <- paste("C", 1:dim(x)[2], sep="")
  } else names(res) <- paste("F", 1:dim(x)[2], sep="")

 if (option == "permutation") {
  for (i in 1:nboot) {
   rPerm   <- apply(x,2,sample, replace=TRUE)
   if (cor == TRUE)        corY <- cor(rPerm, ...)
   if (cor == FALSE)       corY <- cov(rPerm, ...)
   if (model == "factors") corY <- corFA(corY, method="ginv")
   res[i,] <- eigen(corY, only.values=TRUE)$values
   }
  }

 if (option == "bootstrap") {
  for (i in 1:nboot) {
   rBoot   <- sample(1:dim(x)[1], dim(x)[1], replace=TRUE)
   if (cor == TRUE)        corY <- cor(x[rBoot,], ...)
   if (cor == FALSE)       corY <- cov(x[rBoot,], ...)
   if (model == "factors") corY <- corFA(corY, method="ginv")
   res[i,] <- eigen(corY, only.values=TRUE)$values
   #if (cor == TRUE)  res[i,] <- eigen(cor(x[rBoot,], ...), only.values=TRUE)$values
   #if (cor == FALSE) res[i,] <- eigen(cov(x[rBoot,], ...), only.values=TRUE)$values
   }
  }
  
 res        <- data.frame(t(moreStats(res, quantile=quantile)))
 return(res)
 }

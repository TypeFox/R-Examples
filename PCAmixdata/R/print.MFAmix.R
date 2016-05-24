print.MFAmix<-function (x, ...) 
{
  res.mfa <- x
  if (!inherits(res.mfa, "MFAmix")) 
    stop("non convenient data")
  cat("**Results of the Multiple Factor Analysis for mixed data (MFAmix)**\n")
  cat("The analysis was performed on", nrow(res.mfa$global.pca$rec$X.quanti), 
      "individuals, described by", ncol(res.mfa$global.pca$rec$X), "variables\n")
  cat("*Results are available in the following objects :\n\n")
  res <- array("", c(12, 2), list(1:12, c("name", "description")))
  res[1, ] <- c("$eig", "eigenvalues")
  res[2, ] <- c("$eig.separate", "eigenvalues of the separate analyses")
  res[3, ] <- c("$separate.analyses", "separate analyses for each group of variables")
  res[4, ] <- c("$groups", "results for all the groups")
  res[5, ] <- c("$partial.axes", "results for the partial axes")
  res[6, ] <- c("$ind", "results for the individuals")
  res[7, ] <- c("$ind.partial", "results for the partial individuals")
  res[8, ] <- c("$quanti", "results for the quantitative variables")
  res[9, ] <- c("$levels", "results for the levels of the qualitative variables")
  res[10, ] <- c("$quali", "results for the qualitative variables")
  res[11,] <- c("$sqload", "squared loadings")
  res[12, ] <- c("$global.pca", "results for the global PCA")
  
  print(res)
  
}

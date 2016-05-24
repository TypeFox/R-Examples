print.doact.statis.compromise <- function(x,...){
	
  res.doact.statis.compromise <- x
  if (!inherits(res.doact.statis.compromise, "doact.statis.compromise")) stop ("no convenient data")
  cat("**Results from the analysis of the Compromise via DO-ACT **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(14, 2), list(1:14, c("Name", "Description")))
	
	res[1,] <- c("$compromiseMatrix.1", "Compromise Matrix for Dataset 1")
	res[2,] <- c("$compromise.ci.1", "Contributions for the rows (Dataset 1)")
	res[3,] <- c("$compromise.cj.1","Contribution for the colums (Dataset 1) ")
	res[4,] <- c("$compromise.tau.1","Percent of Variance explained (Dataset 1)")
	res[5,] <- c("$compromise.eigs.1","Eigen Values (Dataset 1)")
	res[6,] <- c("$compromise.eigs.vectors.1","Eigen Vectors (Dataset 1)")
	res[7,] <- c("$compromise.fi.1","Factor Scores (Dataset 1)")
	res[8,] <- c("$compromiseMatrix.2", "Compromise Matrix for Dataset 2")
	res[9,] <- c("$compromise.ci.2", "Contributions for the rows (Dataset 2)")
	res[10,] <- c("$compromise.cj.2","Contribution for the columns (Dataset 2)")
	res[11,] <- c("$compromise.tau.2","Percent of Variance explained (Dataset 2)")
	res[12,] <- c("$compromise.eigs.2","Eigen Values (Dataset 2)")
	res[13,] <- c("$compromise.eigs.vectors.2","Eigen Vectors (Dataset 2)")
	res[14,] <- c("$compromise.fi.2","Factor Scores (Dataset 2)")

	print(res)
}
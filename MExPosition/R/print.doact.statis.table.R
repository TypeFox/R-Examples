print.doact.statis.table <- function(x,...){
	
  res.doact.statis.table <- x
  if (!inherits(res.doact.statis.table, "doact.statis.table")) stop ("no convenient data")
  cat("**Results from the analysis of the Tables via DO-ACT **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(19, 2), list(1:19, c("Name", "Description")))
	
	res[1,] <- c("$m", "Masses")
	res[2,] <- c("$table.eigs.1", "Eigen Values (Dataset 1)")
	res[3,] <- c("$table.eigs.vector.1","Eigen Vectors (Dataset 1) ")
	res[4,] <- c("$table.fi.1","Factor Scores (Dataset 1)")
	res[5,] <- c("$table.ci.1","Contribution of the rows (Dataset 1)")
	res[6,] <- c("$table.cj.1","Contribution of the columns (Dataset 1")
	res[7,] <- c("$table.tau.1","Percent of Variance Explained (Dataset 1)")
	res[8,] <- c("$table.partial.fi.1", "Partial Factor Scores (Dataset 1)")
	res[9,] <- c("$table.partial.fi.array.1", "Array of Parital Factor Scores (Dataset 1)")
	res[10,] <- c("$table.loadings.1","Loadings (Dataset 1)")
	res[11,] <- c("$table.eigs.2", "Eigen Values (Dataset 2)")
	res[12,] <- c("$table.eigs.vector.2","Eigen Vectors (Dataset 2) ")
	res[13,] <- c("$table.fi.2","Factor Scores (Dataset 2)")
	res[14,] <- c("$table.ci.2","Contribution of the rows (Dataset 2)")
	res[15,] <- c("$table.cj.2","Contribution of the columns (Dataset 2")
	res[16,] <- c("$table.tau.2","Percent of Variance Explained (Dataset 2)")
	res[17,] <- c("$table.partial.fi.2", "Partial Factor Scores (Dataset 2)")
	res[18,] <- c("$table.partial.fi.array.2", "Array of Parital Factor Scores (Dataset 2)")
	res[19,] <- c("$table.loadings.2","Loadings (Dataset 2)")

	print(res)
}
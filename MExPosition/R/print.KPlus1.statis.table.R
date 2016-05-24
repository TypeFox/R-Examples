print.KPlus1.statis.table <- function(x,...){
	
  res.KPlus1.statis.table <- x
  if (!inherits(res.KPlus1.statis.table, "KPlus1.statis.table")) stop ("no convenient data")
  cat("**Results from the analysis of the Tables via (K+1) STATIS **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(19, 2), list(1:19, c("Name", "Description")))
	
	res[1,] <- c("$m", "Masses")
	res[2,] <- c("$table.eigs.", "Eigen Values")
	res[3,] <- c("$table.eigs.vector","Eigen Vectors")
	res[4,] <- c("$table.fi","Factor Scores")
	res[5,] <- c("$table.ci","Contribution of the rows")
	res[6,] <- c("$table.cj","Contribution of the columns")
	res[7,] <- c("$table.t","Percent of Variance Explained")
	res[8,] <- c("$partial.fi", "Partial Factor Scores)")
	res[9,] <- c("$partial.fi.array.1", "Array of Parital Factor Scores")
	res[10,] <- c("$Q","Loadings")
	
	print(res)
}
print.covstatis.table <- function(x,...){
	
  res.covstatis.table <- x
  if (!inherits(res.covstatis.table, "covstatis.table")) stop ("no convenient data")
  cat("**Results from the analysis of the Tables via COVSTATIS **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(10, 2), list(1:10, c("Name", "Description")))
	
	res[1,] <- c("$m", "Masses")
	res[2,] <- c("$eigs", "Eigen Values")
	res[3,] <- c("$eigs.vector","Eigen Vectors")
	res[4,] <- c("$fi","Factor Scores")
	res[5,] <- c("$ci","Contribution of the rows")
	res[6,] <- c("$cj","Contribution of the columns")
	res[7,] <- c("$t","Percent of Variance Explained")
	res[8,] <- c("$partial.fi", "Partial Factor Scores")
	res[9,] <- c("$partial.fi.array", "Array of Parital Factor Scores")
	res[10,] <- c("$Q","Loadings")

	print(res)
}
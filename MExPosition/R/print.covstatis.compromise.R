print.covstatis.compromise <- function(x,...){
	
  res.covstatis.compromise <- x
  if (!inherits(res.covstatis.compromise, "covstatis.compromise")) stop ("no convenient data")
  cat("**Results from the analysis of the Compromise via COVSTATIS **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(7, 2), list(1:7, c("Name", "Description")))
	
	res[1,] <- c("$compromise", "Compromise Matrix")
	res[2,] <- c("$compromise.ci", "Contributions for the rows")
	res[3,] <- c("$compromise.cj","Contribution for the colums")
	res[4,] <- c("$compromise.tau","Percent of Variance explained")
	res[5,] <- c("$compromise.eigs","Eigen Values")
	res[6,] <- c("$compromise.eigs.vectors","Eigen Vectors")
	res[7,] <- c("$compromise.fi","Factor Scores")

	print(res)
}
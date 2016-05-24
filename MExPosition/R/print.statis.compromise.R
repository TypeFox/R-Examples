print.statis.compromise <- function(x,...){
	
	res.statis.compromise <- x
	if (!inherits(res.statis.compromise, "statis.compromise")) stop ("no convenient data")
  cat("**Results from the analysis of the Compromise via STATIS **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(7, 2), list(1:7, c("Name", "Description")))
	
	res[1,] <- c("$compromise", "Compromise Matrix")
	res[2,] <- c("$compromise.ci","Contributions of the rows of Compromise")
	res[3,] <- c("$compromise.cj","Contribution of the columns of the Compromise")
	res[4,] <- c("$compromise.eigs.vector","Eigen Vectors of Compromise")
	res[5,] <- c("$compromise.eigs","Eigen Values of Compromise")
	res[6,] <- c("$compromise.fi","Factor Scores of Compromise")
	res[7,] <- c("$compromise.t","% Variance Explained of Compromise")

	print(res)
}
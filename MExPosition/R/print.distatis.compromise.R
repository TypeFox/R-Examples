print.distatis.compromise <- function(x,...){
	
	res.distatis.compromise <- x
	if(!inherits(res.distatis.compromise,"distatis.compromise"))stop('no convenient data')
	cat("**Results for the Compromise via DISTATIS**\n\n")
	cat("The results are available for the following objects:\n\n")
	res <- array("",c(7,2),list(1:7,c("Name","Description")))
	
	res[1,] <- c("$compromise","Compromise Matrix")
	res[2,] <- c("$compromise.eigs","Eigen Values of Compromise")
  	res[3,] <- c("$compromise.eigs.vector","Eigen Vectors of Compromise")
  	res[4,] <- c("$compromise.fi","Factor Scores of Compromise")
  	res[5,] <- c("$compromise.tau", "Percent of Variance Explained")
  	res[6,] <- c("$compromise.ci","Contribution of the rows of Compromise")
  	res[7,] <- c("$compromise.cj", "Contribution of the columns of Compromise")

	print(res)
}
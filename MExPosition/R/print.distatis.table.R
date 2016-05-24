print.distatis.table <- function(x,...){
	
	res.distatis.table <- x
	if(!inherits(res.distatis.table,"distatis.table"))stop("no convenient data")
	cat("**Results for the Tables via DISTATIS**\n\n")
	cat("The results are available for the following objects: \n\n")
	res <- array("",c(10,2),list(1:10,c("Name","Description")))
	
	res[1,] <- c("$m","Masses")
	res[2,] <- c("$eigs","Eigen Values for the Tables")
	res[3,] <- c("$eigs.vector","Eigen Vectors for the Tables")
	res[4,] <- c("$fi","Factor Scores for the Tables")
	res[5,] <- c("$t","Percent of Variance Explained")
	res[6,] <- c("$partial.fi","Partial Factor Scores for the Tables")
	res[7,] <- c("$partial.fi.array","Array of Partial Factor Scores for the Tables")
	res[8,] <- c("$ci","Contribution of the rows of the Tables")
	res[9,] <- c("$cj","Contribution of the columns of the Tables")
	res[10,] <- c("$Q","Loadings of the Tables")

	print(res)	
}
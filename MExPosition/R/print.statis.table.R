print.statis.table <- function(x,...){
	
	res.statis.table <- x
	if(!inherits(res.statis.table,"statis.table"))stop("no convenient data")
	cat("**Results for the Tables via STATIS**\n\n")
	cat("The results are available for the following objects: \n\n")
	res <- array("",c(9,2),list(1:9,c("Name","Description")))
	
	res[1,] <- c("$m","Masses")
	res[2,] <- c("$eigs","Eigen Values for the Tables")
	res[3,] <- c("$t","Inertia for the Tables (tau)")
	res[4,] <- c("$eigs.vector","Eigen Vectors for the Tables")
	res[5,] <- c("$fi","Factor Scores for the Tables")
	res[6,] <- c("$partial.fi","Partial Factor Scores for the Tables")
	res[7,] <- c("$partial.fi.array","Array of Partial Factor Scores for the Tables")
	res[8,] <- c("$cj","Contributions of the Columns for the Tables")
	res[9,] <- c("$ci","Contributions of the rows for the Tables")

	print(res)	
}
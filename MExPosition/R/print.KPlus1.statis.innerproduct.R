print.KPlus1.statis.innerproduct <- function(x,...){
	
	res.KPlus1.statis.innerproduct <- x
	if (!inherits(res.KPlus1.statis.innerproduct, "KPlus1.statis.innerproduct")) stop ("no convenient data")
    cat("**Results from the analysis of the between table structure via (K+1) STATIS **\n")
    cat("*The results are available for the following objects:\n\n")
    res <- array("", c(11, 2), list(1:11, c("Name", "Description")))
	
	res[1,] <- c("$S", "i X i X j array of S Matrices")
	res[2,] <- c("S.star", "i X i X h array of S* Matrices")
	res[3,] <- c("$C","k X k, C* Matrix")
	res[4,] <- c("$rvMatrix","RV* Matrix")
	res[5,] <- c("$eigs.vector","Eigen Vectors of C*")
	res[6,] <- c("$eigs","Eigen Values of C*")
	res[7,] <- c("$fi","Factor Scores of C*")
	res[8,] <- c("$ci","Contributions of the Rows of C*")
	res[9,] <- c("$cj","Contributions of the Columns of C*")
	res[10,] <- c("$t","% Variance Explained of C*")
	res[11,] <- c("$alphaWeights","Alpha Weights")
	
	print(res)
}
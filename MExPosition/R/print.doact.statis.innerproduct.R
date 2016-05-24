print.doact.statis.innerproduct <- function(x,...){
	
	res.doact.statis.innerproduct <- x
	if (!inherits(res.doact.statis.innerproduct, "doact.statis.innerproduct")) stop ("no convenient data")
  cat("**Results from the analysis of the between table structure via DO-ACT **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(12, 2), list(1:12, c("Name", "Description")))
	
	res[1,] <- c("$S.1", "i X i X j array of S Matrices for dataset 1")
	res[2,] <- c("$S.2", "i X i X j array of S Matrices for dataset 2")
	res[3,] <- c("$C","k X k, C Matrix")
	res[4,] <- c("$RVMatrix","RV Matrix")
	res[5,] <- c("$eigs.vector","Eigen Vectors of C")
	res[6,] <- c("$eigs","Eigen Values of C")
	res[7,] <- c("$fi","Factor Scores of C")
	res[8,] <- c("$ci","Contributions of the Rows of C")
	res[9,] <- c("$cj","Contributions of the Columns of C")
	res[10,] <- c("$t","% Variance Explained of C (tau)")
	res[11,] <- c("$alphaWeights","Alpha Weights")
	res[12,] <- c("$betaWeights","Beta Weights")
	
	print(res)
}
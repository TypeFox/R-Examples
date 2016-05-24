print.distatis.innerproduct <- function(x,...){
	
	res.distatis.innerproduct <- x
	if(!inherits(res.distatis.innerproduct,"distatis.innerproduct"))stop('no convenient data')
	cat("**Results for the Inner Product via DISTATIS**\n\n")
	cat("The results are available for the following objects:\n\n")
	res <- array("",c(10,2),list(1:10,c("Name","Description")))
	
	res[1,] <- c("$S","i X i X j array of S Matrices")
	res[2,] <- c("$norm.S","Normaization of S Matrices")
  	res[3,] <- c("$C","k X K, C Matrix")
  	res[4,] <- c("$eigs.vector","Eigen Vectors of C")
  	res[5,] <- c("$eigs", "Eigen Values of C")
  	res[6,] <- c("$t","Percent variance explained (tau)")
  	res[7,] <- c("$fi", "Factor Scores of C")
  	res[8,] <- c("$ci","Contributions of the Rows of C")
	res[9,] <- c("$cj","Contributions of the Columns of C")
	res[10,] <- c("$alphaWeights","Alpha Weights")


	print(res)
}
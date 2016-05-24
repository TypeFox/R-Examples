summary.YPmodel.overall <-
function(object=c(),...)
{

	Data <- object

	n <- Data$length	
	Delta <- Data$Delta
	Z <- Data$Z
	X <- Data$X
	CensorRata <- 1 - (sum(Delta)/n)

	GroupNum <- matrix(c(sum(Z),n-sum(Z)), nrow=1, ncol=2)
	colnames(GroupNum) <- c("Control Group","Treatment Group")
	rownames(GroupNum) <- c("Numbers")


	cat("\n-------------------------------------------------------------------------------------------------------------  \n")
	cat("Overall of short-term and long-term hazard ration model  \n")
	cat("\n Censoring rate:\n")
	cat(paste(round(CensorRata,digits = 4)),"\n")
	cat("\n Total Number:\n")
	cat(paste(n),"\n")
	cat("\n Group Number:\n")
	printCoefmat(GroupNum)

}

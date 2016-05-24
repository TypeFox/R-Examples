print.distatis.overview <- function(x,...){
	
	res.distatis.overview <- x
	if(!inherits(res.distatis.overview,"distatis.overview")) stop('no convenient data')
	cat("**Overview of Results for DISTATIS**\n\n")
	cat("The results are available for the following objects:\n\n")
	res <- array("",c(5,2),list(1:5,c("Name","Description")))
	
	res[1,] <- c("$data", "Data Matrix")
  	res[2,] <- c("$sorting","Is this a sorting task?")
	res[3,] <- c("$normalization","Type of normalization selected")
	res[4,] <- c("$table","Table matrix which identifies the tables")
	res[5,] <- c("$num.groups","Num of groups")	
	
	print(res)
}
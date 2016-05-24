print.KPlus1.statis.overview <- function(x,...){
	
	res.KPlus1.statis.overview <- x
	if (!inherits(res.KPlus1.statis.overview, "KPlus1.statis.overview")) stop ("no convenient data")
    cat("**Oveview of the Analysis via (K+1) STATIS **\n")
    cat("*The results are available for the following objects:\n\n")
    res <- array("", c(7, 2), list(1:7, c("Name", "Description")))
	
	res[1,] <- c("$data", "Preprocessed Data")
	res[2,] <- c("$plus1data", "Preprocessed External Table")
	res[3,] <- c("$column.design","Design for Data")
	res[4,] <- c("$row.preprocess","Type of row preprocessing selected")
	res[5,] <- c("$column.preprocess","Type of column preprocessing selected")
	res[6,] <- c("$table.preprocess","Type of table preprocessing selected")
	res[7,] <- c("$num.groups","Number of groups")

	print(res)
}
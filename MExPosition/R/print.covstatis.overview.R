print.covstatis.overview <- function(x,...){
	
	res.covstatis.overview <- x
	if (!inherits(res.covstatis.overview, "covstatis.overview")) stop ("no convenient data")
  cat("**Oveview of the Analysis via COVSTATIS **\n")
  cat("*The results are available for the following objects:\n\n")
  res <- array("", c(4, 2), list(1:4, c("Name", "Description")))
	
	res[1,] <- c("$data", "Preprocessed Data")
	res[2,] <- c("$normalization", "Selected type of Normalization")
	res[3,] <- c("$table","Design for Data")
	res[4,] <- c("$num.groups","Number of groups")

	print(res)
}
`checkTableList` <-
function(listTables){
	isNotMat <- sapply(listTables, function(x) !is.matrix(x))
	if(any(isNotMat))
		stop("All entries in listTables must be matrices.")
	counts <- sapply(listTables, function(x) any(x<0))
	if(any(counts))
		stop("All values in all matrices must be non-negative integers.")
	d <- dim(listTables[[1]])
	rn <- rownames(listTables[[1]])
	cn <- colnames(listTables[[1]])
	for(i in 2:length(listTables)){
		if(any(dim(listTables[[i]])!=d))
			stop("The dimensions differ between the matrices.")
		if(any(rownames(listTables[[i]])!=rn))
			stop("The row names differ between the matrices.")
		if(any(colnames(listTables[[i]])!=cn))
			stop("The column names differ between the matrices.")
	}
	cn
}


nominalCheck <- function(DATA){
	data.dim <- dim(DATA)
	colSums.data <- colSums(DATA)
	orig.cols <- sum(cumsum(colSums.data) %% data.dim[1]==0)	
	
	if( (!max(cumsum(colSums.data) / data.dim[1]) == orig.cols) || (!sum(cumsum(colSums.data) %% data.dim[1] == 0) == orig.cols)){
		#ESCAPE!
		stop("Data is not nominal (disjunctive).")
	}else{
		return(DATA)
	}
	
}
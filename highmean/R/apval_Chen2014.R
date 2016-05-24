apval_Chen2014 <- function(sam1, sam2, eq.cov = TRUE){
	if(eq.cov){
		out <- apval_Chen2014_samecov(sam1, sam2)
	}else{
		out <- apval_Chen2014_diffcov(sam1, sam2)
	}
	return(out)
}
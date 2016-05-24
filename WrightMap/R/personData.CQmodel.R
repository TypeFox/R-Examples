personData.CQmodel <- function(thetas, ...) {
	
	model <- thetas
	
	if(!is.null(model$p.est)) {
	p.est <- model$p.est
	columns.at <- grep("^est", names(p.est), perl = TRUE)
	thetas <- p.est[columns.at]
	}
	else
		thetas <- 0
		
	names(thetas) <- model$dimensions
	return(thetas)
	
}
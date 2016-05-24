trendMatrix.update <- function(model, Xnew) { 

		colnames(Xnew) <- colnames(model@X)
		rownames(Xnew) <- nrow(model@X) + 1:nrow(Xnew)
		Fnew <- model.matrix(model@trend.formula, data=Xnew)
		Fupdated <- rbind(model@F, Fnew)
		return(Fupdated)	
}

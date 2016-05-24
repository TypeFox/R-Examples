createDefaultDesign <-
function(DATA){
		DESIGN <- matrix(1,dim(DATA)[1],1)
		rownames(DESIGN) <- rownames(DATA)	
		return(DESIGN)
}

require(RColorBrewer)

colour_var <- function(var, cols=NULL){
	var <- as.factor(var)
	n <- max(3,length(levels(var)))
	if (is.null(cols))
	 # cols <- rainbow(n,alpha=0.5,s=.5)
	   cols <- brewer.pal(n,"Set2")
     return(cols[as.numeric(var)])
	
	}
as.mlcm.df <- function(d, ...){
	if (!inherits(d, "data.frame")) 
		stop("Object must have class data.frame!\n")
	if((length(d) < 5) || !(length(d) %% 2))
		stop("Data frame must be at least 5 columns and have an odd number of columns!\n")
		class(d) <- c("mlcm.df", "data.frame")
		d
}
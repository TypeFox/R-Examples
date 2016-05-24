indToClass <- function(x){
	
	counts <- x[,2]
	u <- sort(unique(counts))
	
	class.table <- data.frame()
	for (i in 1:length(u)){
		tmp <- counts[counts == u[i]]
		tmp2 <- c(u[i], length(tmp))
		class.table <- rbind(class.table, tmp2)
	}
	colnames(class.table) <- c("capture.class", "No.Ind")
	return(class.table)

}


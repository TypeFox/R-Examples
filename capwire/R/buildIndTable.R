buildIndTable <- function(counts){
	
	ind.table <- cbind(c(1:length(counts)), counts)
	colnames(ind.table) <- c("ID", "No.Obs")
	return(ind.table)
	
}
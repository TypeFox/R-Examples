build.subtable <-
function(table, column1, column2, separator='///') {
	resultTable = c()
	indexNonEmpty = grep('.+', table[[column2]])
	for(i in indexNonEmpty) {
		values = unlist(strsplit(as.character(table[i,column2]), separator))
		resultTable = rbind(resultTable, cbind(table[i, column1], values))
	}
	resultTable = data.frame(resultTable, stringsAsFactors=F)
	colnames(resultTable) = c(column1, column2)
	return(resultTable)
}

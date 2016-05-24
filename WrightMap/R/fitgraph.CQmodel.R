fitgraph.CQmodel <-
function(fitEst,table = NULL,fit.type = "W",itemLabels =NULL,...) {
	#message(table)
	if(is.null(table))
		table.num <- 1
	else
		table.num <- match(table,names(fitEst$RMP))
	table <- fitEst$RMP[[table.num]]
	estlabel <- paste(fit.type,"fit",sep=".")
	llabel <- paste(fit.type,"Low",sep=".")
	ulabel <- paste(fit.type,"High",sep=".")
	
	if(is.null(itemLabels)) {
		#message(table.num)
		colNames <- unlist(strsplit(names(fitEst$RMP)[table.num],"*",fixed = TRUE))
		#message(colNames)
		itemLabels <- do.call(paste,table[colNames])
	}
	fitgraph(unlist(table[estlabel]),unlist(table[llabel]),unlist(table[ulabel]),itemLabels)
}

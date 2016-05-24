expressionPlot <- function(bisepData=data, gene1, gene2)
{
	if(missing(bisepData)) stop("Need to input expression data matrix")
	if(missing(gene1)) stop("Need to specify gene 1")
	if(missing(gene2)) stop("Need to specify gene 2")

	# Extract objects from input + stop if object incorrect

	if("BISEP" %in% names(bisepData) && "BI" %in% names(bisepData) && "DATA" %in% names(bisepData))
	{
		biIndex <- bisepData$BI
		big.model <- bisepData$BISEP
		data2 <- bisepData$DATA
	}
	else
	{
		stop("Input object isn't from BISEP function")
	}
	# Do some gene formatting
	gene1 <- toupper(gene1)
	gene2 <- toupper(gene2)
	if(length(which(rownames(data2) %in% gene1)) == 0) stop("Gene 1 not recognised")
	if(length(which(rownames(data2) %in% gene2)) == 0) stop("Gene 2 not recognised")
	
	plot(as.numeric(data2[gene1,]), as.numeric(data2[gene2,]), pch=16, main=paste(gene1, "vs.", gene2, "Log2 Gene Expression plot", sep=" "), xlab=gene1, ylab=gene2)
	abline(v=big.model[gene1,1], col="red")
	abline(h=big.model[gene2,1], col="red")
}


waterfallPlot <- function(bisepData=data, mutData=mutData, expressionGene, mutationGene)
{
	if(missing(bisepData)) stop("Need to input BISEP object")
	if(missing(mutData)) stop("Need to input mutation data matrix")
	if(missing(expressionGene)) stop("Need to specify gene 1")
	if(missing(mutationGene)) stop("Need to specify gene 2")

	# Extract objects from input + stop if object incorrect

	if("BISEP" %in% names(bisepData) && "BI" %in% names(bisepData) && "DATA" %in% names(bisepData))
	{
		biIndex <- bisepData$BI
		big.model <- bisepData$BISEP
		data <- bisepData$DATA
	}
	else
	{
		stop("Input object isn't from BISEP function")
	}
	if("WT" %in% unique(mutData[,1]) || "MUT" %in% unique(mutData[,1]))
	{
		
		# Do some gene formatting
		gene1 <- toupper(expressionGene)
		gene2 <- toupper(mutationGene)
		if(length(which(rownames(data) %in% gene1)) == 0) stop("Gene 1 not recognised")
		if(length(which(rownames(mutData) %in% gene2)) == 0) stop("Gene 2 not recognised")
	
		# Sort out your matrices
		mutData2 <- mutData[,which(colnames(mutData) %in% colnames(data))]
		data2 <- data[,which(colnames(data) %in% colnames(mutData2))]
		data2 <- data2[,colnames(mutData2)]
		d1 <- as.data.frame(cbind(colnames(data2), as.numeric(data2[gene1,]), as.character(mutData2[gene2,])), stringsAsFactors=FALSE)
		colnames(d1) <- c("sample_name", gene1, gene2)
		d1[,2] <- as.numeric(d1[,2])
		d1 <- d1[order(d1[,2]),]
		d1$Colour <- "orange"
		d1$Colour[which(d1[,3] == "MUT")] <- "black"
		y <- 2+3*d1[,2] + rnorm(length(d1[,2]))
		d.x <- density(d1[,2])
		d.y <- density(y)
		layout( matrix( c(1,1,2,2,2,2,2,2,2,2),ncol=5) )
		plot(d.y$y, d.y$x, ylim=range(y), xlim=rev(range(d.y$y)), type='l')
		barplot((d1[,2] - big.model[gene1,1]), col=d1[,4], border=NA, main=paste(gene1, " Expression" ," coloured by mutation status of ", gene2, sep=""), space=F, cex.main=2)
		legend("bottomright", c("WT", "MUT"), col=c("orange", "black"), pch=15)
	}
	else
	{
		stop("Input object isn't of required type - containing only 'MUT' or 'WT' values")
	}
}


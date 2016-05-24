plotSuper <- function(x, data, BiclustSet){
	superBiclusterSamples <- c()
	superBiclusterGenes <- c()
	
	iGroupColumns <- BiclustSet@ColumnMembership[x, ]
	iGroupRows <- BiclustSet@GenesMembership[, x] 
	superBiclusterSamples <-  as.logical(colSums(iGroupColumns) == nrow(iGroupColumns))
	superBiclusterGenes  <- as.logical(rowSums(iGroupRows) == ncol(iGroupRows))
	
	plotProfilesWithinBicluster(x = data[superBiclusterGenes, superBiclusterSamples], 
			sampleNames = colnames(data)[superBiclusterSamples], 
			main = "gene profiles")
	return(new("Biclust", Number=1, RowxNumber=matrix(superBiclusterGenes, ncol=1), NumberxCol=matrix(superBiclusterSamples, nrow=1),
					info = list(Call="super bicluster",size = length(x)),
					Parameters = list(Call=paste("Obtained from ",length(x),"biclusters" )) ))
}

plotSuperAll <- function(x, data, BiclustSet){
	superBiclusterSamples <- c()
	superBiclusterGenes <- c()
	
	iGroupColumns <- BiclustSet@ColumnMembership[x, ]
	iGroupRows <- BiclustSet@GenesMembership[, x] 
	superBiclusterSamples <-  as.logical(colSums(iGroupColumns) == nrow(iGroupColumns))
	superBiclusterGenes  <- as.logical(rowSums(iGroupRows) == ncol(iGroupRows))

	plotProfilesAcrossAllSamples(x = data,superBiclusterGenes, superBiclusterSamples)
}

plotProfilesAcrossAllSamples <- function(x, coreBiclusterGenes, coreBiclusterSamples){
	namesS <- colnames(x)
	grp <-  rep(1, length(namesS))
	grp[coreBiclusterSamples] <- 2
	matplot(y = t(x[coreBiclusterGenes , order(grp, decreasing=TRUE) ]), type = "n", xlab="",ylab="Expression", axes=FALSE,
			pch = rep(1, ncol(x)))
	matlines(y = t(x[coreBiclusterGenes , order(grp, decreasing=TRUE) ]), type = "l", 
			lty = rep(1, nrow(x)) , col = "black", lwd = 1, pch = 1)
	matlines(y = t(x[coreBiclusterGenes ,(grp==2) ]), type = "l",  
			lty = rep(1, nrow(x)) , col = "red", lwd = 1, pch = 1)
	# axis(1, 1:length(namesS), namesS[order(grp, decreasing = TRUE)], cex = 1, las=2)
	axis(2)
	box()
}


plotProfilesWithinBicluster <- function(x, main = "", sampleNames, geneNames = NULL){
	
	if (!length(sampleNames) == ncol(x))
		stop("'sampleNames' should be of length equal to the number of columns of 'x'")
	
	
	
	matplot(y = t(x), type = "n", xlab = "", ylab = "Expression", axes=FALSE,
			pch = rep(1, ncol(x)))
	title(main = main)
	axis(1, 1: length(sampleNames), sampleNames, cex=1, par(las=2))
	axis(2)
	box() #- to make it look "as usual"
	
	if (length(geneNames)){
		
		if (!length(geneNames) == nrow(x))
			stop("'geneNames' should be of length equal to the number of rows of 'x'")
		
		matlines(y=t(x), type = "l", lty=rep(1,5000) , col=1:length(geneNames),lwd = 1, pch = NULL)
		legend(x="top", legend= geneNames , col = 1:length(geneNames),
				border="none",ncol=3, lty=1, lwd=1,merge = TRUE)
	} else matlines(y=t(x), type = "l", lty=rep(1,5000) , col=rep(1,5000),lwd = 1, pch = rep(1,5000))
	
}
# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### SET COLOR VECTOR
#' color vector
#' @export
COLORS <- as.character(c(
				"red"
				,"darkgreen"
				,"blue"
				,"darkmagenta"
				,"darkorange"
				,"darkblue"
				,"black"
				,"brown"
				,"darkgrey"
				,"red"
				,"brown"
				,rev(colors())) ### if that's not enough
)

#' @export
#' @import Biobase
#' @importFrom graphics plot par grid lines points legend
.idOverviewPlots <- function(userOptions=userOptions
		,esetNorm=esetNorm
		,fileType=fileType
		,sqaPeptide=sqaPeptide
		,sqaProtein=sqaProtein){
	
	######################## OVERVIEW PLOT
	fdr <-userOptions$fdrCutoff
	nbPSM <- sum(!fData(esetNorm)$isFiltered,na.rm=T)
	cex=1
	par(cex.names=1.5, cex.axis= 1.25, cex.lab= 1.25, mar=c(5.1,3.1,4.1,1.1))
	if(userOptions$verbose) cat("IDENTIFICATIONS OVERVIEW PLOT \n")
	plot(0,0,type="n",xlim=c(0,8),ylim=c(-1,1.5), main=paste("\nIDENTIFICATIONS OVERVIEW\n (FDR ",fdr,")",sep=""), axes=F, xlab="", ylab="", cex=cex)
#text(0.5,2,paste("FDR",fdr),cex=cex, pos=4)
	xPos <- -0.5
	yPos <- 1
	
	if(fileType != "ProgenesisProtein"){
		text(xPos,yPos,paste(nbPSM," PSM"),cex=cex, pos=4)
		yPos <- yPos - 0.5
				
	}
		
	if(exists("sqaPeptide")){
		isMod <- nchar(as.character(fData(sqaPeptide$eset)$ptm)) > 0
		nbPeptides <- sum(!fData(sqaPeptide$eset)$isFiltered,na.rm=T)
		nbUnModPeptides <- sum(!fData(sqaPeptide$eset[!isMod,])$isFiltered,na.rm=T)
		nbModPeptides <- sum(!fData(sqaPeptide$eset[isMod,])$isFiltered,na.rm=T)
		nbProteins <-  length(unique(fData(sqaPeptide$eset)[!fData(sqaPeptide$eset)$isFiltered,]$proteinName))
		
		text(xPos,yPos,paste(nbPeptides, " PEPTIDES"),cex=cex, pos=4)
		yPos <- yPos - 0.25
		text(xPos+0.5,yPos,paste(nbUnModPeptides," UN-MOD."),cex=cex, pos=4)
		yPos <- yPos - 0.25
		text(xPos+0.5,yPos,paste(nbModPeptides," MOD."),cex=cex, pos=4)
		yPos <- yPos - 0.5
		
	
		
	}else if("peptide" %in% names(fData(esetNorm))) { ### TMT export
		text(xPos,yPos,paste(length(unique(fData(esetNorm)$peptide))," PEPTIDES" ),cex=cex, pos=4)
		yPos <- yPos - 0.5
	}
	
	if(exists("sqaProtein")){
		nbProteins <- sum(!fData(sqaProtein$eset)$isFiltered,na.rm=T)
	}
	text(xPos,yPos,paste(nbProteins," PROTEINS"),cex=cex, pos=4)
	par(mar=c(5.1,4.1,4.1,2.1))
	######################## OVERVIEW PLOT END
	
	### charge state
	if("charge" %in% names(fData(esetNorm))){
		if(userOptions$verbose) cat("CHARGE STATE PLOT \n")
		chargeTable <- table(fData(esetNorm)$charge[!fData(esetNorm)$isFiltered] )
		barplot2(chargeTable, xlab="Charge State", ylab="PSM Counts", col="blue", plot.grid = TRUE, grid.col="lightgrey")
			
	} 
	if(exists("sqaPeptide")){
		if(userOptions$verbose) cat("INFO: NB. MIS-CLEAVAGES PLOT \n")
		### mis-cleavage
		sel <- !fData(esetNorm)$isFiltered
		nbRows <- sum(sel,na.rm=T)
		nbSel <-min(c(nbRows,500))
		nbMCTable <- (table(getNbMisCleavages(fData(esetNorm)$peptide[sel][sample(nbRows,nbSel)] ))/nbSel)*100
		barplot2(nbMCTable, xlab="Nb. Mis-cleavages", ylab="Peptide Counts (%)", col="blue", plot.grid = TRUE, grid.col="lightgrey")
		
		### don't plot if many (more than 10) NA motifs (NA motid due to wrongly specified fasta)
		if("motifX" %in% names(fData(sqaPeptide$eset)) & (sum(is.na(fData(sqaPeptide$eset)$motifX)) < 10 ) ){ 
			
			if(userOptions$verbose) cat("INFO: MOTIF-X PLOT \n")
			motifTable <- table(.getUniquePtmMotifs(sqaPeptide$eset,format=(fileType == "ScaffoldTMT")+1)$ptm)
			
			if(nrow(motifTable) > 0){ # make sure some non NA motifs were found
				bp <- barplot2(motifTable, ,ylab="Modif. Site Counts", col="blue", plot.grid = TRUE, xaxt="n", grid.col="lightgrey")
				mtext(names(motifTable),side=1,at=bp[,1], line=0.2, las=2,cex=0.6)
			}	
					
		}else{
			
			if(fileType == "ScaffoldTMT"){
				
				ptmTag <- as.character(fData(sqaPeptide$eset)$ptm)[!fData(sqaPeptide$eset)$isFiltered]
				ptmTag[(nchar(ptmTag) == 0)] <- "Unmod" 
				ptmTag <- gsub("[0-9]","",unlist(strsplit(ptmTag,"\\, ")))
				ptmTable <- table(ptmTag[!fData(sqaPeptide$eset)$isFiltered] )
				
			}else{
				### ptm PROGENSIS 
				ptmTag <- as.character(fData(sqaPeptide$eset)$ptm)[!fData(sqaPeptide$eset)$isFiltered]
				ptmTag[nchar(ptmTag) == 0] <- "Unmod" 
				ptmTag <- gsub("\\[[0-9]*\\] {1,}","",unlist(strsplit(ptmTag,"\\|")))
				ptmTable <- table(ptmTag[!fData(sqaPeptide$eset)$isFiltered] )
				
			}
			if(userOptions$verbose) cat("PTM PLOT \n")	
			#barplot2(ptmTable, ylab="Peptide Counts", col="blue", plot.grid = TRUE, las=2, grid.col="lightgrey", cex.names=0.7)
			bp <- barplot2(ptmTable, ylab="Peptide Counts", col="blue", plot.grid = TRUE, xaxt="n", grid.col="lightgrey")
			mtext(names(ptmTable),side=1,at=bp[,1], line=0, cex=0.6, las=2)
			
		}
		
		if("nbPtmsPerPeptide"  %in% names(fData(sqaPeptide$eset))){
			if(userOptions$verbose) cat("INFO: NB.PTMS PER PEPTIDE PLOT \n")	
			### nb. ptms per peptide
			nbPtmPerPeptideTable <- table(fData(sqaPeptide$eset)$nbPtmsPerPeptide[!fData(sqaPeptide$eset)$isFiltered])
			barplot2(nbPtmPerPeptideTable, xlab="Nb. PTM per Peptide", ylab="Peptide Counts", col="blue", plot.grid = TRUE, grid.col="lightgrey")
		}
		
	}else if("peptide" %in% names(fData(esetNorm))){
		if(userOptions$verbose) cat("INFO: NB. MIS-CLEAVAGES PLOT 2 \n")	
		### mis-cleavage
		#nbMCTable <- table(getNbMisCleavages(fData(esetNorm)$peptide, protease="trypsin")[!fData(esetNorm)$isFiltered] )
		#barplot2(nbMCTable, xlab="Nb. Mis-cleavages", ylab="PSM Counts", col="blue", plot.grid = TRUE, grid.col="lightgrey")
		### mis-cleavage
		sel <- !fData(esetNorm)$isFiltered
		nbRows <- sum(sel,na.rm=T)
		nbSel <-min(c(nbRows,500))
		nbMCTable <- (table(getNbMisCleavages(fData(esetNorm)$peptide[sel][sample(nbRows,nbSel)] ))/nbSel)*100
		barplot2(nbMCTable, xlab="Nb. Mis-cleavages", ylab="Peptide Counts (%)", col="blue", plot.grid = TRUE, grid.col="lightgrey")
	
	}
	
	#### peptides per protein
	if(userOptions$verbose) cat("PEPTIDES PER PROTEIN PLOT\n")	
	if(exists("sqaProtein")){
		peptidesPerProtein <- fData(sqaProtein$eset)$nbPeptides[!fData(sqaProtein$eset)$isFiltered]
	}else{
		peptidesPerProtein <- fData(sqaPeptide$eset)$nbPeptides[!fData(sqaPeptide$eset)$isFiltered]
		peptidesPerProtein <- peptidesPerProtein[unique(names(peptidesPerProtein))]
	}
		
	### discard filtered out proteins
	peptidesPerProtein <- peptidesPerProtein[peptidesPerProtein > 0]
	#counts <- max(c(min(peptidesPerProtein,na.rm=T),1),na.rm=T):max(c(max(peptidesPerProtein,na.rm=T),2))
	xPeptides <- min(peptidesPerProtein,na.rm=T):max(c(max(peptidesPerProtein,na.rm=T),10))
	yCount <- unlist(lapply(xPeptides,function(t){ 
						sum(unlist(peptidesPerProtein) == t,na.rm=T) }))

	yCount[yCount == 0] <- NA 
	plot(xPeptides,yCount,type="n", log="x",xlab="Peptides Per Protein", ylab="Protein Counts")
	grid()
	lines(xPeptides,yCount,type="h",col="blue",lwd=2.5)
	
	### Id's ves Retention Time
	if(exists("sqaPeptide") && ("retentionTime"  %in% names(fData(sqaPeptide$eset)))) plotNbIdentificationsVsRT(sqaPeptide$eset)
	

}


### some quality control plots
#' @export
.qcPlots <- function(eset,selection=1:7,nbFeatures=500, userOptions=userOptions, ...){
	
	if(nrow(eset) < nbFeatures) nbFeatures <- nrow(eset) 
	sel <- sample(nrow(eset),nbFeatures,replace=F)
	
	par(mfrow=c(2,2), mar=c(5.1,4.1,4.1,2.1))
	if(1 %in% selection) {
		barplotMSSignal(eset, ...)
		plotMSSignalDistributions(log2(exprs(eset)), col=as.character(.getConditionColors(eset)[pData(eset)$condition,]), lwd=1.5, ...)
		boxplot(getAllCV(eset)*100,  col=as.character(.getConditionColors(eset)[pData(eset)$condition,]))
	
	}
	
	par(mfrow=c(1,1), mar=c(5.1,4.1,4.1,2.1))
	
	if(3 %in% selection)pairsAnnot(log2(exprs(eset)[sel,order(pData(eset)$condition)]), ...)
		
	if((4 %in% selection) && ("pMassError" %in% names(fData(eset)))){
			plotPrecMassErrorDistrib(eset,pMassTolWindow=userOptions$precursorMassFilter)
	}
	
	### retention time normalization plot
	if((5 %in% selection) && ("retentionTime" %in% names(fData(eset))) ){
		plotRTNormSummary(getRTNormFactors(eset, minFeaturesPerBin=100))
	}
	#if(6 %in% selection) invisible(hClustHeatMap(eset[sel & !fData(eset)$isFiltered,],main=paste("SELECTED NB. FEATURES:", nbFeatures), ...))
} 

### some quality control plots
#' @export
.idPlots <- function(eset,selection=1:7, qvalueThrs=0.01,...){
	
	if(!is.null(fData(eset)$idScore)){
		isDec <- isDecoy(fData(eset)$proteinName)
		scores <- fData(eset)$idScore
		qvalues <- fData(eset)$idQValue
	
		if(1 %in% selection) {
			plotScoreDistrib(scores[!isDec],scores[isDec], ...)
			
#			### score cutoffs abline
#			if("protein level"){
#				cutOff <- min(fData(eset)$idScore[fData(eset)$idQValue < qvalueThrs],na.rm=T )
#				if(is.finite(cutOff)) abline(v=cutOff, col="blue")
#			}else{
#				isMod <- nchar(as.character(fData(eset)$ptm)) > 0
#				modCutOff <- min(fData(eset[isMod,])$idScore[fData(eset[isMod,])$idQValue < qvalueThrs],na.rm=T )
#				noModCutOff <- min(fData(eset[!isMod,])$idScore[fData(eset[!isMod,])$idQValue < qvalueThrs],na.rm=T )
#				if(is.finite(modCutOff)) abline(v=modCutOff, col="blue")
#				if(is.finite(noModCutOff))abline(v=noModCutOff, "darkgreen")
#				legend("right",c("unmod. cutoff","mod. cutoff"),col=c("blue","darkgreen"), lty=1)
#			}
					
		}
		
		if(2 %in% selection) plotIdScoreVsFDR(scores,qvalues,qvalueThrs, ...)
		if(3 %in% selection) plotROC(qvalues,qvalueThrs,...)
	}
} 

#' @export
.getConditionColors <- function(eset){
	return(data.frame(colors=as.character(COLORS[1:length(levels(pData(eset)$condition))]), row.names=levels(pData(eset)$condition)))	
}

### color strip for volcano plot
#' @export
#' @importFrom graphics image mtext axis
.dotColorstrip <- function(colors,minSignal=0, maxSignal = maxSignal,lab= "C.V. (%)"  )
{
	bottom <- 1	
	count <- length(colors)
	m <- matrix(1:count, count, 1)
	image(m, col=colors, ylab="", axes=FALSE)
	
	labels <- round(seq(minSignal,maxSignal,length.out=3))
	
	### round off to nearest 5 
	labels <- as.character( round(labels / 5) *5 )
	
	at <- seq(0,1,length.out=3)
	
	axis(bottom,tick=TRUE, labels=labels , at=at )
	mtext(lab, bottom, adj=0.5, line=2)
}

### called from plotVolcano. Creates volcano plot form data.frame input
#' @export
#' @importFrom graphics plot par grid lines points
.plotVolcano <- function(d
		, ratioCutOffAbsLog2=0
		, absLog10pValueCutOff=2
		, xlim = range(d[,1],na.rm=T)
		, ylim = range(abs(        log10(d[,2])[is.finite(  log10(d[,2])  )]    )      ,na.rm=T) # pValue or qValue
		, maxSignal = max(d[,3][is.finite(d[,3])],na.rm=T)	
		#, higlightSel = rep(F,nrow(d))
		, controlCondition = "control"
		, caseCondition = "case"
		,...

){
	
	colorPalette <- rev(rich.colors(32))
	dotColors <- colorPalette[round(1+ d[,3]/maxSignal *31)] 
	
	par(fig=c(0,1,0.18,1))
	plot(0,0
			, xlab= paste("log2(",caseCondition,"/",controlCondition,")", sep="" )
			, ylab= paste("-log10(",names(d)[2],")",sep="")
			, type="n", ylim=ylim,xlim=xlim	
			,...
	)
	
	grid()
	
	### d[,2] qValues or pValues
	points(d[,1], abs(log10(d[,2])), col=dotColors, pch=20)
	
	### draw valid squares contianing valid proteins/peptides 
	lines(c(-ratioCutOffAbsLog2,-ratioCutOffAbsLog2),c(absLog10pValueCutOff,1000), col="grey")
	lines(c(ratioCutOffAbsLog2,ratioCutOffAbsLog2),c(absLog10pValueCutOff,1000), col="grey")
	lines(c(-1000,-ratioCutOffAbsLog2),c(absLog10pValueCutOff,absLog10pValueCutOff), col="grey")
	lines(c(ratioCutOffAbsLog2,1000),c(absLog10pValueCutOff,absLog10pValueCutOff), col="grey")
	
	# @TODO remove this feature
	#points(d$ratio[higlightSel],abs(log10(d[,2]))[higlightSel],col="darkgrey")
	
	### add heat color bar
	par(fig=c(0,1,0,0.3), new=TRUE)
	.dotColorstrip(colorPalette, maxSignal=maxSignal*100)
	par(mfrow=c(1,1))
	
}

#' @export
#' @importFrom graphics arrows
.errorBar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
		stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#' @export
#' @importFrom graphics legend plot par
.outerLegend <- function(...) {
	opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
			mar=c(0.5, 0.5, 0.5, 1), new=TRUE)
	on.exit(par(opar))
	plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
	legend(...)
}

#' @export
#' @import corrplot
#' @importFrom graphics legend
.correlationPlot <- function(d, textCol="black", labels=colnames(d),... ){
	
		
	corrplot(cor(d,use="complete")^2
			, type = "upper"
			, col = colorRampPalette(c("white","white","white","lightblue","blue"))(100)
			, cl.lim = c(0, 1)
			, tl.col = textCol
			#, method="circle"
			, method="pie"
			, addgrid.col="white" 
			,...
	)
	legend("bottomleft"
			, labels
			, fill=unique(textCol)
			#, fill=textCol
			, box.col="white"
	)
	legend("topright"
			, c("","",as.expression(bquote(R^2*"        ")))
			#, fill=unique(textCol)
			#, fill=textCol
			, box.col="white"
	)
}

#' @export
.allpValueHist <- function(sqa, col=as.character(.getConditionColors(sqa$eset)[names(sqa$pValue ),])){
	
	i <- 1
	for(cond in names(sqa$pValue ) ){
		hist(sqa$pValue[,cond],breaks=seq(0,1,length=50), main=cond, xlab="pValue",col=col[i])
		i <- i +1 
	}
}

#' Plots volcano, data points colored by max cv of the 2 compared conditions
#' @param obj safeQuantAnalysis object or data.frame
#' @param adjusted TRUE/FALSE plot qValues or pValues on y-axis
#' @param ratioThrs default 1
#' @param pValueThreshold default 0.01
#' @param ... see plot
#' @import Biobase gplots
#' @export
#' @note  No note
#' @details data.frame input object should contain 3 columns (ratio,qValue,cv)
#' @references NA
#' @examples print("No examples")
plotVolcano <- function(obj
		,ratioThrs=1
		,pValueThreshold=0.01
		,adjusted = T
		
		,...
){
	
	ratioCutOffAbsLog2 <- abs(log2(ratioThrs))
	absLog10pValueCutOff <- abs(log10(pValueThreshold))
	
	### plot volcanos for all case control comparisons
	if(class(obj) == "safeQuantAnalysis"){
		
		### need at least two condiotions
		if(length(unique(pData(obj$eset)$condition)) == 1){
			return(warning("INFO: Only one condition no plotVolcano \n"))
		}
				
		# ensure the same range on all volcano plots
		xlim <- range(obj$ratio, na.rm=T)
		
		if(adjusted){
			ylim <- range(abs(log10(obj$qValue)),na.rm=T)
		}else{
			ylim <- range(abs(log10(obj$pValue)),na.rm=T)
		}
		# ensure same scale on color legend
		cvMax <- max(as.vector(unlist(obj$cv)),na.rm=T)
		
		### avoid crash if no pValues (when no replicates)
		if(sum(!is.finite(c(ylim,cvMax))) > 0 ){
			cvMax <- 100
			ylim <- c(0,1)
		}
		
		controlCondition <- .getControlCondition(obj$eset)
		for(caseCondition in colnames(obj$pValue)){
			
			### create data.frame listing all data points (ratio,pvalue,cv)
			if(adjusted){
				
				d <- data.frame( ratio=obj$ratio[,caseCondition]
						,qValue=obj$qValue[,caseCondition]
						,cv=apply(obj$cv[,c(caseCondition,controlCondition)],1,max,na.rm=T)
				)
			}else{
				d <- data.frame( ratio=obj$ratio[,caseCondition]
						,pValue=obj$pValue[,caseCondition]
						,cv=apply(obj$cv[,c(caseCondition,controlCondition)],1,max,na.rm=T)
				)
			}
			
			.plotVolcano(d
					, ratioCutOffAbsLog2=ratioCutOffAbsLog2
					, absLog10pValueCutOff=absLog10pValueCutOff
					, xlim = xlim
					, ylim = ylim # pValue or qValue
					, maxSignal = cvMax	
					#, higlightSel = higlightSel
					, controlCondition = controlCondition
					, caseCondition = caseCondition
					
					,...)
		} 
	}else{
		.plotVolcano(obj
				, ratioCutOffAbsLog2=ratioCutOffAbsLog2
				, absLog10pValueCutOff=absLog10pValueCutOff
				#, higlightSel = higlightSel
				
				,...)
	}
	
}

#' Display experimental design, high-lighting the control condition 
#' @param eset ExpressionSet
#' @param condColors condition colors
#' @param version version number
#' @import Biobase
#' @importFrom graphics text plot
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotExpDesign <- function(eset, condColors=.getConditionColors(eset),  version="X"){
	
	### plot ctrl at the bottom
	pData(eset) <- rbind(pData(eset)[pData(eset)$isControl,],pData(eset)[!pData(eset)$isControl,])
	
	conditionNames <- as.character(unique(pData(eset)$condition))
	nbConditions <- length(conditionNames)
	controlCondition = .getControlCondition(eset)
	
	sampleNames <- rownames(pData(eset))
	nbSamples <- nrow(pData(eset))
	
	xlim <- c(-1,6)
	ylim <- c(-2,nbSamples+2)
	
	plot(0,0,type="n", xlim=xlim, ylim=ylim, main="Experimental Design", axes=FALSE, xlab="", ylab="")
		
	condYPosStep <- (nbSamples+2)/(nbConditions+1)
	sampleNb <- 1
	
	for(condNb in 1:length(conditionNames)){
		
		condName <- conditionNames[condNb]
		condCol = as.character(condColors[condName,])
		
		### control condition in  and underline
		if(condName == controlCondition){
			text(1,(condNb)*condYPosStep,bquote(underline(bold(.(condName)))), col=condCol)
		}else{
			text(1,(condNb)*condYPosStep,condName, col=condCol)
			
		}
		
		for(i in sampleNames[as.character(pData(eset)$condition) == condName]){
			text(4,sampleNb,paste(paste(sampleNb,":",sep=""), sampleNames[sampleNb]), col=condCol)
			sampleNb <- sampleNb + 1
		}
	}
	
	text(-1,-2,paste("SafeQuant v.", version), pos=4)
	
}

#' Plot lower triangle Pearsons R^2. Diagonal text, upper triangle all against all scatter plots with lm abline
#' @param data data.frame
#' @param textCol text color
#' @param diagText diagnoal text
#' @param col dot col
#' @param isHeatCol heat colors
#' @param cexTxt cex txt 
#' @param ... see plot 
#' @note  No note
#' @export
#' @importFrom graphics pairs abline par points
#' @importFrom stats cor
#' @details No details
#' @references NA
#' @examples print("No examples")
pairsAnnot<-
		function(data,textCol=rep(1,ncol(data)), diagText=colnames(data),col= rgb(0,100,0,50,maxColorValue=255),isHeatCol=F,cexTxt=2,...) {
	
	### we need at least two samples
	if(ncol(data.frame(data)) < 2 ){
		return(warning("INFO: Only one sample no pairsAnnot \n"))
	}
	
	
	### cex as a function of numbers of columns
	cex <- 0.8
	if(ncol(data) < 6){
		cex <- 2
	}else if(ncol(data) < 9){
		cex <- 1.4
	}
	else if(ncol(data) < 12){
		cex <- 1
	}
	
	count <- 1

	panel.lm <-
			function (x, y, col = par("col"), bg = NA, pch = par("pch"),
					#	cex = 1, col.lm = "red", lwd=par("lwd"), ...)
					col.lm = "red", lwd=par("lwd"))
	{
				
		if(isHeatCol){
			df <- data.frame(x,y)
			## Use densCols() output to get density at each point
			densityCol <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
			df$densityCol <- col2rgb(densityCol)[1,] + 1L
			## Map densities to colors
			cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
							"#FCFF00", "#FF9400", "#FF3100"))(256)
			
			# HACK to please CRAN CHECK "rollUp: no visible binding for global variable hCol"
			hCol <- NULL
			df$hCol <- cols[df$densityCol]
			
			points(y~x, data=df[order(df$densityCol),], pch = 20, col =hCol, bg = bg)
			col <<- "black"
			
		}else{
			points(x, y, pch = 20, col =col, bg = bg)
		}
		
		
		ok = is.finite(x) & is.finite(y)
		if (any(ok)){
			abline(lm(y~x,subset=ok), col = "blue", lwd=1.5)
			abline(coef=c(0,1),lty=2)
		}
	}
	
	panel.sse <-
			function(y, x, digits=2,...)
	{
		
		usr = par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		
		### discard non-fininte values
		ok = is.finite(x) & is.finite(y)
		x <- x[ok]
		y <- y[ok]
		
		model = summary(lm(y~x))
		r2= model$r.squared
		
		txt = round(r2, digits)
		txt = bquote(R^2*"" == .(txt))
		
		text(0.5, 0.5, txt,cex=cexTxt)
		
	}
	
	#panel.txt <- function(x, y, labels, cex, font,digits=2, ...){
	panel.txt <- function(x, y, labels, font,digits=2,...){

		txt = diagText[count]
		text(0.5, 0.6, txt, cex=cexTxt, col=textCol[count])
		
		count <<-count+1
	}
	
	pairs(data,lower.panel=panel.sse,upper.panel=panel.lm, text.panel=panel.txt,col=col,...)

	
}

#' Plot Percentage of Features with with missing values
#' @param eset ExpressionSet
#' @param col col
#' @param cex.axis cex.axis
#' @param cex.lab cex.lab
#' @param ... see plot 
#' @import Biobase
#' @importFrom graphics mtext
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
missinValueBarplot <- function(eset, col=as.character(.getConditionColors(eset)[pData(eset)$condition,]), cex.axis=1.25, cex.lab=1.25, ...){
	
	#par(mar=c(7.1,4.1,4.1,2.1))
	
	eset <- eset[!fData(eset)$isFiltered,]
	
	d <- apply(exprs(eset),2, function(t){ (sum(is.na(t)) / length(t))*100 } )
	ylim <- c(0,range(d)[2])
	if(max(d,na.rm=T) == 0) ylim <- c(0,100)
	
	bp <- barplot2( d
			, las=2
			, ylab="% Missing Values"
			, col=col
			, cex.axis = cex.axis
			, cex.lab=cex.lab
			, ylim=ylim
			, xaxt="n"
			, plot.grid=T
			, grid.col="lightgrey"		
			,...
	)

	mtext(names(d),side=1,at=bp[,1], las=2, line=0.3,cex=0.9)
	#par(mar=c(5.1,4.1,4.1,2.1))
	
}

### 
#' Barplot of ms-signal per column
#' @param eset expressionSet
#' @param col default condition colors
#' @param method c("median","sum","sharedSignal")
#' @param cex.lab default 1.25
#' @param cex.axis default 1.25
#' @param cex.names default 0.9
#' @param labels labels
#' @param ... see plot 
#' @import Biobase
#' @export
#' @importFrom graphics mtext
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
barplotMSSignal <- function(eset, col = as.character(.getConditionColors(eset)[pData(eset)$condition,])
	,method=c("sum","sharedSignal")
	,cex.lab=1.25
	,cex.axis=1.25
	,cex.names=0.9
	,labels=rownames(pData(eset))
	,...){
	
	sel <- !fData(eset)$isFiltered

	### only use feature qunatified in all runs for normalization
	if("sharedSignal" %in% method) sel <- sel & (as.vector(apply(is.finite(exprs(eset)),1,sum) == ncol(eset)))
	
	eset <- eset[sel,]			
		
	if("median" %in% method){
		profile <- apply(exprs(eset),2,median,na.rm=T)
		ylab <- "Median MS-Signal (Scaled)"
	}else{
		profile <- apply(exprs(eset),2,sum,na.rm=T)
		ylab <- "Summed MS-Signal (Scaled)"
	}
	profile <- profile/max(profile)
	
	bp <- barplot2(profile,las=2, col=col,ylab=ylab, xaxt="n"
			, plot.grid=T
			, grid.col="lightgrey"
			, cex.lab=cex.lab
			, cex.axis=cex.axis
			,...)
	mtext(labels,side=1,at=bp[,1], las=2, line=0.3,cex=cex.names)
	
}

#' C.V. boxplot 
#' @param eset ExpressionSet
#' @param col col
#' @param cex.names default 0.9
#' @param cex.axis default 1.25
#' @param cex.lab default 1.25
#' @param ylab C.V.
#' @param ... see plot 
#' @note  No note
#' @export
#' @importFrom graphics boxplot grid par
#' @details No details
#' @references NA
#' @examples print("No examples")
cvBoxplot <- function(eset,col=as.character(.getConditionColors(eset)[unique(pData(eset)$condition),]), cex.names=0.9,cex.axis=1.25,cex.lab=1.25,ylab="C.V. (%)",...){
	
	eset <- eset[!fData(eset)$isFiltered,]
	
	cv <- getAllCV(eset)
	### avoid crach when not enough repliocates
	if(sum(!is.na(cv)) > 0){
		boxplot(cv*100,yaxt="n",xaxt="n",...)

		grid(nx=NA, ny=NULL) #grid over boxplot
		par(new=TRUE)
		
		boxplot(cv*100,  col=col,las=2, ylab=ylab
				, cex.names=cex.names
				, cex.axis=cex.axis
				, cex.lab=cex.lab
				, textcolor=0
				, ...)
		
	}
}

#' Plot ms.signal distributions
#' @param d matrix of ms-signals
#' @param col color
#' @param ylab default "Frequnecy" 
#' @param xlab default "MS-Signal"
#' @param ... see plot 
#' @note  No note
#' @export
#' @importFrom graphics plot lines
#' @importFrom graphics plot hist abline legend grid mtext
#' @details No details
#' @references NA
#' @examples print("No examples")
plotMSSignalDistributions <- function(d, col=1:100,ylab="Frequnecy", xlab="MS-Signal",... ){
	

	### create a sms-signal histogram per column,
	breaks <- seq(min(d,na.rm=T),max(d,na.rm=T),length=40)
	mids <-  hist(d[,1],breaks=breaks,plot=F)$mids
	countPerCol <- data.frame(row.names=mids)
	ncol(countPerCol)
	
	for(c in colnames(d)){
		count <- hist(d[,c],breaks=breaks,plot=F)$count
		countPerCol <- cbind(countPerCol,count)
	}
	names(countPerCol) <-  colnames(d)
	
	#par(mar = c(5, 4, 4.1, 7.5))
	### plot histogram trend lines
	plot(mids,mids,ylim=range(countPerCol)
			, type="n"
			,xlab=xlab
			,ylab=ylab
			,...)
	for(i in 1:ncol(countPerCol)){
		lines(mids,countPerCol[,i],col=col[i],...)
	}
	#.outerLegend("right", names(countPerCol),lty=1 , col=col, bg="white")
	#par(mar = c(5, 4.1, 4.1, 2.1))
	
	#legend("topleft", names(countPerCol),lty=1 , col=col)
	
}


#' Hierarchical clustering heat map, cluster by runs intensity, features by ratio and display log2 ratios to control median
#' @param eset ExpressionSet
#' @param conditionColors data.frame of colors per condition
#' @param breaks default seq(-2,2,length=20)
#' @param dendogram see heatmap.2 gplots
#' @param legendPos see legend
#' @param ... see plot
#' @import gplots
#' @importFrom graphics legend
#' @importFrom stats hclust as.dist as.dendrogram
#' @export
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
hClustHeatMap <- function(eset
		,conditionColors =.getConditionColors(eset)
		,breaks=seq(-2,2,length=20)
		,dendogram = "column"
		,legendPos="left"
		,...
){
	
	# do not plot filtered
	eset <- eset[!fData(eset)$isFiltered,]
	
	### we need at least two conditions
	if(ncol(eset) == 1){
			return(cat("INFO: Only one condition no hClustHeatMap \n"))
	}
	
	#d <- log2(exprs(eset))
	### log2 ratios to median of control condition
	log2RatioPerMsRun <- log2(exprs(eset)) - log2(getSignalPerCondition(eset,method="median")[,.getControlCondition(eset)])
	
	feature.cor = cor(t(log2RatioPerMsRun), use="pairwise.complete.obs", method="pearson")
	feature.cor.dist = as.dist(1-feature.cor)
	feature.cor.dist[is.na(feature.cor.dist)] <- 0
	feature.tree = hclust(feature.cor.dist, method="ward")
	
	### !!!! clustering runs based on ratios is not a good idea as ratio corrrealtion of control runs is not expected to be high
	msrun.cor.pearson = cor(log2(exprs(eset)), use="pairwise.complete.obs", method="pearson")
	msrun.cor.pearson.dist = as.dist(1-msrun.cor.pearson)
	### to avoid error when replicates of the same condition are identical, DOES THIS EVER HAPPEN?
	msrun.cor.pearson.dist[is.na(msrun.cor.pearson.dist)] <- 0
	msrun.tree = hclust(msrun.cor.pearson.dist, method="ward")
	
	### sample colors
	samplecolors =  as.vector(unlist(conditionColors[pData(eset)$condition,]))
	labRow <- rownames(log2RatioPerMsRun)
	
	### do not display feature names if too many
	if(nrow(log2RatioPerMsRun) > 50){
		labRow <- rep("",(nrow(log2RatioPerMsRun)))
	}
	
	hm <- heatmap.2(as.matrix(log2RatioPerMsRun)
			, col=colorRampPalette(c(colors()[142],"black",colors()[128]))
			, scale="none"
			, key=TRUE
			, symkey=FALSE
			, trace="none"
			, cexRow=0.5
			, cexCol = 0.7
			,ColSideColors=samplecolors
			,labRow = labRow
			,Rowv=as.dendrogram(feature.tree)
			,Colv=as.dendrogram(msrun.tree)
			,dendrogram=dendogram
			,density.info="density"
			#,KeyValueName="Prob. Response"
			,breaks=breaks
			,na.rm=T
			, ...
	)
	
	legend(legendPos,levels(pData(eset)$condition), fill=as.character(conditionColors[,1]), cex=0.7, box.col=0)
	
}



### 
#' Plot Total Number of diffrentially Abundant Features (applying ratio cutoff) vs. qValue/pValue for all conditions
#' @param sqa SafeQuantAnalysis Object
#' @param upRegulated TRUE/FALSE select for upregulated features 
#' @param log2RatioCufOff log2 ratio cut-off
#' @param pvalCutOff pValue/qValue cut-off 
#' @param pvalRange pValue/qValue range
#' @param isLegend TRUE/FALSE display legend
#' @param isAdjusted TRUE/FALSE qValues/pValue on x-axis
#' @param ylab default  Nb. Features
#' @param ... see plot
#' @note  No note
#' @export
#' @importFrom graphics plot lines legend grid
#' @details No details
#' @references NA
#' @examples print("No examples")
plotNbValidDeFeaturesPerFDR <- function(sqa,upRegulated=T,log2RatioCufOff=log2(1),pvalRange=c(0,0.3)
		,pvalCutOff=1, isLegend=T,isAdjusted=T,ylab="Nb. Features", ... ){
	
	### we need at least two conditions
	if(length(unique(pData(sqa$eset)$condition)) == 1){
		return(cat("INFO: Only one condition no plotNbValidDeFeaturesPerFDR \n"))
	}
	
	if(isAdjusted){
		pvaluesPerCond <- sqa$qValue
		xlab= "False Discovery Rate (qValue)"
	}else{
		pvaluesPerCond <- sqa$pValue
		xlab= "pValue"
	}
	
	ratiosPerCond <- sqa$ratio
	conditionColors <- .getConditionColors(sqa$eset)
	
	pvalCutOffs <- seq(pvalRange[1],pvalRange[2], length.out=10)
	conditions <- names(pvaluesPerCond)
	
	### create data farme of roc curve per condition
	### store curves
	nbPassingCutOffsPerCond <- data.frame(row.names=pvalCutOffs)	
	
	### iterate over all conditions
	for(cond in conditions){
		
		pvals <- pvaluesPerCond[cond]	
		ratios <- ratiosPerCond[cond]
		
		#iterate over all cutoffs
		nbPassingCutOffs <- c()
		for(qCutOff in pvalCutOffs){
			
			if(upRegulated){
				nbPassingCutOffs <- c(nbPassingCutOffs, sum( (pvals < qCutOff) & (ratios >  log2RatioCufOff) ,na.rm=T) )
			}else{
				nbPassingCutOffs <- c(nbPassingCutOffs, sum( (pvals < qCutOff) & (ratios <  -log2RatioCufOff) ,na.rm=T ) )
			}
		}
		
		nbPassingCutOffsPerCond <- data.frame(nbPassingCutOffsPerCond, nbPassingCutOffs)
		
	}
	
	names(nbPassingCutOffsPerCond) <- conditions
	
	# plot roc curves	
	plot(0,0, ylim=c(0,max(nbPassingCutOffsPerCond)), xlim= c(min(pvalCutOffs) , max(pvalCutOffs)), type="n",ylab=ylab,xlab=xlab,  ...)
	grid()
	for(cond in conditions){
		lines(pvalCutOffs,nbPassingCutOffsPerCond[,cond], col= as.character(conditionColors[cond ,]), lwd=2)
	}
	abline(v=pvalCutOff,col="grey",lwd=1.5)
	if(isLegend){
		legend("topleft", conditions, fill=as.character(conditionColors[conditions,]), cex=0.6)
	}
	

}

#' Plot identifications target decoy distribution
#' @param targetScores target Scores
#' @param decoyScores decoy Scores
#' @param xlab default "Identification Score"
#' @param ylab default "Counts"
#' @param ... see plot
#' @note  No note
#' @export
#' @importFrom graphics plot hist points legend grid
#' @details No details
#' @references NA
#' @examples print("No examples")
plotScoreDistrib <-function(targetScores,decoyScores,xlab="Identification Score",ylab="Counts", ...){
	
	if(length(targetScores) >0 & length(decoyScores) >0){
		
		breaks <- hist(c(targetScores,decoyScores),plot=FALSE,breaks=100)$breaks
		
		targetHist = hist(targetScores, breaks=breaks, plot=FALSE)
		decoyHist = hist(decoyScores, breaks=breaks, plot=FALSE)
		
		ylim = c(0,max(c(decoyHist$counts,targetHist$counts),na.rm=T))
		
		plot(0,0,type='n', ylim=ylim, xlab=xlab, ylab=ylab, xlim=range(targetScores,decoyScores,na.rm=T),...)
		grid()
		points(targetHist$mids[targetHist$counts > 0], targetHist$counts[targetHist$counts > 0], col=1, type="h", lwd=4)
		points(decoyHist$mids[decoyHist$counts > 0], decoyHist$counts[decoyHist$counts > 0], col=2, type="h", lwd=5)
		
	}else{
		
		if(length(targetScores) >0){
			targetHist = hist(targetScores, breaks=100, plot=FALSE)
			plot(targetHist$mids, targetHist$counts, xlab=xlab, ylab=ylab,type="n",... )
			points(targetHist$mids[targetHist$counts > 0], targetHist$counts[targetHist$counts > 0], col=1, type="h", lwd=4)
		}
		
		cat("plotScoreDistrib: Not enough target or decoy scores\n")
	}
	
	legend("topright",c("target","decoy"),fill=c(1,2))
}



#' Plot FDR vs. identification score
#' @param idScore vector of identification scores 
#' @param qvals vector of q-valres
#' @param qvalueThrs threshold indicated by horizontal line
#' @param xlab default Identification Score
#' @param ylab default False Discovery Rate
#' @param ... see plot
#' @export
#' @importFrom graphics plot grid abline
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotIdScoreVsFDR <-function(idScore,qvals,qvalueThrs=0.01, ylab="False Discovery Rate", xlab="Identification Score",...){
	plot(sort(idScore),rev(sort(qvals)),type="l",ylab=ylab,xlab=xlab, ... )
	grid()
	abline(h=qvalueThrs,col="grey")
}



### plot identification scores ROC-curve, fdr vs. # valid identifications 
#' Plot Number of Identifications vs. FDR 
#' @param qvals vector of q-values
#' @param qvalueThrs threshold indicated by vertical line
#' @param breaks see breaks for hist function
#' @param xlab default "False Discovery Rate"
#' @param ylab default "Nb. Valid Identifications"
#' @param xlim default c(0,0.1)
#' @param col default blue
#' @param lwd default 1.5
#' @param ... see plot
#' @export
#' @importFrom graphics plot abline grid
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotROC <- function(qvals
		,qvalueThrs=0.01
		,xlab="False Discovery Rate"
		,ylab="Nb. Valid Identifications"
		,xlim=c(0,0.1)
		,breaks=100
		,col="blue"
		,lwd=1.5	
		,... ){
	
	if(length(breaks) == 1 ){
		breaks = seq(xlim[1],xlim[2],length=breaks)
	}else{  ### breaks override xlim
		xlim <- c(min(breaks),max(breaks))
	}
	
	validIds = c()
	for(fdr in breaks){
		validIds = c(validIds,sum(qvals < fdr, na.rm=T))
	}
	
	plot(breaks,sort(validIds), ylab=ylab,xlab=xlab, xlim=xlim,type="l",col=col,lwd=lwd,...)
	grid()
	abline(v=qvalueThrs,lwd=lwd,col="grey")
}



### Plot Precursor Mass Error Distribution 
#' Plot Precursor Mass Error Distribution 
#' @param eset ExpressionSet 
#' @param pMassTolWindow Precursor Mass Error Tolerance Window
#' @param ... see plot
#' @export
#' @importFrom graphics plot hist abline legend grid mtext axis
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotPrecMassErrorDistrib <- function(eset,pMassTolWindow=c(-10,10), ...){
	
#	par(mar = c(5, 4, 4.1, 3.5))
#	plot(rnorm(50), rnorm(50), col=c("steelblue", "indianred"), pch=20)
#	outerLegend("right", legend=c("Foo", "Bar"), pch=20, 
#			col=c("steelblue", "indianred"),
#			horiz=F, bty='n', cex=0.8)
	
	#par(mar=c(5.1,4.1,4.1,4.1))
	isDec <- isDecoy(fData(eset)$proteinName)
	pMassError <- fData(eset)$pMassError
	
	xRange 	<- range(pMassError,pMassTolWindow,na.rm=T)
	breaks <- seq(xRange[1],xRange[2],length=50)
	
	### decoy hist
	decoyMdHist 	<- hist(pMassError[isDec],breaks=breaks,plot=F)
	
	### target hist
	targetMdHist 	<- hist(pMassError[!isDec],breaks=breaks,plot=F)
	
	yRange <- range(decoyMdHist$counts,targetMdHist$counts)
	
	### plot histogram contours
	plot(decoyMdHist$mids,decoyMdHist$counts,type="l", col="red",ylim=yRange,xlab="mass diff. (ppm)",ylab="PSM Frequnecy",lwd=1.5, ...)
	grid()
	lines(decoyMdHist$mids,targetMdHist$counts,lwd=1.5)
	abline(v=pMassTolWindow[1],col="grey",lwd=2)
	abline(v=pMassTolWindow[2],col="grey",lwd=2)
	
	ratio <- (decoyMdHist$counts+1)/(targetMdHist$counts+1)
	ratio[ratio > 1] <- 1
	ratio <- ratio*yRange[2]
	
	lines(decoyMdHist$mids,ratio,lwd=1.5,col="blue",lty=2)
	
	### add ratio ratio line
	#axis(4,at=seq(0,yRange[2],length=5), labels=seq(0,1,length=5),col="blue", col.ticks="blue")
	axis(4,col="blue", col.ticks="blue", labels=NA, at=seq(0,yRange[2],length=5))
	mtext(seq(0,1,length=5),side=4,at=seq(0,yRange[2],length=5),col="blue", line=0.5)
	mtext("decoy-target PSM count ratio ",side=4, col="blue",line=1.5)
	
	legend("left",c("target","decoy"),  fill= c("black","red"))
	
	#par(mar=c(5.1,4.1,4.1,2.1))
	
}

### 
#' Plot precursorMass error v.s score highlighting decoy and displaying user specified user specified precursor mass filter
#' @param eset ExpressionSet
#' @param pMassTolWindow Precursor Mass Error Tolerance Window
#' @param ... see plot
#' @export
#' @importFrom graphics plot points abline legend grid
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
plotPrecMassErrorVsScore <- function(eset, pMassTolWindow=c(-10,10) ,...){
	
	pMassError <- fData(eset)$pMassError
	idScore <- fData(eset)$idScore
	isDec <- isDecoy(fData(eset)$proteinName) 
	
	withinTol <- (pMassError >= pMassTolWindow[1]) & (pMassError <= pMassTolWindow[2])
	
	plot(pMassError, idScore 
			, type="n"		
			#, col = isDecoy+1
			,ylab="score", xlab="mass diff. (ppm)", xlim=range(c(pMassTolWindow,pMassError),na.rm=T)
			,...
	)
	
	grid()
	
	points(pMassError[withinTol], idScore[withinTol], col="black", pch=20 )
	points(pMassError[!withinTol], idScore[!withinTol], col="grey", pch=20 )
	points(pMassError[isDec], idScore[isDec], col="red", pch=20 )
	
	abline(v=pMassTolWindow[1],col="lightgrey",lwd=1)
	abline(v=pMassTolWindow[2],col="lightgrey",lwd=1)
	
	legend("topright",c("target-kept","decoy","target-discarded"), pch=20, col= c("black","red","grey"))
	
}

#' Scatter plot with density coloring
#' @param x number vector
#' @param y number vector
#' @param isFitLm fit linear model
#' @param disp  c("abline","R","Rc") display options 
#' @param legendPos see legend
#' @param ... see plot
#' @import epiR
#' @importFrom graphics plot lines abline legend
#' @importFrom grDevices col2rgb colorRampPalette colors densCols rgb
#' @importFrom stats lm
#' @note  No note
#' @export
#' @references NA
#' @examples print("No examples")
plotXYDensity <- function(x,y,isFitLm=T,legendPos="bottomright",disp=c("abline","R","Rc"),  ...){
	
	df <- data.frame(x,y)
	
	## Use densCols() output to get density at each point
	densityCol <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
	df$densityCol <- col2rgb(densityCol)[1,] + 1L
	
	## Map densities to colors
	cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
					"#FCFF00", "#FF9400", "#FF3100"))(256)
	df$col <- cols[df$densityCol]
	
	## Plot it, reordering rows so that densest points are plotted on top
	plot(y~x, data=df[order(df$densityCol),], pch=20, col=col, cex=2, ...)
	
	### disp linerar model
	if(isFitLm){
		ok <- is.finite(x) & is.finite(y)
		fit <- lm(y[ok] ~ x[ok])
		
		if("abline" %in% disp){
			abline(coef=c(0,1),lty=2)
			abline(fit)
		}
		
		if("lowess" %in% disp){
			lines(lowess(x[ok],y[ok],...))
		}
		
		
		legd <- c()
		
		if("R" %in% disp) legd <- c(legd,as.expression(bquote(R^2*"" == .(round(summary(fit)$r.squared,2)))))
			
		if("Rc" %in% disp) legd <- c(legd,as.expression(bquote(R[c]*"" == .(round(epi.ccc(x,y)$rho.c$est,2)))))
		
		if(length(legd) > 0 ){
			legend(legendPos
					,legend= legd
					,text.col=1, box.col="transparent", cex=2
			)
		}
		
		return(fit)
		
	}
	
	return(NA)
	
}

#' Plot absolut Estimation calibration Curve
#' @param fit simple log-linear model
#' @param dispElements c("formula","lowess","stats")
#' @param cex.lab expansion factor for axis labels
#' @param cex.axis expansion factor for axis 
#' @param cex.text expansion factor for legend
#' @param cex.dot expansion factor for plotted dots
#' @param predictorName predictorName
#' @param text add names beside each dot
#' @param xlab xlab
#' @param ylab ylab
#' @param main main
#' @param ... see plot
#' @note  No note
#' @export
#' @importFrom graphics par plot lines mtext abline legend
#' @importFrom stats coef predict median
#' @references NA
#' @examples print("No examples")
plotAbsEstCalibrationCurve <- function(fit
		,dispElements = c("formula","lowess","stats")
		,xlab="Conc. (CPC) "
		,ylab="Pred. Conc. (CPC) "
		,predictorName = paste("log10(",names(coef(fit))[2],")",sep="")
		,text=F
		,cex.lab=1
		,cex.axis=1
		,cex.text=1
		,cex.dot=1
		,main = ""
		,...){
	x <- predict(fit) +  fit$residuals 
	y <- predict(fit) 			
	
	
	### some extra margin for axis labels
	par(mar=c(5.5,5.5,4.1,2.1))
	plot(10^x, 10^y,log="xy",xlab="",ylab="",main=main,cex.axis=cex.axis,cex=cex.dot,... )
	
	if(text){
		text(10^x, 10^y,rownames(fit$qr$qr))
	}
	
	
	abline(coef=c(0,1),lty=2)
	### add axis labels
	mtext(side=1,xlab,las=1, line=4,cex=cex.lab, ...)
	mtext(side=2,ylab,las=3, line=4,cex=cex.lab, ...)
	
	if( "formula" %in% dispElements){
		legend("bottomright"
				,paste("log10(",ylab,")"," = ", signif(coef(fit)[1],2)," + ",signif(coef(fit)[2],2)," * ",predictorName ,sep="")		
				,box.lwd=0
				,box.col="white"
				,cex=cex.text
				,...
		)
	}
	
	if( "lowess" %in% dispElements){
		lws <- lowess(y ~ x)
		lines(10^lws$x, 10^lws$y,col="red",...)
	}
	
	if( "stats" %in% dispElements){
		
		df <- data.frame(cpc =  x,signal = y)
		medianFoldError <- median(abs(getLoocvFoldError(df)[,1]),na.rm=T)
		
		legend("topleft"
				,legend=c(as.expression(bquote(R^2*"" == .(round(summary(fit)$r.squared,2))))
						,paste("Median Fold Error = ",round(medianFoldError,2))
				)
				#,text.col=c(1,2)
				,box.lwd=0
				,box.col="white"
				,cex=cex.text
				,...
		)
	}
	### reset margins
	par(mar=c(5.1,4.1,4.1,2.1))
}


#' Plot all retention time normalization profiles
#' @param eset ExpressionSet
#' @param ... see plot 
#' @param col condition colors 
#' @note  No note
#' @export
#' @importFrom graphics plot lines abline
#' @details No details
#' @seealso  \code{\link{getRTNormFactors}}
#' @references In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23589346} 
#' @examples print("No examples")
plotRTNormSummary <- function(eset, col = as.character(.getConditionColors(eset)[pData(eset)$condition,1]),...){
	
	rtNormFactors <- getRTNormFactors(eset, minFeaturesPerBin=100)
	ylim <- range(rtNormFactors,na.rm=T)
	runNames <- names(rtNormFactors)
	
	parDef <- par()
	par(mar = c(5, 4, 4.1, 7.5))
	plot(rownames(rtNormFactors),rtNormFactors[,1], type="n", ylab="log2(Ratio)",xlab="Retention Time (min)",ylim=ylim, ...)
	for(i in 1:ncol(rtNormFactors)){
		lines(rownames(rtNormFactors),rtNormFactors[,i], col=col[i], ...   )
	}
	abline(h=0, lty=2)
	.outerLegend("right", runNames, lty=1, col=col, lwd=2, bg="white")
	#.outerLegend("right", runNames, lty=1, col=col, lwd=1.5)
	par(parDef)
}

#' Plot all retention time profile overalying ratios
#' @param rtNormFactors data.frame of normalization factor per r.t bin and sample, obtained by getRTNormFactors
#' @param eset  ExprsssionSet
#' @param samples specify samples (sample numbers) to be plotted
#' @param main main
#' @param ... see plot see plot
#' @note  No note
#' @export
#' @importFrom graphics plot lines abline
#' @details No details
#' @seealso  \code{\link{getRTNormFactors}}
#' @references In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23589346} 
#' @examples print("No examples")
plotRTNorm <- function(rtNormFactors,eset,samples=1:ncol(rtNormFactors),main="", ...){
	
	### select for anchor proteins
	eset <- eset[fData(eset)$isNormAnchor,]
	
	# get all ratios to sample 1
	# @TODO How to select reference run?
	ratios <- log2(exprs(eset)) - log2(exprs(eset)[,1])
	#ratios <- log2(exprs(eset)) - log2(apply(exprs(eset),1,median,na.rm=T))
	
	for(samplesNb in samples){
		plot(fData(eset)$retentionTime,	ratios[,samplesNb]
				, col="lightgrey"
				, xlab="Retention Time (min)"
				, ylab="log2(Ratio)"
				, main=paste(main,names(rtNormFactors)[samplesNb]), ...)
		lines(as.numeric(rownames(rtNormFactors)),as.vector(unlist(rtNormFactors[,samplesNb])), ... )
		abline(h=0,lty=2,...)
	}
}


#' Plot the number of identified Features per Reteintion Time minute.
#' @param eset ExpressionSet
#' @param ... see plot see plot
#' @param cex.axis default 1.25
#' @param cex.lab default 1.25
#' @param col default "blue"
#' @param lwd default 2
#' @note  No note
#' @export
#' @importFrom graphics plot
#' @references NA
#' @examples print("No examples")
plotNbIdentificationsVsRT <- function(eset, cex.axis=1.25,cex.lab=1.25, col="blue", lwd=2, ...){
	rtTable <- table(round(fData(eset)$retentionTime))
	plot(as.numeric(names(rtTable)),rtTable[names(rtTable)], type="h", xlab="Retention Time (min)"
			, ylab="# Identified Features", cex.axis=cex.axis,cex.lab=cex.lab, col=col, lwd=lwd,...)
}


#' Plot qValue vs pValue
#' @param sqa SafeQuantAnalysis Object
#' @param lim x-axis and y-axis range
#' @param ... see plot
#' @note  No note
#' @export
#' @importFrom graphics plot lines abline legend
#' @details No details
#' @references NA
#' @examples print("No examples")
plotQValueVsPValue <- function(sqa, lim=c(0,1), ...){
	conditionColors <- .getConditionColors(sqa$eset)
	conditions <- names(sqa$pValue)
	
	plot(0,0,xlab="pValue",ylab="False Discovery Rate (qValue)", type="n", xlim=lim,ylim=lim)
	grid()
	for(cond in conditions){
		
		pVal <- sqa$pValue[,cond]
		qVal <- sqa$qValue[,cond]
		o <- order(pVal)
		lines(pVal[o],qVal[o],col= as.character(conditionColors[cond ,]), lwd=2)
	}
	
	abline(coef=c(0,1),lty=2)
	legend("bottomright", conditions, fill=as.character(conditionColors[conditions,]), cex=0.6)
}



### @TODO add unit test
#' @export
#' @importFrom graphics par plot lines mtext abline legend
.plotCalibrationCurve <- function(fit
		,dispElements = c("formula","lowess","stats")
		,xlab="Protein Copies/Cell Measured using SID"
		,ylab="Protein Copies/Cell Estimated using iBAQ"
		,cex=1.5
		,...){
	x <- predict(fit) + fit$residuals 						
	y <- predict(fit)		
	
	
	### some extra margin for axis labels
	par(mar=c(6.3,6.3,4.1,2.1))
	plot(10^x, 10^y
		,log="xy"
		,xlab=""
		,ylab=""
#		,yaxt="n"
#		,xaxt="n"
		,cex=cex
		#,lwd=cex
		,pch=19
		,las=2,cex.axis=cex
		,cex.main=cex
		,... )
	
#	axis(1,las=2,cex.axis=cex)
#	axis(2,las=2,cex.axis=cex)

	abline(coef=c(0,1),lty=2)
	### add axis labels
	mtext(side=1,xlab,las=1, line=4.8, cex=cex, ...)
	mtext(side=2,ylab,las=3, line=4.8, cex=cex, ...)
	
	if( "formula" %in% dispElements){
		legend("bottomright"
				,paste("log10(Est. CPC)"," = ", signif(coef(fit)[1],2)," + ",signif(coef(fit)[2],2)," * log10(iBAQ)" ,sep="")		
				,box.lwd=0
				,box.col="white"
				,cex=cex-0.4
				,...
		)
	}
	
	if( "lowess" %in% dispElements){
		lws <- lowess(y ~ x)
		lines(10^lws$x, 10^lws$y,col="red",...)
	}
	
	if( "stats" %in% dispElements){
		
		df <- data.frame(cpc =  x,signal = y)
		medianFoldError <- median(abs(getLoocvFoldError(df)[,1]),na.rm=T)
		linRc <- as.vector(unlist(epi.ccc(x,y)$rho.c)[1])
	
		legend("topleft"
				,legend=c(as.expression(bquote(R^2*"" == .(round(summary(fit)$r.squared,2))))
							#,paste("Lin's Rc  = ",round(linRc,2))
							,as.expression(bquote(R[c]*"" == .(round(linRc,2))))
							,paste("Median Fold Error = ",round(medianFoldError,2))
				)
				#,text.col=c(1,2)
				,box.lwd=0
				,box.col="white"
				,cex=cex
				,...
		)
	}
	### reset margins
	par(mar=c(5.1,4.1,4.1,2.1))
}


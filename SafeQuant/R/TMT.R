# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

#' Get Thermo TMT impurity matrix
#' @param plexNb integer, 6 or 10 plex 
#' @return impurity matrix matrix
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getImpuritiesMatrix <- function (plexNb=6){
#getImpuritiesMatrix <- function (plexNb=6, test=F){
	
#	if(test){
#		
#		#### test from MSnbase
#		tmt126 <- c(0,0,6.1,0)
#		tmt127 <- c(0,0.5,6.7,0)
#		tmt128 <- c(0,1.1,4.2,0)
#		tmt129 <- c(0,1.7,4.1,0)
#		tmt130 <- c(0,1.6,2.1,0)
#		tmt131 <- c(0.2,3.2,2.8,0)
#		
#		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0),c(0,tmt127,0),c(0,tmt128,0),c(0,tmt129,0),c(0,tmt130,0),c(0,tmt131,0)) ) )/100
#		rownames(reporterIonIsotopicDistrib) <- 1:6
#		
#	}else
	if(plexNb == 10){
		
		### FROM THERMO PRODUCT DATA SHEET
		
		# 10-plex
		tmt126 <- c(0,0,4.69,0)
		tmt127N <- c(0,0.4,6.5,0)
		tmt127C <- c(0,0.2,4.6,0.3)
		tmt128N <- c(0,0.9,4.7,0.2)
		tmt128C <- c(0.1,0.53,2.59,0)
		tmt129N <- c(0,0.73,2.49,0)
		tmt129C <- c(0,1.3,2.5,0)
		tmt130N <- c(0,1.2,2.8,2.7)
		tmt130C <- c(0.1,2.9,2.9,0)
		tmt131 <- c(0,2.36,1.43,0)
		
		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0)
			,c(0,tmt127N,0),c(0,tmt127C,0)
			,c(0,tmt128N,0),c(0,tmt128C,0)
			,c(0,tmt129N,0),c(0,tmt129C,0)
			,c(0,tmt130N,0),c(0,tmt130C,0)
			,c(0,tmt131,0)) ) )/100
		rownames(reporterIonIsotopicDistrib) <- 1:10
		
	}else if(plexNb == 6){
		
		### FROM THERMO PRODUCT DATA SHEET
		# 6-plex old
#		tmt126 <- c(0,0,9.4,0.6)
#		tmt127 <- c(0.1,0.8,8.6,0.4)
#		tmt128 <- c(0.1,1.4,7.2,0.5)
#		tmt129 <- c(0.1,1.6,5.7,0.3)
#		tmt130 <- c(0.2,2.1,5.1,0.2)
#		tmt131 <- c(0.1,4.1,4.7,0.1)
		
#		tmt126 <- c(0,0,5.85,0.0)
#		tmt127 <- c(0.0,0.4,6.5,0)
#		tmt128 <- c(0.34,1.2,3.59,0)
#		tmt129 <- c(0,1.59,3.48,0)
#		tmt130 <- c(0.17,2.9,2.18,0)
#		tmt131 <- c(0.35,3.58,2.1,0)
		
		tmt126 <- c(0,   0  , 5.6  ,0)
		tmt127 <- c(0,   0.4, 5    ,0)
		tmt128 <- c(0,   1.2, 5.3  ,0.1)
		tmt129 <- c(0.2, 1.5, 4.1  ,0)
		tmt130 <- c(0.1, 2.6, 2.5  ,0)
		tmt131 <- c(0.1 ,3.3, 2.6  ,0)
		
		#### test from MSnbase
		#tmt126 <- c(0,0,6.1,0)
		#tmt127 <- c(0,0.5,6.7,0)
		#tmt128 <- c(0,1.1,4.2,0)
		#tmt129 <- c(0,1.7,4.1,0)
		#tmt130 <- c(0,1.6,2.1,0)
		#tmt131 <- c(0.2,3.2,2.8,0)
		
		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0),c(0,tmt127,0),c(0,tmt128,0),c(0,tmt129,0),c(0,tmt130,0),c(0,tmt131,0)) ) )/100
		rownames(reporterIonIsotopicDistrib) <- 1:6
		
	}else if(plexNb == -6){	 ### test
		
		tmt126 <- c(0,   0  , 0  ,0)
		tmt127 <- c(0,   0, 0    ,0)
		tmt128 <- c(0,   0, 0  ,0)
		tmt129 <- c(0, 0, 4.1  ,0)
		tmt130 <- c(0, 0, 0  ,0)
		tmt131 <- c(0 ,3.3, 0  ,0)
		#tmt131 <- c(0 ,0, 0  ,0)
		
		plexNb <- 6
		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0),c(0,tmt127,0),c(0,tmt128,0),c(0,tmt129,0),c(0,tmt130,0),c(0,tmt131,0)) ) )/100
		rownames(reporterIonIsotopicDistrib) <- 1:6
	
	}else if(plexNb == -60){	 ### test
		
		tmt126 <- c(0,   0, 0  ,0)
		tmt127 <- c(0,   0, 5  ,0)
		tmt128 <- c(0,   0, 0  ,0)
		tmt129 <- c(0, 1.5, 4.1,0)
		tmt130 <- c(0, 0, 0  ,0)
		tmt131 <- c(0 ,3.3, 0  ,0)
		#tmt131 <- c(0 ,0, 0  ,0)
		
		plexNb <- 6
		reporterIonIsotopicDistrib <-  t( as.matrix( data.frame(c(0,tmt126,0),c(0,tmt127,0),c(0,tmt128,0),c(0,tmt129,0),c(0,tmt130,0),c(0,tmt131,0)) ) )/100
		rownames(reporterIonIsotopicDistrib) <- 1:6
		
	}
	
	###convert  reporterIonIsotopicDistrib -> impuritiesMatrix
	diagonals <- 1-apply(reporterIonIsotopicDistrib,1,sum)
	impuritiesMatrix <- diag(diagonals)
	
	for(i in 1:plexNb){
		
		### get affected channels
		affectedChannels <- c(-3,-2,-1,1,2,3) + i
		
		for(j in 1:6){
			affectedChannel <- 	affectedChannels[j]	
			if((affectedChannel > 0) & (affectedChannel <= plexNb)){
				impuritiesMatrix[i,affectedChannel] <- reporterIonIsotopicDistrib[i,j]
			}
		}
	}
	
	
	
	#diag(impuritiesMatrix) <- 1
	impuritiesMatrix <- t(impuritiesMatrix) 
	return(impuritiesMatrix)
	
}


#' Correct channel intensities based on Reporter ion Isotopic Distributions 
#' @param tmtData data.frame containing tmt channel intensities
#' @param impurityMatrix correction matrix
#' @return data.frame of corrected tmt intensities
#' @export
#' @note  No note
#' @details Same method as MSnbase, and described in Breitwieser et al. 2012 (Book Chapter)
#' @references NA 
#' @examples print("No examples")
purityCorrectTMT <- function(tmtData,impurityMatrix=impurityMatrix ){
	tmtDataCorrected <- matrix(nrow=nrow(tmtData),ncol=ncol(impurityMatrix))
	### solve linear system (spectrum by spectrum)
	for(i in 1:nrow(tmtData)){
		correctedSpectrumIntensities <- solve(impurityMatrix, tmtData[i,])
		tmtDataCorrected[i,1:ncol(impurityMatrix)] <- correctedSpectrumIntensities
	}
	
	tmtDataCorrected[!is.na(tmtDataCorrected) & (tmtDataCorrected < 0)] <- NA
	
#	if(invalidReplace  == "allZero" ){
#		tmtDataCorrected[!is.na(tmtDataCorrected) & (tmtDataCorrected < 0)] <- 0
#	}else if(invalidReplace  == "allNA"){
#		tmtDataCorrected[!is.na(tmtDataCorrected) & (tmtDataCorrected < 0)] <- NA
#	}else if(invalidReplace  == "allOrg"){ ### keep original values
#		tmtDataCorrected[!is.na(tmtDataCorrected) & (tmtDataCorrected <= 0)] <- tmtData[!is.na(tmtDataCorrected) & (tmtDataCorrected <= 0)] 
#	}### else Negative values are allowed
	
	colnames(tmtDataCorrected) <- colnames(tmtData)
	
	return(tmtDataCorrected)
}


#' Create Experimental Design
#' @param tag user input tag e.g. 1,2,3:4,5,6 indicating two condition with 3 reps each
#' @param nbPlex tmt 6  or 10 plex
#' @return expDesign data.frame
#' @export
#' @details The first listed condition is always the control condition
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
createExpDesign <- function(tag,nbPlex){
	
	sampleOrder <- as.numeric(unlist(strsplit(tag,"[\\,\\:]")))
	
	# make sure no duplicates, withing range etc.
	if(is.na(sampleOrder[1]) 
			| (max(table(sampleOrder))>1) 
			| (length(sampleOrder) > 10)
			| (max(sampleOrder) > 10)
			| (min(sampleOrder) < 1)
			| (length(sampleOrder) > nbPlex)
			){
		cat("ERROR: getExpDesign, INVALID EXPERIMENTAL DESIGN",tag,"\n")
		quit("no")
	}
	
	expDesign <- data.frame(row.names=sampleOrder,condition=rep("Ctrl",length(sampleOrder)),isControl=rep(F,length(sampleOrder) ))
	
	# create vector describing condition grouping e.g. c(3,3) ctrl 3 samples, 3 samples cond1
	condGrouping <- c()
	for(cond in unlist(strsplit(tag,":"))){
		condGrouping <- c(condGrouping,length(unlist(strsplit(cond,","))))
	}	
	
	expDesign$condition <- as.character(expDesign$condition)
	
	#browser()
	# update expDesign data.frame
	condNb <- 1
	sampleStartIdx <- 1 
	for(sampleCount in condGrouping ){
		
		if(condNb == 1){
			expDesign$condition[1:sampleCount] <- rep("Ctrl",sampleCount)
			expDesign$isControl[1:sampleCount] <- rep(T,sampleCount)
		}else{
			expDesign$condition[sampleStartIdx:(sampleStartIdx+sampleCount-1)] <- rep(paste("cond",condNb-1,sep="_"),sampleCount)
			expDesign$isControl[sampleStartIdx:(sampleStartIdx+sampleCount-1)] <- rep(F,sampleCount)
		}
		condNb <- condNb+1
		sampleStartIdx <-  sampleStartIdx+sampleCount
	}	
	
	return(expDesign)
}


#' Sum up raw intensities per protein and channel. keep track of number of summed spectra and unique peptides
#' @param intData data.frame of intensities per channel
#' @param proteinACs vector of protein accession numbers
#' @param peptides vector of peptide sequneces
#' @param minNbPeptPerProt minimal number of peptides per protein
#' @return list containing  3 objects 1) data.frame of channel intensities per protein ac, 2) vector listing number of summed spectra per protein, 3) vector listing number of summed peptides per protein
#' @export
#' @details NA
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getIntSumPerProtein  <- function(intData,proteinACs,peptides,minNbPeptPerProt=1){
	uniqueProteinACs <- unique(proteinACs)
	nbProteins <- length(uniqueProteinACs)
	
	perProteinIntSum <- matrix(nrow=nbProteins, ncol=ncol(intData))
	
	c <-1
	### init progressio bar
	pbSum <- txtProgressBar(min = 0, max = nbProteins, style = 3)
	
	spectraPerProtein <- c()
	peptidesPerProtein <- c()
	for(ac in uniqueProteinACs){
		
		# disp/increment progression bar
		setTxtProgressBar(pbSum, c)
		
		
		sel <- proteinACs %in% ac 
		sset 	<- intData[sel,]
		uPept <- unique(peptides[sel])
		
		if(sum(sel) > 1){
			perProteinIntSum[c,] <- apply(sset,2,sum)
		}else{
			perProteinIntSum[c,] <- as.vector(unlist(sset))
		}
		
		spectraPerProtein <- c(spectraPerProtein,nrow(sset))
		peptidesPerProtein <- c(peptidesPerProtein,length(unique(peptides[sel])))
		
		c <- c+1
	}
	# terminate progession bar
	close(pbSum)
	cat("\n")
	
	perProteinIntSum  <- data.frame(perProteinIntSum)
	rownames(perProteinIntSum) <- uniqueProteinACs
	names(perProteinIntSum) <- paste("channel_",1:ncol(intData),sep="")
	
	### FILTER based on min peptides per protein
	ok <- peptidesPerProtein >= minNbPeptPerProt
	if(sum(ok) < 2 ){
		cat("ERROR: getIntSumPerProtein. Not enough Features, Review filtering parameter minNbPeptPerProt ", minNbPeptPerProt,"\n")
		quit(status=-1)
	}
		
	perProteinIntSum <- perProteinIntSum[ok,]
	spectraPerProtein <- spectraPerProtein[ok]
	peptidesPerProtein <- peptidesPerProtein[ok]

	ret <- list()
	ret$perProteinIntSum <- perProteinIntSum
	ret$spectraPerProtein <- spectraPerProtein
	ret$peptidesPerProtein <- peptidesPerProtein
	
	return(ret)
	
}




## get index of max (created for data.table)
#' get index of max in vecotr of numeric values
#' @param v vector
#' @export
getMaxIndex <- function(v){
	if(all(is.na(v))) return(1)
	return(which(max(v,na.rm=T) == v)[1])	
}

## return first entry per column
#' @export
.getFirstEntry <- function(x){
	return(x[1])
}

#' @export
.log2Exprs <- function(eset){
	exprs(eset) <- log2(exprs(eset))
	return(eset)
}

#' @export
.exp2Exprs <- function(eset){
	exprs(eset) <- 2^(exprs(eset))
	return(eset)
}

#' @export
.getControlCondition <-function(eset){
	return(unique(as.character(pData(eset)$condition[pData(eset)$isControl]))[1])
}

#' Create ExpressionSet object
#' @param expressionMatrix matrix of expression signals per feature and sample
#' @param expDesign experimental design data.frame
#' @param featureAnnotations data.frame including e.g: Protein Description, Id score etc.
#' @return ExpressionSet object
#' @import Biobase 
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
createExpressionDataset <- function(expressionMatrix=expressionMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations){
	
	### make sure that only one unique condition is specified as control
	if(length(unique(as.character(expDesign$condition[expDesign$isControl])) ) > 1){
		stop("Invalid experimental design")	
	}
	
	### phenoData: stores expDesign
	# display with pData(eset): 
	#	condition    isControl
	#	A_rep_1         A T
	#	A_rep_2         A T
	#	B_rep_1         B    F
	#	B_rep_2         B    F
	#	C_rep_1         C    F
	#	C_rep_2         C    F
	#pData <- new("AnnotatedDataFrame", data=expDesign)
	pData <- AnnotatedDataFrame(data=expDesign)
	
	### featureData:add more data to each feature. E.g: Protein Description, Id score etc.
	
	return(ExpressionSet(assayData = expressionMatrix
					, phenoData =  pData							### yeah this is weird, but gives error if rolling up already 
					# rolled up eset unless colindices are explicitly specified		
					#, featureData = new("AnnotatedDataFrame", data= featureAnnotations[,1:ncol(featureAnnotations)])  
					, featureData = AnnotatedDataFrame(data= featureAnnotations[,1:ncol(featureAnnotations)])  
			)
	)
}


#' Perform statistical test (mderated t-test), comparing all case to control
#' @param eset ExpressionSet
#' @param method c("pairwise","all") 
#' @param adjust TRUE/FALSE adjust for multiple testing using Benjamini & Hochberg  (1995) method 
#' @param log T/F log-transform expression values
#' @return ExpressionSet object
#' @export
#' @import limma Biobase
#' @importFrom stats model.matrix p.adjust
#' @note  No note
#' @details No details
#' @references Empirical Bayes method, Smyth (2004), \url{http://www.ncbi.nlm.nih.gov/pubmed/16646809} 
#' @seealso \code{\link[limma]{eBayes}}
#' @examples print("No examples")
getAllEBayes <- function(eset=eset, adjust=F, log=T, method="pairwise"){
	
	controlCondition <- .getControlCondition(eset)
	caseConditions <- setdiff(unique(pData(eset)$condition) ,controlCondition)
	
	if(log)(exprs(eset) <- log2(exprs(eset)))
	
	#### NON-PAIR-WISE COMPARISONS
	if("all" %in% method){ 

		# Example 4,5,6 are control		
		#		> design
		#		(Intercept) f2 f3 f4 f5 f6
		#		1            1  1  0  0  0  0
		#		2            1  1  0  0  0  0
		#		3            1  1  0  0  0  0
		#		4            1  0  0  0  0  0
		#		5            1  0  0  0  0  0
		#		6            1  0  0  0  0  0
		#		7            1  0  1  0  0  0
		#		8            1  0  1  0  0  0
		#		9            1  0  1  0  0  0
		#		10           1  0  0  1  0  0
		#		11           1  0  0  1  0  0
		#		12           1  0  0  1  0  0
		#		13           1  0  0  0  1  0
		#		14           1  0  0  0  1  0
		#		15           1  0  0  0  1  0
		#		16           1  0  0  0  0  1
		#		17           1  0  0  0  0  1
		#		18           1  0  0  0  0  1
		#		attr(,"assign")
		#		[1] 0 1 1 1 1 1
		#		attr(,"contrasts")
		#		attr(,"contrasts")$f
		#		[1] "contr.treatment"
		
		condToNb <- data.frame(row.names=unique(pData(eset)$condition))
		condToNb[as.character(pData(eset)$condition[pData(eset)$isControl])[1],1] <- 1 
		condToNb[unique(as.character(pData(eset)$condition[!pData(eset)$isControl])),1] <- 2:nrow(condToNb)
		nbToCond <- data.frame(row.names=condToNb[,1],rownames(condToNb))
		
		f <- factor(condToNb[eset$condition,])
		design <- model.matrix(~f)
		#### calculate modified t-statistic, empirical Bayes method, Smyth (2004) 
		fit <- eBayes(lmFit(eset,design))
		pvalues <- data.frame(fit$p.value[,2:ncol(fit$p.value)])
		
	}else{ 	#### PAIR-WISE COMPARISONS "pairwise" %in% method
		pvalues <- data.frame(row.names=featureNames(eset))
		
		for(cC in caseConditions){
			
			### at least one replicate of one condition requires
			selCol <- unlist(pData(eset)$condition %in% c(controlCondition,cC))
			if(sum(selCol) > 2 ){
				
				esetPair <- eset[,selCol]
				
				f <- factor(as.character(esetPair$condition))
				design <- model.matrix(~f)
				#### calculate moderated t-statistic, empirical Bayes method, Smyth (2004) 
				fit <- eBayes(lmFit(esetPair,design))
				
				p <- fit$p.value[,2]
				
			}else{
				p <- rep(NA,nrow(eset))
			}
			pvalues <- cbind(pvalues,p)
		}
		names(pvalues) <- caseConditions
			
	}
		
	if(adjust){ ### adjust for multiple testing using Benjamini & Hochberg (1995) method 
		for(i in 1:ncol(pvalues)){
			pvalues[,i] <-p.adjust(pvalues[,i],method="BH")
		}
	}
	
	return(pvalues)	
	
}

#' Calculate ratios, comparing all case to control
#' @param eset ExpressionSet 
#' @param method median or mean
#' @param log2 transform
#' @return ExpressionSet object
#' @export
#' @import Biobase 
#' @note  No note
#' @details No details
#' @references NA
#' @examples print("No examples")
getRatios <- function(eset, method="median", log2=T){
	
	controlCondition <- .getControlCondition(eset)
	caseConditions <- setdiff(unique(pData(eset)$condition) ,controlCondition)
	
	ratios <- data.frame(row.names=featureNames(eset))
	
	allConditions <-  pData(eset)$condition
	
	for(cC in caseConditions){
		
		mCase <- exprs(eset)[,allConditions == cC]
		mControl <- exprs(eset)[,allConditions == controlCondition]
		
		if(method == "median"){
			if(sum(allConditions == cC) > 1){
				mCase <- unlist(apply(exprs(eset)[,allConditions == cC],1,median, na.rm=T))
			}
			if(sum(allConditions == controlCondition) > 1){
				mControl <- unlist(apply(exprs(eset)[,allConditions == controlCondition],1,median, na.rm=T))
			}
		}else if(method == "mean"){
			if(sum(allConditions == cC) > 1){
				mCase <- unlist(apply(exprs(eset)[,allConditions == cC],1,mean, na.rm=T))
			}
			if(sum(allConditions == controlCondition) > 1){
				mControl <- unlist(apply(exprs(eset)[,allConditions == controlCondition],1,mean, na.rm=T))
			}
		}else{
			stop("Unknown method ",method)
		}
		ratios <- cbind(ratios,mCase/ mControl)
	}
	
	names(ratios) <- caseConditions
	
	if(log2){
		return(log2(ratios))
	}else{
		return(ratios)
	}
}

#' Calculate Coefficiant of Variance per feature (Relative standard Deviation) per Condition
#' @param eset ExpressionSet
#' @return data.frame of CVs per condition
#' @import Biobase
#' @export
#' @note  No note
#' @details CV = sd / mean 
#' @references NA
#' @seealso  \code{\link{getCV}}
#' @examples print("No examples")
getAllCV <- function(eset){
	
	allConditions <- unique(pData(eset)$condition)
	cv <- data.frame(row.names=featureNames(eset))
	for(cond in allConditions){
		
		cvTmp <- rep(NA,nrow(cv))
		colMatch <- pData(eset)$condition ==  cond
		
		### if replicates
		if(sum(colMatch) > 1){
			cvTmp <- getCV(exprs(eset)[,colMatch])
		}
		
		cv <- cbind(cv,cvTmp)
	}
	names(cv) <- allConditions
	return(cv)
}

#' Get signal at zscore x (x standard deviations below mean)
#' @param intensities refrence run signals
#' @param promille baseline value set as specified promille
#' @return baseline value
#' @export
#' @importFrom stats quantile
#' @note  No note
#' @references NA
#' @examples print("No examples")
getBaselineIntensity <- function(intensities , promille = 5){
	
	### 
	#intensities <- sample(intensities[!is.na(intensities)],1000,replace=T)
	if((promille < 1) |  (promille > 1001)){
		stop("Invalid percentile: ", promille)
	}
	
	intensities <- intensities[is.finite(intensities)]
	
	suppressWarnings(return(quantile(intensities,probs = seq(0,1,0.001),na.rm=T)[promille+1]))
	
}

#' Calculate Coefficiant of Variance per feature (Relative standard Deviation) 
#' @param data data.frame of replicate signals
#' @return vector of CVs
#' @export
#' @importFrom stats sd
#' @note  No note
#' @details CV = sd / mean 
#' @references NA 
#' @examples print("No examples")
getCV <- function(data){
	return(apply(data,1,sd,na.rm=T)/apply(data,1,mean,na.rm=T))
}

#' Get retentiontime base normalization factors
#' @param eset ExpressionSet
#' @param minFeaturesPerBin  minumum number of features per bin. If nb. features are < minFeaturesPerBin -> include neighbouring bins.
#' @return data.frame normalization factors per retention time bin (minute)
#' @import limma Biobase
#' @export
#' @note  No note
#' @details No details
#' @references In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23589346} 
#' @examples print("No examples")
getRTNormFactors <- function(eset, minFeaturesPerBin=100){
	
	### check if eset contain necessary retentionTime colummn
	if(is.null(fData(eset)$retentionTime)){
		stop("retentionTime not found")
	}
	
	### check if isNormAnchor and isFiltered columns are defiend, if -> get normalization factors from nonFiltered anchor proteins
#	if(!is.null(fData(eset)$isNormAnchor) & !is.null(fData(eset)$isFiltered)){
#		sel <- fData(eset)$isNormAnchor & !fData(eset)$isFiltered
#		if(sum(sel) == 0){
#			return(stop("Invalid Anchor Protein Selection"))
#		}
#		eset <- eset[sel,]
#	}
	# No big difference if using all features. Avoids down stream bug when applying norm factors
	
	### make sure minFeaturesPerBin is not > tot. nb. features
	minFeaturesPerBin <- min(c(minFeaturesPerBin,nrow(eset)))
	
	# get all ratios to sample 1
	# @TODO How to select reference run?
	ratios <- log2(exprs(eset)) - log2(exprs(eset)[,1])
	#ratios <- log2(exprs(eset)) - log2(apply(exprs(eset),1,median,na.rm=T))
	
	### get median ratio per minute bin
	roundedRT <- round(fData(eset)$retentionTime)
	rtBin <- sort(unique(round(fData(eset)$retentionTime)))
	rtBin <-  rtBin[!is.na(rtBin)] 
	
	# normalization factors per retention time bin bin
	rtNormFactors <- data.frame(row.names=rtBin)
	
	# iterate over samples
	for(i in 1:ncol(ratios)){
		ratio <- ratios[,i]
		rtNFac <- data.frame(row.names=rtBin)
		
		# iterate over all retention time bins
		for(rt in rtBin){
			
			match <- roundedRT %in% rt
			
			### while number of features are < minFeaturesPerBin -> include neighbouring bins
			rtBins <- rt
			while(sum(match) < minFeaturesPerBin ){
				rtBins <- ((min(rtBins)-1):(max(rtBins)+1))
				match <- roundedRT %in% rtBins
			}
			# median or sum?
			rtNFac[as.character(rt),1] <- median(ratio[match],na.rm=T)
		}
		
		# store rt bin norm factors
		rtNormFactors <- cbind(rtNormFactors,rtNFac)
		
	}
	names(rtNormFactors) <- rownames(pData(eset))
	return(rtNormFactors)
}

#' Normalization data per retention time bin
#' @param eset ExpressionSet
#' @param rtNormFactors  obtained using getRTNormFactors
#' @return data.frame normalization factors per retention time bin (minute)
#' @export
#' @import limma Biobase
#' @note  No note
#' @details Normalize for variations in elelctrospray ionization current.
#' @seealso  \code{\link{getRTNormFactors}}
#' @references In Silico Instrumental Response Correction Improves Precision of Label-free Proteomics and Accuracy of Proteomics-based Predictive Models, Lyutvinskiy et al. (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23589346} 
#' @examples print("No examples")
rtNormalize <- function(eset,rtNormFactors){
	
	### check if eset contain necessary retentionTime colummn
	if(is.null(fData(eset)$retentionTime)){
		stop("retentionTime not found")
	}
		
	roundedRT <-  round(fData(eset)$retentionTime)
	
	# FIX
	if(sum(!(as.character(roundedRT)[!is.na(roundedRT)] %in% row.names(rtNormFactors))) > 0){
		stop("ERROR: Missing R.T. Norm Factors.")
	}
		
	for(i in 1:ncol(eset)){
		# normalize 
		
		nF <- rtNormFactors[as.character(roundedRT),i]
		nF[is.na(nF)] <- 0
		
		exprs(eset)[,i]  <- 2^(log2(exprs(eset)[,i]) - nF)
	}
	
	return(eset)
}


#' Get normalization factors. calculated as summed/median signal per run (column) over summed/median of first run. 
#' @param eset ExpressionSet
#' @param method c("sum","median)
#' @return vector of normalization factors
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @examples print("No examples")
getGlobalNormFactors <- function(eset, method="sum"){
	
	sel <- rep(T,nrow(eset))
	
	### check if isNormAnchor and isFiltered columns are defined, if -> get normalization factors from nonFiltered anchor proteins
	if(!is.null(fData(eset)$isNormAnchor) & !is.null(fData(eset)$isFiltered)){
		
		### only use feature qunatified in all runs for normalization
		isAllSignal <- as.vector(apply(is.finite(exprs(eset)),1,sum) == ncol(eset))
		
		sel <- fData(eset)$isNormAnchor & !fData(eset)$isFiltered & isAllSignal
	
		if(sum(sel) == 0){
			stop("Error: getGlobalNormFactors -> Invalid Anchor Protein Selection ")
		}
	}
	
	if(method == "sum"){
		rawDataIdx = apply(data.frame(exprs(eset)[sel,]),2, FUN=sum, na.rm=T)
	}else if(method == "median"){
		rawDataIdx = apply(data.frame(exprs(eset)[sel,]),2, FUN=median, na.rm=T)
	}else{
		stop("Error: Unknown Global Normalization Method", method)		
	}
	
	normFactors = as.numeric(rawDataIdx[1]) / as.numeric(rawDataIdx)
	
	return(normFactors)
	
}

#' Normalize, Norm factors calculated as median signal per run (column) over median of first run. 
#' @param eset ExpressionSet
#' @param globalNormFactors globalNormFactors
#' @return eset ExpressionSet
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @seealso getGlobalNormFactors
#' @examples print("No examples")
globalNormalize <- function(eset,globalNormFactors){
	
	normData <- data.frame(matrix(t(apply(exprs(eset), 1, function(.rawData)mapply(globalNormFactors, .rawData, FUN="*"))),ncol=ncol(exprs(eset)),dimnames=list(row.names(exprs(eset)),names(exprs(eset)))))
	names(normData) <- sampleNames(eset) 
	
	return(createExpressionDataset(as.matrix(normData),pData(eset),fData(eset)))
	
}


#' Normalize 
#' @param eset ExpressionSet
#' @param	method c("global","rt","quantile") 
#' @return eset ExpressionSet
#' @export
#' @import limma
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @seealso getGlobalNormFactors, getRTNormFactors
#' @examples print("No examples")
sqNormalize <- function(eset, method="global"){
	
	esetNorm <- eset
	
	if("global" %in% method){
		
		globalNormFactors <- getGlobalNormFactors(esetNorm)
		
		### add normalization factors to ExpressionSet
		pData(esetNorm)$globalNormFactors <- globalNormFactors
		esetNorm <- globalNormalize(esetNorm,globalNormFactors)
	}
	if("rt" %in% method){
			
			rtNormFactors <- getRTNormFactors(esetNorm, minFeaturesPerBin=100)
			esetNorm <- rtNormalize(eset,rtNormFactors)
			
	}
	if("quantile" %in% method){
		# limma
		exprs(eset) <- normalizeQuantiles(exprs(eset))	
	}
	
	### @ experimental
#	if("vsn" %in% method){
#		library(vsn)
#		esetNorm <- justvsn(eset)
#		exprs(esetNorm) <- 2^exprs(esetNorm)
#	}
	
#	if(identical(esetNorm,eset)){
#		warning("No normalization performed")
#		print( apply(exprs(esetNorm),2,median,na.rm=T))
#	}
	
	return(esetNorm)

}

#' Summarize replicate signal per condition (min)
#' @param eset ExpressionSet
#' @param method median (default), mean, max, min, sd
#' @return data.frame of per condition signals
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getSignalPerCondition <- function(eset,method="median"){
	
	#conditionNames <- levels(pData(eset)$condition)
	conditionNames <- unique(as.character(pData(eset)$condition))
	
	perCondSignal <- data.frame(row.names=rownames(eset))
	
	if(method=="median"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,median, na.rm=T))
			}else{
				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
			}
		}
	}else if(method=="mean"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,mean, na.rm=T))
			}else{
				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
			}
		}
	}
	else if(method=="max"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,max, na.rm=T))
			}else{
				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
			}
		}
	}else if(method=="min"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,min, na.rm=T))
			}else{
				perCondSignal <- cbind(perCondSignal,exprs(eset)[,condMatch])
			}
		}
	}else if(method=="sd"){
		for(cond in conditionNames){
			condMatch <-  cond== pData(eset)$condition
			if(sum(condMatch) > 1){
				perCondSignal <- cbind(perCondSignal,apply(exprs(eset)[,condMatch],1,sd, na.rm=T))
			}else{
				perCondSignal <- cbind(perCondSignal,NA)
			}
		}
	}
	else{
		stop("Unknown method: method ")
	}
	names(perCondSignal) <- conditionNames
	
	return(perCondSignal)
}



#' Calculate Mean of X most intense features
#' @param entryData data.frame listing feature intensities of one entry. Typically rows corresponds to Peptide entries of one protein
#' @param topX best X flyers
#' @return vector of topX intensities per column (sample)
#' @export
#' @note  No note
#' @details No details
#' @references Absolute quantification of proteins by LCMSE: A virtue of parallel MS acquisition, Silva (2006), \url{http://www.ncbi.nlm.nih.gov/pubmed/16219938},  Critical assessment of proteome-wide label-free absolute abundance estimation strategies. Ahrne (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23794183}
#' @examples print("No examples")
getTopX <- function(entryData,topX=3){
	
	# if number of rows are fewer than topX  
	topX <- min(c(topX,nrow(entryData)))
	
	### get row order based on decreasing row sum
	if( !is.null(ncol(entryData))  && (ncol(entryData) > 1)  ){
		o <- order(apply(entryData,1,sum,na.rm=T),decreasing=T)[1:topX]
		
		if(topX==1){
			
			return(entryData[o,])
		}
		return(apply(entryData[o,],2,mean,na.rm=T))
	}
	
	# only one column
	o <- order(entryData,decreasing=T)[1:topX]
	return(mean(entryData[o],na.rm=T))	
	
}

#' Calculate intensity-based absolute-protein-quantification (iBAQ) metric per protein 
#' @param eset protein level ExpressionSet
#' @param proteinDB list protein sequneces
#' @param peptideLength peptide length interval (to get number of peptides used for normalization)
#' @param nbMiscleavages number of mis-cleavages allowed when digesting protein sequneces in silico (to get number of peptides used for normalization)
#' @param proteaseRegExp protease Reg Exp cleavage rule
#' @return ExpressionSet
#' @export
#' @note  No note
#' @details No details
#' @references Global quantification of mammalian gene expression control, Schwanhausser (2011), \url{http://www.ncbi.nlm.nih.gov/pubmed/21593866}, 
#' Critical assessment of proteome-wide label-free absolute abundance estimation strategies. Ahrne (2013), \url{http://www.ncbi.nlm.nih.gov/pubmed/23794183}
#' @examples print("No examples")
getIBAQEset <- function(eset
		, proteinDB=NA
		, peptideLength = c(5,36)
		, nbMiscleavages = 0
		, proteaseRegExp=.getProteaseRegExp("trypsin")
){
	
	# get number of detectable peptides per protein
	nbPeptides <- vector(length=nrow(eset))
	i <- 0
	for(i in 1:nrow(eset)){
		
		ac <- as.character(fData(eset)$proteinName[i])
		
		### keep first listed entry sp|P04049|RAF1_HUMAN;sp|P15056|BRAF_HUMAN  -> sp|P04049|RAF1_HUMAN
		#ac <- gsub("\\;.*","",ac)
		
		nbPep <- NA
		if( !is.null(proteinDB[[ac]]) ){
			nbPep <- getNbDetectablePeptides(getPeptides(proteinDB[[ac]],proteaseRegExp=proteaseRegExp,nbMiscleavages=nbMiscleavages),peptideLength=peptideLength)
		}else{
			warning("WARN: ",ac," NOT FOUND IN PROTEIN DATABASE")
			#cat("WARN: ",ac," NOT FOUND IN PROTEIN DATABASE\n")
		}
		nbPeptides[i] <- nbPep
	}
	
	### create new "absquant" eset
	esetIBAQ <- eset
	# scale protein intensity by number of detectable peptides 
	exprs(esetIBAQ) <- exprs(eset) / nbPeptides
	
	### store normalization factors
	fData(esetIBAQ)$nbTrypticPeptides <- nbPeptides
	
	return(esetIBAQ)
	
}

### require colnames "signal", "cpc"
#' Leave-One-Out Cross Validate Qunatification Model
#' @param df data.frame  of two columns 1) "signal" - ms metric 2) "cpc" absolute quantity
#' @return data.frame of fold errors per (left-out) protein
#' @export
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @seealso NA
#' @examples print("No examples")
getLoocvFoldError <- function(df){
	
	ok <- is.finite(df[,1]) & is.finite(df[,2])
	df <- df[ok,]
	
	foldError <- vector(length=nrow(df))
	
	for (i in 1:nrow(df)) {
		fit <- lm(cpc ~ signal, data=df[-i,] )
		foldError[i] <- 10^predict(fit,newdata=df[i,]) / 10^df[i,]$cpc
		
	}
	
	### if foldErro < 1, take neg. of inverse	
	foldError[foldError < 1] <- -1/foldError[foldError < 1]
	
	### return fold error per protein
	foldError <- data.frame(foldError,row.names=rownames(df))
	
	return(foldError)
	
}

#' Set value to NA if it deviatves with more than 1.5 * IQR from lower/upper quantile
#' @param x vector numeric
#' @param na.rm logical indicating whether missing values should be removed.
#' @param ... qunatile args 
#' @export
#' @importFrom stats quantile IQR
#' @note  No note
#' @details No details
#' @keywords normalization
#' @references NA 
#' @seealso NA
#' @examples print("No examples")
removeOutliers <- function(x, na.rm = TRUE, ...){
	qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
	H <- 1.5 * IQR(x, na.rm = na.rm)
	y <- x
	y[x < (qnt[1] - H)] <- NA
	y[x > (qnt[2] + H)] <- NA
	
	return(y)
}


#' Roll up feature intensites per unique colum combination
#' @param eset ExpressionSet
#' @param featureDataColumnName vector of column names e.g. peptide or proteinName
#' @param method "sum", "mean" or "top3" 
#' @return ExpressionSet object
#' @export
#' @details featureDataColumnName = c("peptide","charge","ptm"), method= c("sum"), sums up intensities per peptie modification charge state
#' @import Biobase data.table
#' @note  No note
#' @references No references
#' @examples print("No examples")
rollUp <- function(eset, method = "sum", 	featureDataColumnName =  c("proteinName") ){
	
	# HACK to please CRAN CHECK "rollUp: no visible binding for global variable "idScore""
	idx <- idScore <- V1 <- allAccessions <- NULL

	### apply filter
	eset <- eset[!fData(eset)$isFiltered,] 
	
	selectedColumns <- names(fData(eset)) %in% featureDataColumnName
	allIndexTags <- as.vector(unlist(apply(data.frame(fData(eset)[,selectedColumns]),1,function(t){
								return(paste(as.vector(unlist(t)),collapse="_"))
							})))
	
	DT <- data.table(
			idx = allIndexTags
			,exprs(eset)
			,fData(eset)	
	)
	setkey(DT,idx)		
	
	if(method =="sum"){
		rDT <- DT[, lapply(.SD, sum, na.rm=TRUE), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
		rolledAssayData <- data.frame(rDT,row.names=rDT[,1], check.names=F)
	}else if(method =="median"){
		rDT <- DT[, lapply(.SD, median, na.rm=TRUE), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
		rolledAssayData <- data.frame(rDT,row.names=rDT[,1], check.names=F)
		
	}else if(method =="mean"){
		rDT <- DT[, lapply(.SD, mean, na.rm=TRUE), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
		rolledAssayData <- data.frame(rDT,row.names=rDT[,1], check.names=F)
	}else if(method =="top3"){
		rDT <- DT[, lapply(.SD, getTopX, topX=3), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
		rolledAssayData <- data.frame(rDT,row.names=rDT[,1], check.names=F)
	}else if(method =="top1"){
		rDT <- DT[, lapply(.SD, getTopX, topX=1), by=idx, .SDcols=c(2:(ncol(eset)+1)) ]
		rolledAssayData <- data.frame(rDT,row.names=rDT[,1], check.names=F)
	}
	
	# idScore
	if("idScore" %in% names(DT) ){
		
		### get fData of highest scoring row per rollUP level
		indices <- DT[, .I[getMaxIndex(idScore) ], by=idx]$V1
		rolledFData <- data.frame(data.frame(DT)[indices,names(fData(eset))],row.names=rownames(rolledAssayData))
		# replace idScore colum, by summed score
		sumDT <- DT[, sum(idScore,na.rm=T), by=idx ]
		rolledFData$idScore = sumDT[,V1]
	}else{
		
		# @TODO inefficient solution -> speed up	
		# allColumns but idx,intensity data and idScore
		#selColOld <- which(!(names(DT)  %in%  c(colnames(exprs(eset)),"idx","idScore")))
		selCol <- which((names(DT)  %in% names(fData(eset))))
		
		#rolledFDataOld <- data.frame(DT[, lapply(.SD, .getFirstEntry), by=idx, .SDcols=selCol],row.names=rownames(rolledAssayData))[,2:(length(selCol)+1)]
		rolledFData <- data.frame(DT[, lapply(.SD, .getFirstEntry), by=idx, .SDcols=selCol],row.names=rownames(rolledAssayData))[,names(fData(eset))]
		
#		print(names(rolledFData))
#		print("")
#		print(names(rolledFDataOld))	
	}
	
	### concatenate allAccessions
	if("allAccessions" %in% names(DT)){
		rolledFData$allAccession <- DT[, list( allAccessionsTMP  = paste(unique(unlist(strsplit(paste(allAccessions,collapse=";"),";"))),collapse=";") ), by = key(DT)]$allAccessionsTMP	
	}
		
	rolledAssayData <- as.matrix(rolledAssayData)
	rolledAssayData[rolledAssayData == 0 ] <- NA
	#names(rolledFData) <- names(fData(eset))
	
	if(!is.null(rolledFData$isNormAnchor)){
		rolledFData$isNormAnchor <- T
	}
	# reset filter
	if(!is.null(rolledFData$isFiltered)){
		rolledFData$isFiltered <- F
	}

	### set peptides per protein
	rolledFData$nbPeptides <- getNbPeptidesPerProtein(eset)[as.character(rolledFData$proteinName)]
	
	### if 


	return(createExpressionDataset(expressionMatrix=rolledAssayData,expDesign=pData(eset),featureAnnotations=rolledFData))
}

#' Per Feature Normalization
#' @param eset ExpressionSet
#' @param normFactors matrix normalization factors (logged) (row names are proteins)
#' @return ExpressionSet object
#' @export
#' @details Example Usage: Normalize phospho peptide signals for Protein Changes 
#' @note  No note
#' @references No references
#' @examples print("No examples")
perFeatureNormalization <- function(eset,normFactors){
	
	coveredPeptideSel <- fData(eset)$proteinName %in% rownames(normFactors)
	
	if(sum(coveredPeptideSel) == 0){
		warning("No shared entries between target data and normalisation data")
		return(eset)
	}
	
	# make sure target data and norm data have same dimensions
	if( !all(colnames(normFactors) %in% pData(eset)$condition) ){
		stop("Invalid norm factors")
	}
	
	# normalise log ratios
	exprs(eset)[coveredPeptideSel,]	<- exprs(eset)[coveredPeptideSel, ] - normFactors[as.character(fData(eset)[coveredPeptideSel,]$proteinName),pData(eset)$condition]

	
	return(eset)	
	
}

#' Standardise data
#' @param d vector or data.frame or matrix
#' @return vector or data.frame or matrix
#' @export
#' @note  No note
#' @details No details
#' @examples print("No examples")
standardise <- function(d){
	
	d[!is.finite(as.matrix(d))] <- NA
	
	if((class(d) == "data.frame") | (class(d) == "matrix")){
		for(i in 1:ncol(d)){
			d[,i] <- d[,i]- mean(d[,i],na.rm=T)
			d[,i] <- d[,i]/ sd(d[,i],na.rm=T)
		}
		return(d)
	}else{
		return( (d - mean(d,na.rm=T)) / sd(d,na.rm=T) )
	}
}



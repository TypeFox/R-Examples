# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

### @TODO 
parseCSV <- function(file=file,expDesign=expDesign){}


################################## LFQ ##################################

###  get column indices of intensity data
#.getProgenesisCsvIntColIndices <- function(file){
#	
#	con <- file(file) 
#	open(con);
#	line1 <- readLines(con, n = 1)
#	close(con)
#	intStartCol <- grep("Normalized abundance",unlist(strsplit(line1,",")))
#	intEndCol <- grep("Raw abundance",unlist(strsplit(line1,",")))-1
#	nbSamples <- intEndCol - intStartCol + 1
#	return((intStartCol:intEndCol)+nbSamples)
#}

##  get column indices of spectral count data
#' @export
.getProgenesisCsvExpressionColIndices <- function(file, method="auc"){
	
	con <- file(file) 
	open(con);
	line1 <- readLines(con, n = 1)
	close(con)
	
	if(grepl("Raw abundance",line1)){
		### non-fractionated data
		intStartCol <- grep("Normalized abundance",unlist(strsplit(line1,",")))
		intEndCol <- grep("Raw abundance",unlist(strsplit(line1,",")))-1
		nbSamples <- intEndCol - intStartCol + 1
	}else{
		### fractionated data
		intStartCol <- grep("Normalized abundance",unlist(strsplit(line1,",")))
		intEndCol <- grep("ectral counts",unlist(strsplit(line1,",")))-1
		nbSamples <- intEndCol - intStartCol + 1
		
		intStartCol <- intStartCol - nbSamples
		intEndCol <- intEndCol - nbSamples
	}
	
	if(method =="auc"){
		return(return((intStartCol:intEndCol)+nbSamples))
	}else if(method =="spc"){
		
		startCol <- grep("ectral counts",unlist(strsplit(line1,",")))
		endCol <- startCol+nbSamples-1
		return((startCol:endCol))
		
	}else{
		stop("Unknown quantification indices\n")
	}
	
	#@TODO fractionated data. Exports only contain Normalized Abundance data
	
}


# Get Experimental 	
#' Parse Experimental Design from Progenesis Csv Export
#' @param file path to progenesis csv file
#' @param expressionColIndices default .getProgenesisCsvExpressionColIndices(file)
#' @return data.frame describing experimental design
#' @export
#' @importFrom utils read.csv
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
getExpDesignProgenesisCsv <- function(file, expressionColIndices = .getProgenesisCsvExpressionColIndices(file)){
	
	header <-  read.csv(file,nrows=2, check.names=F)[,expressionColIndices]
	header <- data.frame(header) ## in case only one cond one sample
	
	conditionLine <- as.character(unlist(header[1,]))
	currentCond <- as.character(conditionLine[1])
	condition <- c()
	
	for(i in 1:length(conditionLine)){
		
		if(nchar(conditionLine[i]) > 0 ){
			currentCond <- conditionLine[i]
		}
		condition <- c(condition,currentCond) 
		
	}
	
	## set control condition
	isControl <- rep(F,length(condition))
	isControl[condition[1] == condition] <- T
	
	### check if sample names are unique, if not, -> add running number 
	# substitute - -> .
	#sampleNames <- gsub("\\-",".",as.character(as.vector(unlist(header[2,]))))
	sampleNames <- as.character(as.vector(unlist(header[2,])))
	if(T %in% (table(sampleNames) > 1)) sampleNames <- paste(sampleNames,1:length(sampleNames),sep="_") 
	
	expDesign <- data.frame(
			condition=condition
			,isControl=isControl
			,row.names=sampleNames
	)
	return(expDesign)
}

#' Parse Progenesis Protein Csv
#' @param file path to Progenesis Protein csv file
#' @param expDesign experimental design data.frame
#' @param method auc (area under curve) or spc (spectral count)
#' @return ExpressionSet object
#' @export
#' @importFrom utils read.csv
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
parseProgenesisProteinCsv <- function(file=file,expDesign=expDesign, method="auc"){
	
	res <- read.csv(file,skip=2,allowEscapes=T, check.names=F)
	expMatrix <- as.matrix(res[,.getProgenesisCsvExpressionColIndices(file, method=method)])
	
	# set 0 features to NA
	expMatrix[expMatrix == 0] <- NA
	# discard features where all intensities are NA (i.e. 0)
	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))
	
	row.names(expMatrix) <- res$Accession
	
#	[1] "Accession"                      "Peptide count"                 
#	[3] "Peptides used for quantitation" "Confidence score"              
#	[5] "Anova (p)*"                     "Max fold change"               
#	[7] "Highest mean condition"         "Lowest mean condition"         
#	[9] "Description"                    "A11-03216"                     
	
	featureAnnotations <- data.frame(
			proteinName=res$Accession
			,proteinDescription=res$Description
			,idScore=res[,"Confidence score"]
			#	,nbPeptides=res[,"Peptides used for quantitation"] ### old Progenesis 
			,nbPeptides=res[,which(names(res) == "Unique peptides" | names(res) == "Peptides used for quantitation")[1]] # Progensis QI
			,isNormAnchor=rep(T,nrow(expMatrix))
			,isFiltered=rep(F,nrow(expMatrix))
			,row.names=res$Accession)
	
	featureAnnotations <- featureAnnotations[!allColNA,]
	
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
	
	### re-order and exclude channels 
	expMatrix <- as.matrix(expMatrix[!allColNA ,rownames(expDesign)])
	colnames(expMatrix) <- rownames(expDesign)
	
	return(createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations))
}

#' Parse Progenesis Feature Csv Export
#' @param file path to Progenesis Feature csv file
#' @param expDesign experimental design data.frame
#' @param method auc (area under curve) or spc (spectral count)
#' @return ExpressionSet object
#' @export
#' @import Biobase
#' @importFrom utils read.csv
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
parseProgenesisFeatureCsv <- function(file=file,expDesign=getExpDesignProgenesisCsv(file), method="auc"){
	
	### stop if not all samples labelled with a given condition are assigned as control
	# Example:
	#			condition isControl
	#	A11.09066 Condition 1      TRUE
	#	A11.09067 Condition 1     FALSE
	#	A11.09068 Condition 1     FALSE
	#	A11.09070 Condition 2     FALSE
	#	A11.09071 Condition 2     FALSE
	#	A11.09072 Condition 2     FALSE
	if(length(unique(expDesign$condition))  == length(unique(expDesign[!expDesign$isControl,]$condition))){
		stop("Invalid Exp. Design")
	}
	
	# read csv file
	res <- read.csv(file,skip=2,allowEscapes=T,check.names=F)
	expMatrix <- as.matrix(res[,.getProgenesisCsvExpressionColIndices(file, method=method)])
	
	# set 0 features to NA
	expMatrix[expMatrix == 0] <- NA
	# discard features where all intensities are NA (i.e. 0)
	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))
	
	### Mass filed missing in old Progenesis Feature export
	if(is.null(res$Mass)){res$Mass <- rep(NA,nrow(expMatrix))}
	
	# added
	ptm <- res[,"Variable modifications ([position] description)"]
	# strip ; character, which appears sometimes (Proteome Discoverer data?) example: modifAnnot: [5] Phospho; (S)|[7] Phospho (S)
	ptm <- gsub("\\;","",ptm)
	
	nbPtmsPerPeptide <- unlist(lapply(ptm,function(t){
						t <- as.character(t)
						return(sum(unlist(gregexpr("\\|",t)[[1]]) > 0) + (nchar(t) > 0)  )}))
	
	
#	"#"                                              
#	[2] "m/z"                                            
#	[3] "Retention time (min)"                           
#	[4] "Retention time window (min)"                    
#	[5] "Mass"                                           
#	[6] "Charge"                                         
#	[7] "Max fold change"                                
#	[8] "Highest mean condition"                         
#	[9] "Lowest mean condition"                          
#	[10] "Anova"                                          
#	[11] "Maximum CV"                                     
#   [60] "Notes"                                          
#	[61] "Score"                                          
#	[62] "Mass error (u)"                                 
#	[63] "Mass error (ppm)"                               
#	[64] "Protein"                                        
#	[65] "Sequence"                                       
#	[66] "Variable modifications ([position] description)"
#	[67] "Description"  
	
	featureAnnotations <- data.frame(
			proteinName=res$Protein
			,proteinDescription=res$Description
			,peptide=res$Sequence
			,idScore=as.numeric(as.character(res$Score))
			,mass=res$Mass
			,pMassError=res[,"Mass error (ppm)"]
			,mz=res[,"m/z"]
			,retentionTime=res[,"Retention time (min)"]
			,charge=as.numeric(res$Charge)
			,ptm=ptm
			,isNormAnchor=rep(T,nrow(expMatrix))
			,isFiltered=rep(F,nrow(expMatrix))
			#		,row.names=res$Accession
			# added
			,nbPtmsPerPeptide = nbPtmsPerPeptide
	)
	
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
	## @TODO what if "X001_Yewtivya" "001_Yewtivya"
	#grepl("X[0-9]",colnames(expMatrix))
	
	### re-order and exclude channels  
	#expMatrix <- expMatrix[,rownames(expDesign)]
	
	#return(createExpressionDataset(expressionMatrix=expMatrix[!allColNA,],expDesign=expDesign,featureAnnotations=featureAnnotations[!allColNA,]))
	
	#@TODO 
	#allColNA <- rep(F,length(allColNA))
	


	# discard non peptide annotated rows
	isPep <- !is.na(featureAnnotations$idScore) 
	
	
#	print(head(featureAnnotations))
#	print(nrow(featureAnnotations))
#	print(length(allColNA))
#	print(length(isPep))
#	print(head(expMatrix))
#	print(nrow(expMatrix))
#	print(rownames(expDesign))
	
	featureAnnotations <- data.frame(featureAnnotations)[!allColNA & isPep,]
		
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
		
	### re-order and exclude channels  
	if(ncol(expMatrix) > 1){
		expMatrix <- as.matrix(expMatrix[!allColNA & isPep ,rownames(expDesign)])
	}else{ # to avoid crash when only one run
		expMatrix <- as.matrix(expMatrix[!allColNA & isPep,])
	}
	
	colnames(expMatrix) <- rownames(expDesign)
	
	return(createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations))
	
}



#' Parse Progenesis Peptide Measurement Csv Export
#' @param file path to Progenesis Peptide Measurement csv file
#' @param expDesign experimental design data.frame
#' @param method auc (area under curve) or spc (spectral count)
#' @param expressionColIndices default .getProgenesisCsvExpressionColIndices()
#' @return ExpressionSet object
#' @export
#' @import Biobase
#' @importFrom utils read.csv
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
parseProgenesisPeptideMeasurementCsv <- function(file,expDesign=expDesign,	method="auc", 	
		expressionColIndices = .getProgenesisCsvExpressionColIndices(file, method=method) ){
	
	# HACK to please CRAN CHECK "rollUp: no visible binding for global variable "Sequence""
	Sequence <- Score <- resDTIndex <- proteinScore <- Accession <- NULL
	
	### stop if not all samples labelled with a given condition are assigned as control
# Example:
#	condition isControl
#	A11.09066 Condition 1      TRUE
#	A11.09067 Condition 1     FALSE
#	A11.09068 Condition 1     FALSE
#	A11.09070 Condition 2     FALSE
#	A11.09071 Condition 2     FALSE
#	A11.09072 Condition 2     FALSE
	if(length(unique(expDesign$condition))  == length(unique(expDesign[!expDesign$isControl,]$condition))){
		stop("Invalid Exp. Design")
	}
	
	# read csv file
	res <- read.csv(file,skip=2,allowEscapes=T,check.names=F)
	#	[1] "#"                                     
	#	[2] "Retention time (min)"                  
	#	[3] "Charge"                                
	#	[4] "m/z"                                   
	#	[5] "Measured mass"                         
	#	[6] "Mass error (u)"                        
	#	[7] "Mass error (ppm)"                      
	#	[8] "Score"                                 
	#	[9] "Sequence"                              
	#	[10] "Modifications"                         
	#	[11] "Accession"                             
	#	[12] "All accessions (for this sequence)"    
	#	[13] "Grouped accessions (for this sequence)"
	#	[14] "Shared accessions (for this sequence)" 
	#	[15] "Description"                           
	#	[16] "Use in quantitation"                   
	#	[17] "A15-02164"                             
	#	[18] "A15-02164"                             
	#	[19] "A15-02164"                             
	#	[20] "Mass"    
	
	################################## PROTEIN INFERENCE ################################## 
	
	# OCCAMS RAZOR IMPLEMENTATION (HIGHEST SCORING PROTEIN TAKES IT ALL)
	# NOT COMPATIBLE WITH STANDARD iBAQ AND TOP3 ABSOLUT QUANTIFICATION! # EXAMPLE: Peptide Measurement Output
	# Accession	sp|O43707|ACTN4_HUMAN
	# All accessions (for this sequence)	sp|O43707|ACTN4_HUMAN;sp|P12814|ACTN1_HUMAN;sp|P35609|ACTN2_HUMAN;sp|Q08043|ACTN3_HUMAN
	# Grouped accessions (for this sequence)	sp|P35609|ACTN2_HUMAN;sp|Q08043|ACTN3_HUMAN
	# Shared accessions (for this sequence)	sp|P12814|ACTN1_HUMAN
		
		### 
		
	#	Feature	Score	Sequence	Accession				All accessions (for this sequence)	Grouped accessions (for this sequence)	Shared accessions (for this sequence)
	#	1770	62.91	AALSALESFLK	sp|P78527|PRKDC_HUMAN	sp|P78527|PRKDC_HUMAN		
		
	# A "Grouped accessions (for this sequence)" never appears in the "Accession" column as ProteinScore(GroupedProtein) <= ProteinScore(Accession)
	# allGrouped <-  as.vector(unlist(lapply(resDT$"Grouped accessions (for this sequence)",function(t){strsplit(as.character(t),";")})))
	# sum(allGrouped %in% as.character(resDT$Accession))
		
	# OCCAMS RAZOR IMPLEMENTATION (HIGHEST SCORING PROTEIN/PROTEIN GROUP TAKES IT ALL)
	# NOT COMPATIBLE WITH STANDARD iBAQ AND TOP3 ABSOLUT QUANTIFICATION! 
	
	
	cat("INFO: ASSIGNING PEPTIDES TO PROTEINS APPLYING OCCAM'S RAZOR \n" )	
	
	#Check that file includes column "Grouped accessions (for this sequence)"
	if(!("Grouped accessions (for this sequence)" %in% names(res))) stop("Column \"Grouped accessions (for this sequence)\" missing in file ",file)
	names(res)[1] <- "Feature"
	res$Score <- as.numeric(as.character(res$Score))
	
	resDT <-  data.table(res,key="Accession")
	resDT$"Grouped accessions (for this sequence)" <- as.character(resDT$"Grouped accessions (for this sequence)")
	resDT$Accession <- as.character(resDT$Accession)
	
	# 1) FIND PROTEIN GROUPS.  A PROTEIN GROUP IS A SET OF PROTEINS THAT A) SHARE THE EXACT SAME PEPTIDES.
	# B) All proteins of the group ONLY map to this set of peptides
	# In the Progensis Export "Peptide Measurments", the column "Grouped accessions (for this sequence)"
	# is defined according to Condition A) only.
	# "leading" refers to column "Accession" 
	isGroup <- nchar(resDT$"Grouped accessions (for this sequence)") > 0
	groupedAccessions <- unique(resDT$Accession[isGroup])
	nonGroupedAccessions <- unique(resDT$Accession[!isGroup])
	onlyGrouped <- groupedAccessions[ !(groupedAccessions %in% nonGroupedAccessions) ]
	isOnlyGrouped <- resDT$Accession %in% onlyGrouped
	
	leadingGroupAcTable <- table(unique(data.frame(resDT$Accession[isOnlyGrouped],resDT$"Grouped accessions (for this sequence)"[isOnlyGrouped]))[,1])
# leading AC only part of one group
	leadingACTrueGroup <- names(leadingGroupAcTable[leadingGroupAcTable == 1])
	
# add column myProteinGroup
	resDT <- cbind(resDT,myProteinGroup= resDT$Accession )
	isTrueGroup <- resDT$Accession %in% leadingACTrueGroup
	resDT$myProteinGroup[isTrueGroup] <- paste(resDT[isTrueGroup,]$Accession,resDT[isTrueGroup,]$"Grouped accessions (for this sequence)",sep=";")
	
# 2) ASSIGN BEST SCORING PEPTIDE TO EACH FEATURE
# Note that a peptide can be listed multiple times per 
# Feature, but assigned different ACs (i.e. the peptide is shared)
# A) Roll-up on feature level, keeping the best scoring peptide
# B) list all featureNB_peptides to be kept
# C) filter resDT keeping only these featureNB_peptides
	setkey(resDT,"Feature")
	
# A) Roll-up on feature level, keeping the best scoring peptide
	featureDTTmp <- resDT[, list( peptide  = Sequence[order(Score,decreasing=T)[1]]), by = key(resDT)]
# B) list all featureNB_peptides to be kept
	keptFeaturePeptides <- unique(paste(featureDTTmp$Feature,as.character(featureDTTmp$peptide),sep=""))
# C) filter resDT keeping only these featureNB_peptides
	resDT <- resDT[paste(resDT$Feature,as.character(resDT$Sequence),sep="") %in% keptFeaturePeptides,]
	
# 3)  Score all proteins
	setkey(resDT,"Accession")
	proteinDT <- resDT[, list(proteinScore = sum(Score,na.rm=T)), by = key(resDT)]
	rownames(proteinDT) <- proteinDT$Accession
	
# 4) Add proteinScore to resDT 
#resDT <- cbind(resDT,proteinScore=proteinDT[resDT$Accession,]$proteinScore, tmpAC = proteinDT[resDT$Accession,]$Accession) # check
	resDT <- cbind(resDT,proteinScore=proteinDT[resDT$Accession,]$proteinScore)
	
# 5) Roll-up on feature level, keeping OCCAMS myProteinGroup per Feature
	setkey(resDT,"Feature")
	resDT <- cbind(resDT,resDTIndex=1:nrow(resDT))
	featureDT <- resDT[resDT[, list(  selectedResDTIndex  = resDTIndex[order(proteinScore,decreasing=T)[1]]), by = key(resDT)]$selectedResDTIndex,]
	
	#@TODO add alternative/subset accessions 
	# 6) Roll-up on peptide level concatenating all proteins matching a given peptide 
	# create peptide to protein dict
	setkey(resDT,"Sequence")
	peptideDT <- resDT[, list(allAccessions = paste(unique(Accession),collapse=";")), by = key(resDT)] # could also include groups..
	rownames(peptideDT) <- peptideDT$Sequence 
	# get allProteins of a given peptide
	featureDT <- cbind(featureDT,peptideDT[as.character(featureDT$Sequence),]) 
		
	################################## PROTEIN INFERENCE END ################################## 
	
	expMatrix <- as.matrix(featureDT[,expressionColIndices,with=FALSE])
	
	# set 0 features to NA
	expMatrix[expMatrix == 0] <- NA
	# discard features where all intensities are NA (i.e. 0)
	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))
	
	### Mass filed missing in old Progenesis Feature export
	if(is.null(featureDT$Mass)){featureDT$Mass <- rep(NA,nrow(expMatrix))}
	
	# added
	ptm <- featureDT$Modifications
	# strip ; character, which appears sometimes (Proteome Discoverer data?) example: modifAnnot: [5] Phospho; (S)|[7] Phospho (S)
	ptm <- gsub("\\;","",ptm)
	
	nbPtmsPerPeptide <- unlist(lapply(as.character(ptm),function(t){
						return(sum(unlist(gregexpr("\\|",t)[[1]]) > 0) + (nchar(t) > 0)  )}))
	
	suppressWarnings(score <- as.numeric(as.character(featureDT$Score)))
	### non scored entries are assigned a score of zero
	score[is.na(score)] <- 0
	
	featureAnnotations <- data.frame(
			proteinName=featureDT$myProteinGroup
			,proteinDescription=featureDT$Description
			,peptide=featureDT$Sequence
			,idScore= score
			,mass=featureDT$"Measured mass"
			,pMassError=featureDT$"Mass error (ppm)"
			,mz=featureDT$"m/z"
			,retentionTime=featureDT$"Retention time (min)"
			,charge=as.numeric(featureDT$Charge)
			,ptm=ptm
			,isNormAnchor=rep(T,nrow(expMatrix))
			,isFiltered=rep(F,nrow(expMatrix))
			#		,row.names=res$Accession
			# added
			,nbPtmsPerPeptide = nbPtmsPerPeptide
			,allAccessions = featureDT$allAccessions 	#@TODO add alternative/subset accessions 
	)
	
	# discard non peptide annotated rows
	isPep <- score > 0  
	featureAnnotations <- data.frame(featureAnnotations)[!allColNA & isPep,]
	
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
	
	### re-order and exclude channels  
	if(ncol(expMatrix) > 1){
		expMatrix <- as.matrix(expMatrix[!allColNA & isPep,rownames(expDesign)])
	}else{ # to avoid crash when only one run
		expMatrix <- as.matrix(expMatrix[!allColNA & isPep,])
	}
	
	colnames(expMatrix) <- rownames(expDesign)
	
	return(createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations))
}

#parseProgenesisPeptideMeasurementCsv2 <- function(file,expDesign=expDesign,	method="auc" ){
#	
#	### stop if not all samples labelled with a given condition are assigned as control
#	# Example:
#	#	condition isControl
#	#	A11.09066 Condition 1      TRUE
#	#	A11.09067 Condition 1     FALSE
#	#	A11.09068 Condition 1     FALSE
#	#	A11.09070 Condition 2     FALSE
#	#	A11.09071 Condition 2     FALSE
#	#	A11.09072 Condition 2     FALSE
#	if(length(unique(expDesign$condition))  == length(unique(expDesign[!expDesign$isControl,]$condition))){
#		stop("Invalid Exp. Design")
#	}
#	
#	# read csv file
#	res <- read.csv(file,skip=2,allowEscapes=T,check.names=F)
#	#	[1] "#"                                     
##	[2] "Retention time (min)"                  
##	[3] "Charge"                                
##	[4] "m/z"                                   
##	[5] "Measured mass"                         
##	[6] "Mass error (u)"                        
##	[7] "Mass error (ppm)"                      
##	[8] "Score"                                 
##	[9] "Sequence"                              
##	[10] "Modifications"                         
##	[11] "Accession"                             
##	[12] "All accessions (for this sequence)"    
##	[13] "Grouped accessions (for this sequence)"
##	[14] "Shared accessions (for this sequence)" 
##	[15] "Description"                           
##	[16] "Use in quantitation"                   
##	[17] "A15-02164"                             
##	[18] "A15-02164"                             
##	[19] "A15-02164"                             
##	[20] "Mass"    
#	
#	expMatrix <- as.matrix(res[,.getProgenesisCsvExpressionColIndices(file, method=method)])
#	
#	# set 0 features to NA
#	expMatrix[expMatrix == 0] <- NA
#	# discard features where all intensities are NA (i.e. 0)
#	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))
#	
#	### Mass filed missing in old Progenesis Feature export
#	if(is.null(res$Mass)){res$Mass <- rep(NA,nrow(expMatrix))}
#	
#	# added
#	ptm <- res[,"Modifications"]
#	nbPtmsPerPeptide <- unlist(lapply(ptm,function(t){
#						t <- as.character(t)
#						return(sum(unlist(gregexpr("\\|",t)[[1]]) > 0) + (nchar(t) > 0)  )}))
#	
#	if("All accessions (for this sequence)"  %in% names(res)){
#		protein <- res[,"All accessions (for this sequence)"]
#	}else{
#		cat("WARN: All accessions per Peptide were not exported \n")
#		protein <- res[,"Accession"]
#	}
#	
#	suppressWarnings(score <- as.numeric(as.character(res$Score)))
#	### non scored entries are assigned a score of zero
#	score[is.na(score)] <- 0
#	
#	# sort protein accessions
##	protein <- as.character(protein)
##	protein <- as.vector(unlist(lapply(protein, function(prot){
##								
##								prot <- paste(sort(as.vector(unlist(strsplit(prot,";")))),collapse=";")
##								return(prot)
##								
##							}  )))
#	
#	featureAnnotations <- data.frame(
#			proteinName=protein
#			,proteinDescription=res$Description
#			,peptide=res$Sequence
#			,idScore= score
#			,mass=res[,"Measured mass"]
#			,pMassError=res[,"Mass error (ppm)"]
#			,mz=res[,"m/z"]
#			,retentionTime=res[,"Retention time (min)"]
#			,charge=as.numeric(res$Charge)
#			,ptm=ptm
#			,isNormAnchor=rep(T,nrow(expMatrix))
#			,isFiltered=rep(F,nrow(expMatrix))
#			#		,row.names=res$Accession
#			# added
#			,nbPtmsPerPeptide = nbPtmsPerPeptide
#	)
#	
#	### strip off added .1  A11.03216.1 -> A11.03216
#	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
#	## @TODO what if "X001_Yewtivya" "001_Yewtivya"
#	#grepl("X[0-9]",colnames(expMatrix))
#	
#	### re-order and exclude channels  
#	#expMatrix <- expMatrix[,rownames(expDesign)]
#	
#	#return(createExpressionDataset(expressionMatrix=expMatrix[!allColNA,],expDesign=expDesign,featureAnnotations=featureAnnotations[!allColNA,]))
#	
#	# discard non peptide annotated rows
#	isPep <- nchar(as.character(featureAnnotations$peptide)) > 0 
#	featureAnnotations <- data.frame(featureAnnotations)[!allColNA & isPep,]
#	
#	### strip off added .1  A11.03216.1 -> A11.03216
#	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))
#	
#	### re-order and exclude channels  
#	if(ncol(expMatrix) > 1){
#		expMatrix <- as.matrix(expMatrix[!allColNA & isPep,rownames(expDesign)])
#	}else{ # to avoid crash when only one run
#		expMatrix <- as.matrix(expMatrix[!allColNA & isPep,])
#	}
#	
#	colnames(expMatrix) <- rownames(expDesign)
#	
#	return(createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations))
#}


################################## LFQ END ##############################

################################## TMT ##################################

### get line number from which Scaffold file should be read
#' @export
.getSkipLineNb <- function(fileName){
	
	conn<-file(fileName,open="r")
	suppressWarnings(linn<-readLines(conn))
	skip <- 1
	while(regexpr("accession", ignore.case = T,as.character(linn[skip])) == -1 ){
	#while(regexpr("Bio Sample",as.character(linn[skip])) == -1 ){
		skip <- skip +1 
	}
	close(conn)
	
	return(skip-1)
}

### 6-plex or 10-plex
#' @export
.getNbPlex <- function(fileName){
	
	### parse data
	res <- read.csv(fileName,sep="\t",skip=.getSkipLineNb(fileName),nrows=1)
	return(length(grep("^Normalized",names(res))))
}


### get file type
#' @export
.getFileType <- function(file){
	
	if( sum(grepl("Modifications",names(read.csv(file,skip=2, nrows=1))))  > 0 ){
		return("ProgenesisPeptide")
	}else if( sum(grepl("Retention.time..min",names(read.csv(file,skip=2, nrows=1))))  > 0 ){
		return("ProgenesisFeature")
	}else if(sum( grepl("Peptide.count",names(read.csv(file,skip=2, nrows=1)) )) > 0  ){
		return("ProgenesisProtein")
	}else if((grepl("^Identification.Criteria",names(read.csv(file,skip=2, nrows=1))) > 0) &
			(grepl("^Experiment",names(read.csv(file,skip=1, nrows=0))) > 0)){
		return("ScaffoldTMT")
	}else if(sum(grepl("^LFQ",names(read.csv(file,allowEscapes=T, check.names=F,sep="\t",nrows=1)))) > 0
			& 	sum(grepl("Fasta headers",names(read.csv(file,allowEscapes=T, check.names=F,sep="\t",nrows=1)))) > 0
			){ 
		return("MaxQuantProteinGroup") 
	}else if(F){ 
		return("GenericCSV") 
	}else{
		return("Unknown")	
	}
}


#' Parse scaffold output .xls file (RAW export)
#' @param file path to Scaffold file
#' @param expDesign experimental design data.frame
#' @param keepFirstAcOnly TRUE/FALSE If multiple ACs in Accession.Numbers filed. Then keep the first one only
#' @param isPurityCorrect do purity correction
#' @return ExpressionSet object
#' @export
#' @import Biobase
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
parseScaffoldRawFile <- function(file, expDesign=expDesign,keepFirstAcOnly=FALSE,isPurityCorrect=T){
	
	### parse data
	res <- read.csv(file,sep="\t",skip=.getSkipLineNb(file))
	if(keepFirstAcOnly) res$Accession.Numbers <- gsub("\\,.*","",as.character(res$Accession.Numbers))
	res$Accession.Numbers <- gsub("\\[","",as.character(res$Accession.Numbers))
	res$Accession.Numbers <- gsub("\\]","",as.character(res$Accession.Numbers))
	
	###filter out con and decoy
	### HACK, to make sure GFP proteins are not discarded
	res$Accession.Numbers <- gsub(".*P42212.*","sp|P42212|GFP",res$Accession.Numbers)
	
	
	# CREATE EXPRESSION SET	
	nbPlex <- length(grep("^Normalized",names(res)))
	expMatrix <- as.matrix(res[, 10:(10+(nbPlex-1)) ])
	
	if(isPurityCorrect){
		expMatrix <- purityCorrectTMT(expMatrix,getImpuritiesMatrix(nbPlex))
	}
	
	### re-order and exclude channels  
#	expMatrix <- expMatrix[,as.numeric(rownames(expDesign))]
#	colnames(expMatrix) <- rownames(expDesign)
	
	# set 0 features to NA
	expMatrix[expMatrix == 0] <- NA
	
	# discard features where all intensities are NA (i.e. 0)
	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))
	
	### re-order and exclude channels  
	expMatrix <- as.matrix(expMatrix[!allColNA ,as.numeric(rownames(expDesign))])
	colnames(expMatrix) <- rownames(expDesign)
	
	featureAnnotations <- data.frame(
			proteinName=res$Accession.Numbers
			,proteinDescription=res$Protein.Name 
			,peptide=res$Peptide.Sequence
			#,idScore=res$Score
			#,mass=res$Acquired.Plus.One.Mass-1
			#,pMassError=res$Mass.error..ppm.
			,mz=((as.numeric(as.character(res$Acquired.Plus.One.Mass)) -1)+res$Charge)/res$Charge
			#,retentionTime=res$Retention.time..min.
			,charge=as.numeric(res$Charge)
			,spectrumName = res$Spectrum.Name
			#,ptm=res$Variable.modifications...position..description.
			,isNormAnchor=rep(T,nrow(res))
			,isFiltered=rep(F,nrow(res))
	)
	
	
	featureAnnotations <- featureAnnotations[!allColNA,]
	
	### strip off added .1  A11.03216.1 -> A11.03216
	#colnames(expMatrix) <- gsub("\\.1$","",colnames(expMatrix))

	
	eset <- createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations)
	rownames(eset) <- fData(eset)$spectrumName
	return(eset)
}

### maxQuant protein group parser

#' Parse MaxQuant Protein Group Txt
#' @param file path to MaxQuant Protein txt file
#' @param expDesign experimental design data.frame
#' @param method auc (area under curve) or spc (spectral count)
#' @return ExpressionSet object
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @seealso \code{\link[Biobase]{ExpressionSet}}
#' @examples print("No examples")
parseMaxQuantProteinGroupTxt <- function(file=file,expDesign=expDesign, method="auc"){
	
	res <- read.csv(file,allowEscapes=T, check.names=F,sep="\t")
	
	### get
	if(method == "auc"){
		expMatrix <- as.matrix(res[,grepl("^LFQ",names(res))])
	}else if(method == "spc"){
		expMatrix <- as.matrix(res[,grepl("^MS\\/MS Count",names(res))])
	}else{
		stop("Unknown parsing method:", method )
	}
	
	# set 0 features to NA
	expMatrix[expMatrix == 0] <- NA
	# rename runs 1:X
	colnames(expMatrix) <- 1:ncol(expMatrix)
	# discard features where all intensities are NA (i.e. 0)
	allColNA <-  as.vector(apply(expMatrix,1,function(r){ return(sum(!is.na(r)) == 0)}))
	
#	names(res)
#	[1] "Protein IDs"                         
#	[2] "Majority protein IDs"                
#	[3] "Peptide counts (all)"                
#	[4] "Peptide counts (razor+unique)"       
#	[5] "Peptide counts (unique)"             
#	[6] "Protein names"                       
#	[7] "Gene names"                          
#	[8] "Fasta headers"                       
#	[9] "Number of proteins"                  
#	[10] "Peptides"                            
#	[11] "Razor + unique peptides"             
#	[12] "Unique peptides"                     
#	[13] "Peptides 43"                         
#	[14] "Peptides 44"                         
#	[15] "Peptides 45"                         
#	[33] "Razor + unique peptides 43"          
#	[34] "Razor + unique peptides 44"          
#	[35] "Razor + unique peptides 45"          
#	[53] "Unique peptides 43"                  
#	[54] "Unique peptides 44"                  
#	[55] "Unique peptides 45"                  
#	[73] "Sequence coverage [%]"               
#	[74] "Unique + razor sequence coverage [%]"
#	[75] "Unique sequence coverage [%]"        
#	[76] "Mol. weight [kDa]"                   
#	[77] "Sequence length"                     
#	[78] "Sequence lengths"                    
#	[79] "Q-value"                             
#	[80] "Sequence coverage 43 [%]"            
#	[81] "Sequence coverage 44 [%]"            
#	[82] "Sequence coverage 45 [%]"            
#	[100] "Intensity"                           
#	[101] "Intensity 43"                        
#	[102] "Intensity 44"                        
#	[103] "Intensity 45"                        
#	[121] "LFQ intensity 43"                    
#	[122] "LFQ intensity 44"                    
#	[123] "LFQ intensity 45"                    
#	[141] "MS/MS Count 43"                      
#	[142] "MS/MS Count 44"                      
#	[143] "MS/MS Count 45"                      
#	[161] "MS/MS Count"                         
#	[162] "Only identified by site"             
#	[163] "Reverse"                             
#	[164] "Potential contaminant"               
#	[165] "id"                                  
#	[166] "Peptide IDs"                         
#	[167] "Peptide is razor"                    
#	[168] "Mod. peptide IDs"                    
#	[169] "Evidence IDs"                        
#	[170] "MS/MS IDs"                           
#	[171] "Best MS/MS"                          
#	[172] "Oxidation (M) site IDs"              
#	[173] "Oxidation (M) site positions"      
	
	row.names(expMatrix) <- res[,"Protein IDs"]
	
	if( !all(rownames(expDesign) %in% colnames(expMatrix)) ){
		print(expDesign)
		stop("Invalid ExpDesign")
	}
	
	featureAnnotations <- data.frame(
			proteinName=res[,"Protein IDs"]
			,proteinDescription=res[,"Fasta headers"]
			,idScore=res[,"Q-value"]
			#	,nbPeptides=res[,"Peptides used for quantitation"] ### old Progenesis 
			,isDecoy=res[,"Reverse"] == "+"
			,nbPeptides=res[,"Peptides"] # Progensis QI
			,isNormAnchor=rep(T,nrow(expMatrix))
			,isFiltered=res[,"Reverse"] == "+"
			,row.names=res[,"Protein IDs"])
	
	featureAnnotations <- featureAnnotations[!allColNA,]
	
	
	### re-order and exclude channels 
	expMatrix <- as.matrix(expMatrix[!allColNA ,rownames(expDesign)])
	colnames(expMatrix) <- rownames(expDesign)
	
	return(createExpressionDataset(expressionMatrix=expMatrix,expDesign=expDesign,featureAnnotations=featureAnnotations))
}


#' Parse scaffold PTM Spectrum Report
#' @param file path to Scaffold file
#' @return data.frame
#' @export
#' @note  No note
#' @details No details
#' @references NA 
#' @examples print("No examples")
parseScaffoldPTMReport <- function(file){
	
	d <- read.csv(file,sep="\t")
	# Note that peptides with multiple protein assignments are listed once per protein.
	# Howerver, only one peptide per spectrum is listed
	#d <- read.csv(file,sep="\t",skip=.getSkipLineNb(file))[,c(4,5,9,10,16,20)]
	d <- read.csv(file,sep="\t")[,c(3,4,5,9,10,16,20)]
	d <- unique(d)
	rownames(d) <- d$Spectrum.Name
	#	d <- d[,1:5]
	#	names(d) <- c("ptm","ptmLocProb","idScore","ptmLocMascotConfidence","pMassError")
	
	d <- d[,1:(ncol(d) - 1)]
	names(d) <- c("ptmPeptide","ptm","ptmLocProb","idScore","ptmLocMascotConfidence","pMassError")
	
	### add column nbPtmsPerPeptide
#	a <- c("dfgdsg,dfsd,s","sdfsd,we","we","")
#	unlist(lapply(strsplit(a,""),function(t){return(sum(grepl("\\,",t)) + (length(t)> 0))}))
#	
	d$ptm <- as.character(d$ptm )
	#S8 Phospho, S30 Phospho, M35 Oxidation (number of ',' +1 )
	d <- cbind(d,nbPtmsPerPeptide= unlist(lapply(strsplit(d$ptm,""),function(t){return(sum(grepl("\\,",t)) + (length(t)> 0))})) )

	return( d )
}


#' Add scaffold ptm annotaitons to tmt experiment
#' @param eset ExpressionSet
#' @param file path to Scaffold file
#' @return ExpressionSet object
#' @export
#' @import Biobase
#' @note  No note
#' @references No references
#' @examples print("No examples")
addScaffoldPTMFAnnotations <- function(eset,file){
	
	scaffoldPTMFData <- parseScaffoldPTMReport(file)
	
	spectrumNameIntersect <- intersect(fData(eset)$spectrumName,rownames(scaffoldPTMFData))
	# conflicting peptide entries
	misMatch <- fData(eset)[spectrumNameIntersect,]$peptide != toupper(scaffoldPTMFData[spectrumNameIntersect,]$ptmPeptide)
	# discard mismatches from intersect and Scaffold PTM df
	spectrumNameIntersect <- spectrumNameIntersect[!misMatch]
	scaffoldPTMFData <- scaffoldPTMFData[spectrumNameIntersect,]
	
	missingEntriesFrac <- 1- length(spectrumNameIntersect) / nrow(eset)
	if(missingEntriesFrac > 0)cat("INFO: ", round((missingEntriesFrac * 100)), "% of TMT spectra not listed in Scaffold PTM export \n" )
	
	fData(eset) <- data.frame(fData(eset), scaffoldPTMFData[rownames(fData(eset)),], row.names=rownames(fData(eset)))

	#set NA scores 0
	#fData(eset)$idScore[is.na(fData(eset)$idScore)] <- 0
	
	return(eset)
	
}




################################## TMT END ##############################






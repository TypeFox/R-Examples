# TODO: Add comment
# 
# Author: erikahrne
###############################################################################

# set filter, filer = data.frame all TRUE/FALSE 
#' @export
.setFilter <- function(eset=eset, filter=filter){
	
	fData(eset)$isFiltered <- apply(filter,1,function(f){ return(sum(f,na.rm=T) > 0 ) })
	return(eset)
}	

### add ptm coord to eset fData (data frame  listing phospho coordinates and motifX)
#' @export
#' @importFrom seqinr read.fasta
#' @importFrom utils setTxtProgressBar txtProgressBar
.addPTMCoord <- function(eset,proteinDB,motifLength=4, isProgressBar=F,format=1){
	
	ptmCoordDf <- data.frame()
	
	### progress bar
	pbSum <- txtProgressBar(min = 0, max = nrow(eset), style = 3)
	
	for(i in 1:nrow(eset)){
		
		### increment progressbar
		if(isProgressBar) setTxtProgressBar(pbSum, i)
		
		modifAnnot <-  as.character(fData(eset)$ptm[i])
		
		### if peptide is modified
		if((nchar(modifAnnot) > 0) & !is.na(modifAnnot) ){
			peptide <- as.character(fData(eset)$peptide[i])
			proteinAC <- as.character(fData(eset)$proteinName[i])
			
			modifCoord <- NA
			motifX <- NA
			
			proteinSeq <- proteinDB[proteinAC]
			
			if(proteinSeq == "NULL" ){
				warning(proteinAC," NOT FOUND IN PROTEIN FASTA","\n")
#					modifCoord <- "Err"
#					motifX <- "Err"
			}else{
				
				modifCoord <- getModifProteinCoordinates(modifAnnot,peptide,proteinSeq, format=format)
				motifX <- paste(getMotifX(modifCoord,peptide, proteinSeq,motifLength),collapse=",")
				modifCoord <- paste(modifCoord,collapse=",")
			}
			
			
		}else{
			modifCoord <- ""
			motifX <- ""
		}
		
		ptmCoordDf <- rbind(ptmCoordDf,data.frame(modifCoord,motifX))
		
	}
	
	# close progress bar
	setTxtProgressBar(pbSum, i)
	close(pbSum)
	
	fData(eset) <- cbind(fData(eset),ptmCoordDf)
	
	return(eset)
	
}

### more to be added
#' @export
.getProteaseRegExp <- function(protease="trypsin"){
	
	if(protease == "trypsin"){
		return("[KR](?!P)")
	}else if(protease == "lys-c"){
		return("K(?!P)")
	}else{
		stop("Unknown Protease:", protease )
	}
	
}

# @details returns a dataframe lising all unique motifs and the corresponding ptm. multiple entries will be created  for multiply modified peptides
#' @export
.getUniquePtmMotifs <- function(eset, format=1){
	
	eset <- eset[!fData(eset)$isFiltered & (nchar(as.character(fData(eset)$ptm)) > 0),]
	
	### has to be done this way as sometimes the matching AC is not found in the database...
	motifXTag <- c()
	ptmTag <- c()
	
	if(format ==1){
		for(i in 1:nrow(eset)){
			
			ptms <- gsub("\\[[0-9]*\\] {1,}","",unlist(strsplit(as.character(fData(eset)$ptm[i]),"\\|")))
			mX <- unlist(strsplit(as.character(fData(eset)$motifX[i]),"\\,"))
			
			ptmTag <- c(ptmTag,ptms)
			
			#if(length(ptms) > 0){
			if(length(ptms) == length(mX) ){
				motifXTag <- c(motifXTag,mX)
				
			}else{
				motifXTag <- c(motifXTag,rep(NA,length(ptms)))
			}
			#}
		}
	}else{ # scaffold PTM
		for(i in 1:nrow(eset)){
			ptms <- gsub("[0-9]","",unlist(strsplit(as.character(fData(eset)$ptm[i]),"\\, ")))
			mX <- unlist(strsplit(as.character(fData(eset)$motifX[i]),"\\,"))
			
			ptmTag <- c(ptmTag,ptms)
			
			#if(length(ptms) > 0){
			if(length(ptms) == length(mX) ){
				motifXTag <- c(motifXTag,mX)
				
			}else{
				motifXTag <- c(motifXTag,rep(NA,length(ptms)))
			}
			#}
		}
	}
	
	
	
	sel <- !is.na(motifXTag) & grepl("\\*",motifXTag)
	return(unique(data.frame(ptm = ptmTag[sel], motif = motifXTag[sel])))
}

### Detectable/Identifiable peptide lengths (600 - 4000 Da -> 600/110 - 4000/110 ->) 5 - 36 AA's
#' Get number peptides passing defined length criteria
#' @param peptides list of peptides 
#' @param peptideLength vector of two integers defining peptide length range
#' @return integer corresponding to number of detectable peptides
#' @export
#' @note  No note
#' @details No details
#' @examples print("No examples")
getNbDetectablePeptides <- function(peptides, peptideLength=c(5,36)){
	okLength <- (nchar(peptides) >= peptideLength[1]) & (nchar(peptides) <= peptideLength[2])
	return(sum(okLength))
}

#' Digest protein 
#' @param proteinSeq protein sequence 
#' @param proteaseRegExp protease Regular Expression
#' @param nbMiscleavages default 0
#' @return vector of peptides
#' @export
#' @note  No note
#' @details No details
#' @examples print("No examples")
getPeptides <- function(proteinSeq,proteaseRegExp=.getProteaseRegExp("trypsin"),nbMiscleavages=0){
	
	allAA <- as.vector(unlist(strsplit(proteinSeq,"")))
	
	### get all full cleaved peptides without cleaved residue
	fcPeptides <- strsplit(proteinSeq,proteaseRegExp,perl=TRUE)[[1]]
	
	### add cleaved residue
	matchPos <- gregexpr(proteaseRegExp,proteinSeq,perl=TRUE)[[1]]
	separator <- allAA[ matchPos ]
	###if c-term peptide isn't tryptic
	if(length(separator) < length(fcPeptides)) separator <- c(separator,"")
	fcPeptides <- paste(fcPeptides,separator,sep="")
	fcPeptides <- fcPeptides[nchar(fcPeptides) > 0 ]
	
	### if no mis-cleavages that's it
	if(nbMiscleavages == 0){
		return(fcPeptides)
	}
	
	allPeptides <- c()
	### handle miscleavages
	for(i in 1:length(fcPeptides)){
		#cat(fcPeptides[i]," ",i," -------------\n")
		for(j in i:(i+nbMiscleavages)){
			if(j <= length(fcPeptides)){
				pept <- paste(fcPeptides[i:j],collapse="")
				#cat(j," ---",pept,"\n")
				allPeptides <- c(allPeptides,pept)
			}
		}
	}
	
	return(allPeptides[nchar(allPeptides) > 0])
	
}


#' Add identification leve q-values to ExpressionSet (calculated based on target-decoy score distribution)
#' @param eset ExpressionSet
#' @return ExpressionSet object
#' @details if ptm column is part if the ExpressionSet q-values are calculated seperately for modified and non-modified features
#' @export
#' @note  No note
#' @details No details
#' @seealso \code{\link{getIdLevelQvals}}
#' @examples print("No examples")
addIdQvalues <- function(eset=eset){
	
	### only calculate q-values for non-filtered features
	sel <- rep(T,nrow(eset))
	if(!is.null(fData(eset)$isFiltered)){
		sel <- !fData(eset)$isFiltered
	}
	
	# init vector to be added
	idQValue <- rep(1,nrow(eset))
	
	# if ptm column exist seperate calculate qvals sepearately for modif and nomd features
	if(!is.null(fData(eset)$ptm)){
		
		isMod <- nchar(as.character(fData(eset)$ptm)) > 0
		
		# get qvalues for modified features
		if(sum(isMod  & sel ) > 0){
			idQValue[isMod & sel ] <- getIdLevelQvals(fData(eset)$idScore[isMod & sel],isDecoy(fData(eset)$proteinName[isMod & sel]))
		}
		# get qvalues for non-modified features
		if(sum(!isMod  & sel ) > 0){
			idQValue[!isMod & sel] <- getIdLevelQvals(fData(eset)$idScore[!isMod & sel],isDecoy(fData(eset)$proteinName[!isMod & sel]))
		}
	}else{# get qvalues if ptm column is missing
		idQValue[sel] <- getIdLevelQvals(fData(eset)$idScore[sel] ,isDecoy(fData(eset)$proteinName[sel]))
	}
	
	# add q-value vector
	fData(eset)$idQValue <- idQValue
	
	return(eset)
	
}


### check if protein is a contaminant
#' Check if protein is a contaminant entry
#' @param ac vector of protein accession numbers 
#' @return vector TRUE/FALSE
#' @export
#' @note  No note
#' @details contanminants proteins are typically annotated as: CON_P0000
#' @references NA 
#' @examples print("No examples")
isCon <- function(ac){
	ac <- as.character(ac)
	return(regexpr("\\Wcon_",ac,ignore.case=TRUE, perl=TRUE) > -1 | regexpr("^con_",ac,ignore.case=TRUE, perl=TRUE) > -1 | regexpr("_con",ac,ignore.case=TRUE, perl=TRUE) > -1 )
}

### check if protein is a Decoy
#' Check if protein is a decoy entry
#' @param ac vector of protein accession numbers 
#' @return vector TRUE/FALSE
#' @export
#' @note  No note
#' @details decoy proteins are typically annotated as: REV_P0000
#' @references NA 
#' @examples print("No examples")
isDecoy <-function(ac){
	
	return(((regexpr("^rev_",ac,ignore.case=TRUE) > -1)
						| (regexpr("^decoy_",ac,ignore.case=TRUE) > -1)
						| (regexpr("^revipi",ac,ignore.case=TRUE) > -1)
						| (regexpr("^reverse_",ac,ignore.case=TRUE) > -1)))
	
}




### calculate id level qvals based on target-decoy score distributions
#' Calculates identification level q-values based on target-decoy score distributions
#' @param scores peptide/protein identficationscore
#' @param isDecoy vector of TRUE/FALSE 
#' @return vector of q.values
#' @export
#' @note  No note
#' @details q-value = (Nb. Decoy Entries at idScore Threshold S*) / (Nb. Target Entries at idScore Threshold S). (* idScore >= S)
#' @references NA 
#' @examples print("No examples")
getIdLevelQvals <- function(scores, isDecoy){
	
	if(sum(isDecoy) == 0){
		return(rep(0,length(scores)))
	}
	
	targetScores = scores[!isDecoy]
	decoyScores = scores[isDecoy]
	
	qvals = c()
	for(score in scores){
		
		qval = sum(decoyScores >= score) / sum(targetScores >= score)
		qvals = c(qvals,qval)
	}
	
	return(qvals)
	
}

### Get score cutoff for a given fdr cut-off
#' Get score cutoff for a given fdr cut-off
#' @param scores peptide/protein identficationscore
#' @param isDecoy vector of TRUE/FALSE 
#' @param fdrCutOff [0,1]
#' @return scoreCutoff 
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
getScoreCutOff <- function(scores,isDecoy,fdrCutOff = 0.01){
	return(min(scores[(getIdLevelQvals(scores ,isDecoy ) < fdrCutOff)]))
}

### create motif-x peptide annotation
#' Create motif-x peptide annotation
#' @param modifPos vector positions
#' @param peptide peptide sequence
#' @param proteinSeq protein sequence
#' @param motifLength motif flanking sequence
#' @return vector of motifs 
#' @export
#' @note  No note
#' @details motif-x example PGDYS*TTPG
#' @references NA 
#' @examples print("No examples")
getMotifX <- function(modifPos,peptide,proteinSeq, motifLength=4){
	
	#modifPos <- getModifProteinCoordinates(modifAnnot,peptide,proteinSeq)
	
	if( modifPos[1] > -1){
		
		motifX <- c()
		
		for(mp in modifPos){
			
			prefix <- substr(proteinSeq,mp-motifLength,mp-1)
			suffix <- substr(proteinSeq,mp+1,mp+motifLength)
			modifAA <-  paste(substr(proteinSeq,mp,mp),"*",sep="")
			
			mx <- paste(prefix,modifAA,suffix,sep="")
			motifX <- c(motifX,mx)
		}
		
		return(motifX)
		
	}else{
		return(NA)
	}
}

### Get modification coordinates on protein
#' Get modification coordinates on protein
#' @param modifAnnot modifcation as annotated by progenesis. E.g. '[15] Phospho (ST)|[30] Phospho (ST)'
#' @param peptideSeq peptide sequence
#' @param proteinSeq protein sequence
#' @param format c(1,2) 1. progenesis 2. scaffold 
#' @return vector of protein coordinates (mmodification residue number) 
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
getModifProteinCoordinates <- function(modifAnnot,peptideSeq,proteinSeq, format=1){
	
	peptideSeq <- as.character(peptideSeq)
	modifPos <- c()
	
	if(format == 1){
		### parse modifAnnot -> modif position(s) on peptide positions
		modifList <- strsplit(as.character(modifAnnot),"\\|")[[1]]
		
		for(m in modifList){
			
			if(grepl("N\\-term",m)){
				m <- 1
			}else if(grepl("C\\-term",m)){
				m <- nchar(peptideSeq)
			}else{
				m <- gsub("\\[","",m)
				m <- as.numeric(gsub("\\].*","",m))
			}
			modifPos <- c(modifPos,m)
		}
	}else if(format == 2){
		
		modifList <- strsplit(as.character(modifAnnot),"\\, ")[[1]]
		
		for(m in tolower(modifList)){
			
			### @ TODO check how n-term and c-term modificaitons are annotated
			if(grepl("n\\-term",m)){
				m <- 1
			}else if(grepl("c\\-term",m)){
				m <- nchar(peptideSeq)
			}else{
				m <- gsub("[a-z A-Z]","",m)
			}
			modifPos <- c(modifPos,m)
		}
				
		modifPos <- as.numeric(modifPos)
		
	}else{
		stop("ERROR: getModifProteinCoordinates - Unknown PTM format", format)
	}
	
	### get peptide position on protein (first matching position if duplicates)
	pepStartPos <- regexpr(peptideSeq,proteinSeq)[1]
	
	if(pepStartPos < 0){
		cat("\nERROR - getModifProteinCoordinates -:",peptideSeq," not matching current protein sequence (check provided protein sequence db)","\n" )	
		return(pepStartPos)
	}
	
	return((modifPos+pepStartPos-1))
	
}

#' Get number of mis-cleavages perp peptide
#' @param peptide character vector
#' @param protease regular expression
#' @return vector ofintegers
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
getNbMisCleavages <-function(peptide, protease="trypsin"){
	
	return(
			unlist(lapply(peptide,function(t){
								cTermRegExp <- paste(.getProteaseRegExp(protease),"$",sep="")	
								### discard 1 if cTerm matching regExp
								
								#	print(gregexpr(.getProteaseRegExp(protease),t,perl=T)[[1]])
								
								return((sum(gregexpr(.getProteaseRegExp(protease),t,perl=T)[[1]] > -1 )) - grepl(cTermRegExp,t,perl=T))
								
							}))
	)
}

#' Get number of peptides per protein
#' @param eset ExpressionSet
#' @return table
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
getNbPeptidesPerProtein <- function(eset){
	
	sel <- !fData(eset)$isFiltered
	return(sort(table(unique(data.frame(fData(eset[sel,])$proteinName,fData(eset[sel,])$peptide))[,1]),decreasing=T))
}

#' Set nbPeptides coulmn of featureData
#' @param eset ExpressionSet
#' @return eset
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
setNbPeptidesPerProtein <- function(eset){
	fData(eset)$nbPeptides	<- getNbPeptidesPerProtein(eset)[as.character(fData(eset)$proteinName)]
	return(eset)
}

#' Get number of spectra per protein
#' @param eset ExpressionSet
#' @return table
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
getNbSpectraPerProtein <- function(eset){
	
	sel <- !fData(eset)$isFiltered
	return(sort(table(fData(eset[sel,])$proteinName),decreasing=T))
	
}

#' Set nbPeptides coulmn of featureData
#' @param eset ExpressionSet
#' @return eset
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
setNbSpectraPerProtein <- function(eset){
	fData(eset)$nbSpectra	<- getNbSpectraPerProtein(eset)[as.character(fData(eset)$proteinName)]
	return(eset)
}

#' Get modification coordinates on protein
#' @param d numeric vector
#' @param nbSd range spanning number of sd frmo mean
#' @return vector range boundries
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
getMeanCenteredRange <- function(d,nbSd=4){
	return(c(mean(d,na.rm=T) - nbSd*sd(d,na.rm=T),mean(d,na.rm=T) + nbSd*sd(d,na.rm=T)))
}

#' Check if ACs are in "non-stripped" uniprot format e.g. "sp|Q8CHJ2|AQP12_MOUSE" 
#' @param acs accession numbers
#' @return boolean TRUE/FALSE
#' @export
#' @note  No note
#' @details TRUE if less than 10% of ACs contain a "|" character
#' @references NA 
#' @examples print("No examples")
isStrippedACs <-function(acs){
	acs <- as.character(acs)
	return( (sum(grepl("\\|",acs)) / length(acs)) < 0.9 )
}

#' strip uniprot format e.g. "sp|Q8CHJ2|AQP12_MOUSE" ->  Q8CHJ2
#' @param acs accession numbers
#' @return vector character
#' @export
#' @note  No note
#' @details TRUE if less than 10% of ACs contain a "|" character
#' @references NA 
#' @examples print("No examples")
stripACs <- function(acs){
	acs <- gsub("[a-z]{1,4}\\|","",acs)
	acs <- gsub("\\|.*","",acs)
	return(acs)
	
}

#' Get amino acid coordinates on protein
#' @param peptideSeq peptide sequence
#' @param proteinSeq protein sequence
#' @param aaRegExpr target AA reg exp
#' @return vector of protein coordinates (mmodification residue number) 
#' @export
#' @note  No note
#' @details NA
#' @references NA 
#' @examples print("No examples")
getAAProteinCoordinates <- function(peptideSeq,proteinSeq, aaRegExpr="[STY]"){
	
	peptideSeq <- as.character(peptideSeq)
	
	### parse modifAnnot -> modif position(s) on peptide positions
	#modifList <- strsplit(as.character(modifAnnot),"\\|")[[1]]
	modifPos <- as.vector(unlist(gregexpr(aaRegExpr,peptideSeq)[[1]]))
	
	### get peptide position on protein (first matching position if duplicates)
	pepStartPos <- regexpr(peptideSeq,proteinSeq)[1]
	
	if(pepStartPos < 0){
		cat("\nERROR - getModifProteinCoordinates -:",peptideSeq," not matching current protein sequence (check provided protein sequence db)","\n" )	
		return(pepStartPos)
	}
	
	return((modifPos+pepStartPos-1))
	
}
setClass("CLUMP",
	representation = representation(
						strEqcCommand		=	"character",
						rcdCriterion		=	"character",
						colInMarker			=	"character",
						colClump			=	"character",
						numPvalLim			=	"numeric",
						numPosLim			=	"numeric",
						numR2Lim			=	"numeric",
						filePLINK			=	"character",
						fileBfile			=	"character",
						blnAddClumpInfo 	= 	"logical",
						strTag				=	"character"
						),
	prototype = prototype(
						strEqcCommand		=	"",
						rcdCriterion		=	"",
						colInMarker			=	"",
						colClump			=	"",
						numPvalLim			=	5e-8,
						numPosLim			=	500000,
						numR2Lim			=	0.2,
						filePLINK			=	"",
						fileBfile			=	"",
						blnAddClumpInfo 	= 	FALSE,
						strTag				=	""
						)
	#contains = c("EcfReader")
)

setGeneric("setCLUMP", function(object) standardGeneric("setCLUMP"))
setMethod("setCLUMP", signature = (object = "CLUMP"), function(object) {
	
	aEqcSlotNamesIn = c("rcdCriterion","colInMarker","colClump","numPvalLim","numPosLim","numR2Lim","filePLINK","fileBfile","blnAddClumpInfo","strTag")

	objEqcReader <- EqcReader(object@strEqcCommand,aEqcSlotNamesIn)
	
	if(length(objEqcReader@lsEqcSlotsOut) > 0) {
		for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
			tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
			tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
			
			if(all(!is.na(tmpSlotVal))) slot(object, tmpSlot) <- tmpSlotVal
		}
	}
	return(object)
})

#############################################################################################################################
validCLUMP <- function(objCLUMP) {
	
	### Valid with GWADATA?
	
	if(objCLUMP@rcdCriterion == "") 
		warning(paste("EASY WARNING:CLUMP\n No criterion rcdCriterion defined. All data will be used for clumping. Buffering may need some time!" ,sep="" ))
		#cat(paste(" EASY WARNING:CLUMP\n No criterion rcdCriterion defined. All data will be used for clumping.", sep=""))
	if(objCLUMP@colInMarker == "") 
		stop(paste(" EASY ERROR:CLUMP\n No column colInMarker defined. Please set colInMarker.", sep=""))
	if(objCLUMP@colClump == "") 
		stop(paste(" EASY ERROR:CLUMP\n No column colClump defined. Please set colClump.", sep=""))
	
	if(!file.exists(objCLUMP@filePLINK))
		stop(paste("EASY ERROR:CLUMP\n File filePLINK\n ",objCLUMP@filePLINK,"\n does not exist.", sep=""))

	if(!file.exists(paste(objCLUMP@fileBfile,".bed",sep="")))
		stop(paste("EASY ERROR:CLUMP\n File fileBfile\n ",objCLUMP@fileBfile,"\n does not exist.", sep=""))
	
	return(TRUE)
}
CLUMP.GWADATA.valid <- function(objCLUMP, objGWA) {
	
	isAv <- objCLUMP@colInMarker %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:CLUMP\n Defined column colInMarker \n",objCLUMP@colInMarker, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
	isAv <- objCLUMP@colClump %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:CLUMP\n Defined column colClump \n",objCLUMP@colClump, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))
	
}

#############################################################################################################################
CLUMP.run <- function(objCLUMP, objGWA, objREPORT) {

	rcdCriterion 	<- objCLUMP@rcdCriterion
	colInMarker 	<- objCLUMP@colInMarker
	colClump 		<- objCLUMP@colClump
	numPvalLim 		<- objCLUMP@numPvalLim
	numPosLim 		<- objCLUMP@numPosLim
	numR2Lim 		<- objCLUMP@numR2Lim
	filePLINK 		<- objCLUMP@filePLINK
	fileBfile 		<- objCLUMP@fileBfile
	blnAddClumpInfo	<- objCLUMP@blnAddClumpInfo
	strTag 			<- objCLUMP@strTag
	
	objRCD 	<- RCD(rcdCriterion)
	out 	<- RCD.eval(objRCD, objGWA)
	
	numClumpCrit = length(which(out))
	
	objGWA.clumpin <- GWADATA.getrows(objGWA, which(out))
	
	is_na = is.na(GWADATA.getcol(objGWA.clumpin, objCLUMP@colClump))

	objGWA.clumpin@tblGWA <- 	objGWA.clumpin@tblGWA[!is_na,]

	numClumpNA = length(which(is_na))
	
	objGWA.clumpin <- GWADATA.renamecol(objGWA.clumpin, colInMarker, "SNP")
	#objGWA.clumpin <- GWADATA.renamecol(objGWA.clumpin, colInPval, "P")
	
	#### Buffer
	#GWADATA.write(objGWA.clumpin, strSuffix = ".PLINK_TMP")
	
	#fileInTmp	<- paste(objGWA.clumpin@pathOut,"/PLINK_IN.tmp",sep="")
	fileBasePlink = paste(objGWA.clumpin@pathOut,"/",objGWA@fileInShortName,sep="")
	fileBasePlink = sub(".txt","",fileBasePlink)
	fileBasePlink = sub(".gz","",fileBasePlink)
	if(nchar(strTag)>0) fileBasePlink = paste(fileBasePlink,".",strTag,sep="")
	
	fileInPlink	<- paste(fileBasePlink ,".plink_in",sep="")
	write.table(objGWA.clumpin@tblGWA, fileInPlink, row.names=F, quote=F, sep="\t")		
	fileOutPlink	<- paste(fileBasePlink ,".plink_out",sep="")
	
	i = 1
	while(file.exists(paste(fileOutPlink,".clumped",sep=""))) {
		#fileOut <- paste(fileOutBase,".",i,".txt.gz",sep="")
		fileOutPlink <- paste(fileOutPlink,".",i,sep="")
		i = i + 1
	}
	
	#### PLINK clump
	## bp to kb
	numPosLim <- round(numPosLim/1000)
	
	plink_Call = paste(filePLINK,"--noweb","--bfile",fileBfile,"--clump",fileInPlink,"--out",fileOutPlink,"--clump-kb",numPosLim, "--clump-r2",numR2Lim,"--clump-p1",numPvalLim,"--clump-p2",numPvalLim,"--clump-field",colClump,"--clump-verbose")
	system(plink_Call)
	
	## creates three files fileOutPlink.clumped, fileOutPlink.log, fileOutPlink.nosex

	vDir<-scan(file = paste(fileOutPlink,".clumped",sep=""), what=character(0), n = -1, sep = "\n",quiet=TRUE,blank.lines.skip=FALSE)
	
	strSplitClump = "------------------------------------------------------------------"

	iClumpChange = grep(strSplitClump, vDir)

	numClumps = length(iClumpChange)

	cat(paste("Found", numClumps, "independent clumps ... "))
	cat("\n")
	
	tblClumpAll = data.frame()
	tblClumpIndex = data.frame()

	cat("Reformatting PLINK output ... ")
	cat("\n")
	for(i in 1:numClumps) {
		
		if(i == 1) vDirClump = vDir[1:(iClumpChange[i] - 1)]
		else vDirClump = vDir[(iClumpChange[i-1] + 1):(iClumpChange[i] - 1)]
		
		for(j in 1:length(vDirClump)) while(grepl('  ',vDirClump[j])) vDirClump[j]=gsub('  ',' ', vDirClump[j]) 
		
		########### REWRITE TO BLOCKS DEFINED BY which(vDirClump)=="" !!!!!!
		## JUST DISTINGUISH BLOCKS -> DONE
		
		aiBlocks = which(vDirClump=="")
		
		for(j in 1:(length(aiBlocks)-1)) {
			vDirTemp = vDirClump[aiBlocks[j]:(aiBlocks[j+1])]
			vDirTemp = vDirTemp[vDirTemp!=""]
			
			if(grepl(" CHR F SNP BP P TOTAL NSIG S05 S01 S001 S0001", vDirTemp[1], fixed=T)) {
				vDirClumpIndex = vDirTemp
				g=strsplit(vDirClumpIndex," ")
				g=lapply(g,function(x) x[-1])
				d = data.frame(t(g[[2]]),stringsAsFactors=FALSE)
				names(d) = t(g[[1]])
				#d = cbind(d, "aLociTag" = i)
				d = cbind(d, "aLociTag" = i, "aTopHit" = i)
				
				### P will be added later via acolAdd!
				idxP = match("P", names(d))
				d = d[,-idxP]
				
				tblClumpIndex = rbind(tblClumpIndex, d)
			}
			
			if(grepl(" KB RSQ ALLELES F P ", vDirTemp[1], fixed=T)) {
				vDirTempPart2 = vDirClump[aiBlocks[j+1]:(aiBlocks[j+2])]
				vDirTempPart2 = vDirTempPart2[vDirTempPart2!=""]
				vDirClumpOut = c(vDirTemp,vDirTempPart2)
				vDirClumpOut = gsub("(INDEX) ","",vDirClumpOut,fixed=T)
				vDirClumpOut = vDirClumpOut[vDirClumpOut!=""]
				vDirClumpOut[1] = paste(" SNP",vDirClumpOut[1],sep="")
				
				#vDirClumpOut[1] = gsub("ANNOT", paste(strsplit(clump.annot,",")[[1]],collapse=" "), vDirClumpOut[1])
				#vDirClumpOut = gsub(", "," ",vDirClumpOut,fixed=TRUE)
				
				g=strsplit(vDirClumpOut," ")
				g=lapply(g,function(x) x[-1])
				
				for(k in 2:length(g)) {
					if(length(g[[k]]) != length(g[[1]])) {
						ag = g[[k-1]]
						ag[(length(ag)-length(g[[k]])+1):length(ag)] <- g[[k]]
						#ag[(length(ag)-1):length(ag)] <- g[[k]]
						g[[k]] <- ag
					}
				}
				
				d=data.frame()
				for(k in 2:length(g)) d = rbind(d, data.frame(t(g[[k]]),stringsAsFactors=FALSE))
				names(d) = t(g[[1]])
				aTopHit = rep(NA,dim(d)[1])
				aTopHit[1] = i
				d = cbind(d, "aLociTag" = i, "aTopHit" = aTopHit)
				
				### P will be added later via acolAdd!
				idxP = match("P", names(d))
				d = d[,-idxP]
				
				tblClumpAll = rbind(tblClumpAll, d)
					
			}	
		}	
	}
	
	# tblOut <- tblOut[,-match(c("ALLELES","F"),names(tblOut))]

	tblClumpIndex <- tblClumpIndex[,c("SNP","CHR","BP","TOTAL","aLociTag","aTopHit")]
	tblClumpAll <- tblClumpAll[,-match(c("ALLELES","F"),names(tblClumpAll))]
	## merge files (loci with only one SNP are not listed in tblClumpAll !!!
	## merge by SNP/aLociTag/aTopHit
	tblOutAll <- merge(tblClumpAll,tblClumpIndex,all=TRUE)
	
	### add aNumLocusSNPs to tblOutAll
	aNumLocusSNPs <- rep(NA, length(tblOutAll$aLociTag))
	for(locusNum in unique(na.omit(tblOutAll$aTopHit))) {
		aNumLocusSNPs[tblOutAll$aLociTag == locusNum] <- length(which(tblOutAll$aLociTag == locusNum))
	}
	tblOutAll = cbind(tblOutAll, aNumLocusSNPs)
	if(any(grepl(" not found in dataset", vDir))) {
		aRowNotFound = vDir[which(grepl(" not found in dataset", vDir))]
		aSnpNotFound <- sub(" not found in dataset","",aRowNotFound)
		tblNotFound <- as.data.frame(matrix(NA,length(aSnpNotFound),ncol(tblOutAll)))
		names(tblNotFound) <- names(tblOutAll)
		tblNotFound$SNP <- aSnpNotFound
		tblOutAll <- rbind(tblOutAll, tblNotFound)
	}
	#tblClumpIn = read.table(fileIn, sep="\t",header=TRUE,stringsAsFactors=FALSE)
	objGWA.clumpout 			<- GWADATA()
	objGWA.clumpout@tblGWA 		<- tblOutAll
	objGWA.clumpout@aHeader 	<- names(tblOutAll)
	objGWA.clumpout@aClasses	<- sapply(tblOutAll,class)

	iRename = which(names(objGWA.clumpout@tblGWA)%in%names(objGWA.clumpin@tblGWA) & names(objGWA.clumpout@tblGWA)!="SNP")
	if(length(iRename)>0) {
		for(iRenameTmp in iRename)
			objGWA.clumpout <- GWADATA.renamecol(objGWA.clumpout, names(objGWA.clumpout@tblGWA)[iRenameTmp], paste(names(objGWA.clumpout@tblGWA)[iRenameTmp],".ref",sep=""))
	}

	objGWA.clump 	<- GWADATA.merge(objGWA.clumpin, objGWA.clumpout, blnAll.In = TRUE, blnAll.Add = TRUE, strBy.In = "SNP", strBy.Add = "SNP", strSuffix.In = "", strSuffix.Add ="" )
	objGWA.clump 	<- GWADATA.renamecol(objGWA.clump, "SNP", colInMarker)
	objGWA.clump.x 	<- GWADATA.getrows(objGWA.clump, which(!is.na(objGWA.clump@tblGWA$aTopHit)))
	objGWA.clump.x 	<- GWADATA.removecols(objGWA.clump.x, c("KB","RSQ"))
	
	objGWA.clump.nomatch <- GWADATA.getrows(objGWA.clump, which(is.na(objGWA.clump@tblGWA$aLociTag)))
	
	if(nchar(strTag)>0) strTag = paste(strTag,".",sep="") 
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numClumpCrit",sep=""),numClumpCrit)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numClumpNA",sep=""),numClumpNA)
	objREPORT <- REPORT.addval(objREPORT,paste(strTag,"numClumpLoci",sep=""),nrow(objGWA.clump.x@tblGWA))
	
	if(blnAddClumpInfo) {
		# objGWA.clump.tmp <- GWADATA.getcols(objGWA.clump, c(colInMarker , "aLociTag", "aTopHit"))
		# objGWA.clump.tmp <- GWADATA.renamecol(objGWA.clump.tmp, "aLociTag", paste(strTag,"aLociTag",sep=""))
		# objGWA.clump.tmp <- GWADATA.renamecol(objGWA.clump.tmp, "aTopHit", paste(strTag,"aTopHit",sep=""))
		objGWA.clump.tmp <- GWADATA.getcols(objGWA.clump, c(colInMarker , "aLociTag", "aTopHit", "aNumLocusSNPs"))
		objGWA.clump.tmp <- GWADATA.renamecol(objGWA.clump.tmp, "aLociTag", paste(strTag,"aLociTag",sep=""))
		objGWA.clump.tmp <- GWADATA.renamecol(objGWA.clump.tmp, "aTopHit", paste(strTag,"aTopHit",sep=""))
		objGWA.clump.tmp <- GWADATA.renamecol(objGWA.clump.tmp, "aNumLocusSNPs", paste(strTag,"aNumLocusSNPs",sep=""))
		
		if(paste(strTag,"aLociTag",sep="") %in% objGWA@aHeader) 
			stop("EASY ERROR:CLUMP\nAdding column aLociTag to large data-set failed because column\n",paste(strTag,"aLociTag",sep=""),"\n is already present GWADATA. PLease use --strTag to influence the column name.")
		if(paste(strTag,"aTopHit",sep="") %in% objGWA@aHeader) 
			stop("EASY ERROR:CLUMP\nAdding column aTopHit to large data-set failed because column\n",paste(strTag,"aTopHit",sep=""),"\n is already present GWADATA. PLease use --strTag to influence the column name.")
		if(paste(strTag,"aNumLocusSNPs",sep="") %in% objGWA@aHeader) 
			stop("EASY ERROR:CLUMP\nAdding column aNumLocusSNPs to large data-set failed because column\n",paste(strTag,"aNumLocusSNPs",sep=""),"\n is already present GWADATA. PLease use --strTag to influence the column name.")		
		objGWA <- GWADATA.merge(objGWA,objGWA.clump.tmp, 
								strSuffix.In = "", strSuffix.Add = "", 
								blnAll.In = TRUE, blnAll.Add = TRUE, 
								strBy.In = colInMarker, strBy.Add = colInMarker)
	}
	
	return(list(objGWA, objGWA.clump,objGWA.clump.x,objGWA.clump.nomatch, objREPORT))

}

CLUMP <- function(strEqcCommand){ 
	## Wrapper for class definition
	CLUMPout <- setCLUMP(new("CLUMP", strEqcCommand = strEqcCommand))
	validCLUMP(CLUMPout)
	#CLUMPout.valid <- validCLUMP(CLUMPout)
	return(CLUMPout)
}


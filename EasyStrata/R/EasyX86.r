# rm(list = ls(all = TRUE))
# graphics.off()
# closeAllConnections()
# options(show.error.messages = TRUE)
# options(warn=0)

# fileECF		<-	"/d/GENETICS/ScriptsDBDebug/easy_x/140331_debug_easyqc_85/140404_2_garbagecleaning/2_metalevel_qc.gwa_1.ecf"
# pathClasses <- 	"/d/GENETICS/ScriptsDBDebug/easy_x/bin/"

# library(Cairo)
# library(plotrix)

# # base functions and classes
# source(paste(pathClasses,"clsEqcReader4.r",sep=""))
# source(paste(pathClasses,"clsRCD7.r",sep=""))
# #source(paste(pathClasses,"clsGWADATA41.r",sep=""))
# source(paste(pathClasses,"clsGWADATA43.r",sep=""))
# source(paste(pathClasses,"clsEASYMERGE.r",sep=""))
# source(paste(pathClasses,"fnEasyPLOT10.r",sep=""))
# source(paste(pathClasses,"clsREPORT9.r",sep=""))

# # inherited classes
# source(paste(pathClasses,"clsMERGE11.r",sep=""))
# #source(paste(pathClasses,"clsSPLOT19.r",sep=""))
# source(paste(pathClasses,"clsSPLOT21.r",sep=""))

# # evaluation classes
# source(paste(pathClasses,"clsADDCOL11.r",sep=""))
# #source(paste(pathClasses,"clsADJUSTALLELES13.r",sep=""))
# source(paste(pathClasses,"clsADJUSTALLELES14.r",sep=""))
# source(paste(pathClasses,"clsAFCHECK15.r",sep=""))
# source(paste(pathClasses,"clsANNOTATE8.r",sep=""))
# source(paste(pathClasses,"clsBONFERRONI2.r",sep=""))
# source(paste(pathClasses,"clsBOXPLOT.r",sep=""))
# source(paste(pathClasses,"clsCALCPDIFF4.r",sep=""))
# source(paste(pathClasses,"clsCALCPHET2.r",sep=""))
# source(paste(pathClasses,"clsCALCULATE3.r",sep=""))
# source(paste(pathClasses,"clsCHECKCOLS2.r",sep=""))
# source(paste(pathClasses,"clsCLEAN5.r",sep=""))
# source(paste(pathClasses,"clsCLEANDUPLICATES8.r",sep=""))
# source(paste(pathClasses,"clsCLUMP4.r",sep=""))
# source(paste(pathClasses,"clsCRITERION6.r",sep=""))
# source(paste(pathClasses,"clsEDITCOL.r",sep=""))
# source(paste(pathClasses,"clsEVALSTAT7.r",sep=""))
# source(paste(pathClasses,"clsEXTRACTSNPS3.r",sep=""))
# source(paste(pathClasses,"clsFDR3.r",sep=""))
# source(paste(pathClasses,"clsFILTER6.r",sep=""))
# source(paste(pathClasses,"clsFLIPSTRAND3.r",sep=""))
# source(paste(pathClasses,"clsGC5.r",sep=""))
# source(paste(pathClasses,"clsGETCOLS.r",sep=""))
# source(paste(pathClasses,"clsGETNUM2.r",sep=""))
# source(paste(pathClasses,"clsINDEP4.r",sep=""))
# source(paste(pathClasses,"clsJOINTTEST2.r",sep=""))
# # source(paste(pathClasses,"clsMCMARKER6.r",sep=""))
# # source(paste(pathClasses,"clsMCSTRAND4.r",sep=""))
# source(paste(pathClasses,"clsMERGEEASYIN4.r",sep=""))
# # source(paste(pathClasses,"clsMERGESTRAT.r",sep=""))
# source(paste(pathClasses,"clsMETAANALYSIS3.r",sep=""))
# source(paste(pathClasses,"clsMHPLOT14.r",sep=""))
# source(paste(pathClasses,"clsMIAMIPLOT14.r",sep=""))
# source(paste(pathClasses,"clsPZPLOT7.r",sep=""))
# source(paste(pathClasses,"clsQQPLOT17.r",sep=""))
# # source(paste(pathClasses,"clsRBINDTRAITS.r",sep=""))
# source(paste(pathClasses,"clsREMOVECOL.r",sep=""))
# source(paste(pathClasses,"clsRENAMECOL2.r",sep=""))
# #source(paste(pathClasses,"clsRENAMEMARKER9.r",sep=""))
# source(paste(pathClasses,"clsRENAMEMARKER10.r",sep=""))
# source(paste(pathClasses,"clsRPLOT18.r",sep=""))
# source(paste(pathClasses,"clsSTRSPLITCOL.r",sep=""))
# source(paste(pathClasses,"clsWRITE6.r",sep=""))

# if(!file.exists(fileECF))
	# stop(paste("EASY ERROR:\n ECF-file \n ",fileECF,"\n does not exist!!!\n", sep=""))

# connection.Logfile <- file(paste(fileECF,".out",sep=""),open='w')
# sink(connection.Logfile, split=TRUE)
# ##########################################################################################################################
# ##########################################################################################################################
# ##########################################################################################################################

EasyX.run <- function(fileECF){ 
	
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	cat("|       EasyStrata     |     v8.6     |    20/June/2014      |\n")
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	cat("|  (C) 2013 Thomas Winkler, GNU General Public License, v3   |\n")
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	cat("|  For bug-report, please e-mail:                            |\n")
	cat("|  thomas.winkler@klinik.uni-regensburg.de                   |\n")
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	
	cat("\n")
	cat("+++++\n")
	.timeStart <- Sys.time()
	cat(paste("Starting EasyX:", .timeStart,"\n"))
	cat(paste("Running script ",fileECF, "...\n"))
	
	idxDeviceOffset <- as.integer(dev.cur()) - 1
	
	fileECFName <- strsplit(fileECF,"/")[[1]][length(strsplit(fileECF,"/")[[1]])]
	##########################
	#print("Reading EasyX config-data...")
	#aEcfRows	<- scan(file = fileECF, what=character(0), n = -1, comment.char = "#",  strip.white = TRUE,sep = "\n",quiet=TRUE)
	aEcfRows	<- scan(file = fileECF, what=character(0), n = -1, comment.char = "",  strip.white = TRUE,sep = "\n",quiet=TRUE)
	isRemove = substring(aEcfRows,1,1)=="#"
	aEcfRows = aEcfRows[-which(isRemove)]

	#iScript = match(c("START EASYQC", "STOP EASYQC"), aEcfRows)

	iStart = which(aEcfRows == "START EASYGWA" | aEcfRows == "START EASYQC" | aEcfRows == "START EASYX" | aEcfRows == "START EASYPLOT" | aEcfRows == "START EASYSTRATA")
	iMergeEasyin = which(substring(aEcfRows,1,11) == "MERGEEASYIN")
	iMergeStrat = which(substring(aEcfRows,1,10) == "MERGESTRAT")
	iRbindTrait = which(substring(aEcfRows,1,11) == "RBINDTRAITS")
	iStop = which(aEcfRows == "STOP EASYGWA" | aEcfRows == "STOP EASYQC" | aEcfRows == "STOP EASYX" | aEcfRows == "STOP EASYPLOT" | aEcfRows == "STOP EASYSTRATA")

	#if(any(is.na(iScript)) | iScript[2] < iScript[1]) 
	if(!all(length(iStart)==1 & length(iStop)==1 & iStart < iStop))
		stop(paste("EASY ERROR:\n 'START EASYSTRATA' and 'STOP EASYSTRATA' not properly defined in \n  ",fileECFName,"\n !!!\n", sep=""))
	if(length(iMergeEasyin) > 0) {
		if(!all(length(iMergeEasyin) == 1 & iMergeEasyin > iStart & iMergeEasyin < iStop))
			stop(paste("EASY ERROR:\n 'MERGEEASYIN' not properly defined in \n  ",fileECFName,"\n This command may only be stated once !!!\n", sep=""))
	}
	if(length(iMergeStrat) > 0) {
		if(!all(length(iMergeStrat) == 1 & iMergeStrat > iStart & iMergeStrat < iStop))
			stop(paste("EASY ERROR:\n 'MERGESTRAT' not properly defined in \n  ",fileECFName,"\n This command may only be stated once !!!\n", sep=""))
	}
	if(length(iRbindTrait) > 0) {
		if(!all(length(iRbindTrait) == 1 & iRbindTrait > iStart & iRbindTrait < iStop))
			stop(paste("EASY ERROR:\n 'RBINDTRAITS' not properly defined in \n  ",fileECFName,"\n This command may only be stated once !!!\n", sep=""))
	}
	cat("\n")
	cat("+++++\n")
	cat("Getting list of input files ...\n")
	
	#if(iScript[1] == 1)
	#	stop(paste("EASY ERROR:\n No DEFINE and no EASYIN argument found in configuration part of ECF-file \n ",fileECFName,"\n !!!\n", sep=""))
	
	### Get aEcfConfigCommands
	aindConfigRows = c(1,iStart-1)
	aEcfConfigRows <- aEcfRows[aindConfigRows[1]:aindConfigRows[2]]
	strEcfConfigRows=paste(aEcfConfigRows,collapse="\n")
	strEcfConfigRows=gsub("\t"," ",strEcfConfigRows)
	strEcfConfigRows.ByCommand=gsub("\n--"," --",strEcfConfigRows)
	while(grepl("  ", strEcfConfigRows.ByCommand)) strEcfConfigRows.ByCommand <- sub("  "," ",strEcfConfigRows.ByCommand)
	aEcfConfigCommands=strsplit(strEcfConfigRows.ByCommand,"\n")[[1]]
	
	### Get aEqcScriptCommands
	aindScriptRows = c(iStart+1,iStop-1)
	aEqcScriptRows <- aEcfRows[aindScriptRows[1]:aindScriptRows[2]]
	strEqcScriptRows=paste(aEqcScriptRows,collapse="\n")
	strEqcScriptRows=gsub("\t"," ",strEqcScriptRows)
	strEqcScriptRows.ByCommand=gsub("\n--"," --",strEqcScriptRows)
	while(grepl("  ", strEqcScriptRows.ByCommand)) strEqcScriptRows.ByCommand <- sub("  "," ",strEqcScriptRows.ByCommand)
	aEqcScriptCommands=strsplit(strEqcScriptRows.ByCommand,"\n")[[1]]
	
	#### Check correct usage of functions
	astrEasyFuns.used <- unlist(lapply(strsplit(aEqcScriptCommands," "),function(x) x[1]))
	astrEasyFuns.allowed <- c("MERGE", "MERGEEASYIN","SPLOT", "ADDCOL","ADJUSTALLELES",
								"CALCULATE","CLEAN","CLEANDUPLICATES","CRITERION","EDITCOL","EVALSTAT","EXTRACTSNPS","FILTER",
								"FLIPSTRAND","GC","GETCOLS","GETNUM","QQPLOT","REMOVECOL","RENAMECOL","RENAMEMARKER","RPLOT",
								"STRSPLITCOL","WRITE","AFCHECK","PZPLOT",
								"ANNOTATE","BONFERRONI","CALCPDIFF","CALCPHET","CLUMP","FDR","INDEP","JOINTTEST","METAANALYSIS","MHPLOT","MIAMIPLOT")
	isEasyFunMatch = astrEasyFuns.used%in%astrEasyFuns.allowed
	if(any(!isEasyFunMatch)) 
		stop(paste("EASY ERROR:\n Functions \n ",paste(astrEasyFuns.used[which(!isEasyFunMatch)],collapse=","),"\n are not allowed EasyX functions. Please revise or remove the function !!!\n", sep=""))
	
	iMergeCommand = which(substring(aEqcScriptCommands,1,11) == "MERGEEASYIN")
	nEqcCommands = length(aEqcScriptCommands)
	ls_aEqcScriptCommands = list()
	if(length(iMergeCommand)==1) {
		ls_aEqcScriptCommands = list(aEqcScriptCommands[1:iMergeCommand])
		if(nEqcCommands > iMergeCommand) ls_aEqcScriptCommands = c(ls_aEqcScriptCommands, list(aEqcScriptCommands[(iMergeCommand+1):nEqcCommands]))
	} else {
		ls_aEqcScriptCommands = list(aEqcScriptCommands)
	}
	
	################################################################################################################################################
	################################################################################################################################################
	
	container.GWADATA <- container.EASYMERGE <- list()
	objGWADATA.default 	<- GWADATA()
	objEASYMERGE.default 	<- EASYMERGE()
	isEasyMerge = FALSE
	icount.GWADATA <- 0
	
	for(iConfigCommand in 1:length(aEcfConfigCommands)) {

		strConfigCommand = aEcfConfigCommands[iConfigCommand]
		strConfigFun = strsplit(strConfigCommand," ")[[1]][1]
		
		if(strConfigFun == "DEFINE") {
			objGWADATA.default 			<- GWADATA.define(objGWADATA.default, strConfigCommand)
		}
		if(strConfigFun == "EASYIN") {
			objGWADATA 	<- GWADATA.easyin(objGWADATA.default, strConfigCommand)	
			
			if(objGWADATA@fileInType == "METALSCRIPT") 		ls_objGWADATA <- GWADATA.getmetalscriptfiles(objGWADATA)
			else if(objGWADATA@fileInType == "FILELIST") 	ls_objGWADATA <- GWADATA.getfilelistfiles(objGWADATA)
			else if(grepl("*",objGWADATA@fileIn, fixed=T))	ls_objGWADATA <- GWADATA.getfiles(objGWADATA)
			else ls_objGWADATA <- list(objGWADATA)
			
			container.GWADATA = c(container.GWADATA , ls_objGWADATA)
			icount.GWADATA = icount.GWADATA + 1
		}
		# EASYMERGE
		if(strConfigFun == "EASYMERGE") {
			objEASYMERGE <- EASYMERGE.easymerge(objEASYMERGE.default, strConfigCommand, icount.GWADATA)	
			ls_objEASYMERGE <- list(objEASYMERGE)
			container.EASYMERGE = c(container.EASYMERGE , ls_objEASYMERGE)
			isEasyMerge = TRUE
		}
		
	}
	
	nGWA <- length(container.GWADATA)
	if(isEasyMerge) aiEasyMergeIDs = unlist(lapply(container.EASYMERGE, function(x) x@iMergeID))
	
	if(nGWA < 1)
		stop(paste("EASY ERROR:\n No EASYIN argument found in configuration part of ECF-file \n ",fileECFName,"\n PLease specify at least one input file !!!\n", sep=""))
	
	cat("Using:\n")
	for(i in 1:nGWA) {
		#stop()
		container.GWADATA[[i]] <- GWADATA.init(container.GWADATA[[i]])
		cat(paste("   + ",container.GWADATA[[i]]@fileIn,"\n"))
		if(isEasyMerge) {
			if(any(i == aiEasyMergeIDs)) {
				aiMatch = which(aiEasyMergeIDs == i)
				for(iMatch in aiMatch) {
					container.EASYMERGE[[iMatch]] <- EASYMERGE.init(container.EASYMERGE[[iMatch]])
					cat(paste("   + EASYMERGE ",container.EASYMERGE[[iMatch]]@fileIn,"\n"))
				}
			}
		}
	}
	
	### PathOut set to pathOut from first GWADATA
	pathOut <- ifelse(objGWADATA.default@pathOut != getwd(), objGWADATA.default@pathOut, container.GWADATA[[1]]@pathOut)
	
	#pathOut <- container.GWADATA[[1]]@pathOut
	fileOutBody <- paste(pathOut,"/",fileECFName,sep="")
	fileOutBody <- gsub(".ecf", "", fileOutBody)
	
	cat("\n")
	cat("+++++\n")
	cat("Default output path is \n")
	cat(pathOut)
	cat("\n")
		
	#	return(list(objGWADATA.default, container.GWADATA, fileOutBody))

	#}
	
	
	########################
	################################################################################################################################################
	################################################################################################################################################
	
	cat("\n")
	cat("+++++\n")
	cat("Performing validity check on 10 rows from each file :\n")
	#stop()
	#aEqcScriptCommandsTmp=ls_aEqcScriptCommands[[1]]
	
	container.GWADATA.start <- container.GWADATA
	blnAddFileInTagToReport <- length(unique(unlist(lapply(container.GWADATA,function(x) x@fileInTag)))) > 1
	
	################################################################################################################################################
	################################################################################################################################################
	########################
		
	for(iRun in 1:2) {
		container.GWADATA <- container.GWADATA.start
		
		if(iRun == 1) {
			isValidScript = FALSE
		} else {
			cat("\n Passed validity check!\n")
			isValidScript = TRUE
		}
		
		iCode = 1
		
		for(aEqcScriptCommands in ls_aEqcScriptCommands) {
				
			#EasyX.runseq <- function(aEqcScriptCommands, container.GWADATA, objGWADATA.default, fileOutBody, isValidScript) {
			
			if(iCode == 2 & isValidScript) fileOutBody <- paste(fileOutBody, ".merged",sep="")
			
			objREPORT 		<- REPORT(fileOutBody)
			
			container.ReportPlots = list()
			container.MultiPlots = list()
			afileOutMultiPlot = c()
			container.Merge = list()
			container.Boxplots = list() ## [[1]] list of boxplotstat of boxplot 1
			container.Rename = list()
			container.GWADATA.out = list()
			
			nGWA <- length(container.GWADATA)
			
			for(iGWA in 1:nGWA) {
				
				isGarbageCleaning <- ifelse(isValidScript, objGWA@blnGarbageCleaning, FALSE)
				
				iReportPlot <- iMultiPlot <- iMerge <- iBoxplot <- iRename <- 0
				
				objGWA <- container.GWADATA[[iGWA]]
				if(isValidScript) {
					cat("\n")
					cat("+++++\n")
					cat(paste("Processing file:",objGWA@fileInShortName,"\n"))
				} else {
					cat("\n")
					cat(paste("   + ",objGWA@fileInShortName,"-> "))
				}
				
				#if(!(objGWA@blnMergedStrat | objGWA@blnRbindTraits)) {
				if(!(objGWA@blnMergedEasyin)) {
					if(isValidScript) 	objGWA <- GWADATA.read(objGWA)
					else				objGWA <- GWADATA.read.10rows(objGWA)
				}
				
				### EASYMERGE
				if(isEasyMerge & iCode == 1) {
					if(any(iGWA == aiEasyMergeIDs)) {
						aiMatch = which(aiEasyMergeIDs == iGWA)
						for(iMatch in aiMatch) {
							objEASYMERGE <- container.EASYMERGE[[iMatch]]
							if(isValidScript) 	objEASYMERGE <- EASYMERGE.read(objEASYMERGE)
							else				objEASYMERGE <- EASYMERGE.read.10rows(objEASYMERGE)
							EASYMERGE.GWADATA.valid(objEASYMERGE, objGWA)
							objGWA <- EASYMERGE.run(objEASYMERGE, objGWA)
						}
					}
				}
				###
				objREPORT	<-	REPORT.addval(objREPORT,"fileInShortName",objGWA@fileInShortName)
				if(blnAddFileInTagToReport) objREPORT	<-	REPORT.addval(objREPORT,"fileInTag",objGWA@fileInTag)
				#objREPORT	<-	REPORT.addval(objREPORT,"fileInTag",objGWA@fileInTag)
				objREPORT	<-	REPORT.addval(objREPORT,"numSNPsIn",dim(objGWA@tblGWA)[1])
				objREPORT	<-	REPORT.addval(objREPORT,"numSNPsOut","NA")
				
				#### Go through script commands
				
				for(iCommand in 1:length(aEqcScriptCommands)) {
				
					strCommand = aEqcScriptCommands[iCommand]
					strScriptFun = strsplit(strCommand," ")[[1]][1]
					
					if(isValidScript) {
						cat(paste("   + ",gsub("--","\n      --",strCommand),sep=""))
						cat("\n")
					}
					
					if(strScriptFun == "MIAMIPLOT") {
					
						objMIAMIPLOT <- MIAMIPLOT(strCommand)
						MIAMIPLOT.GWADATA.valid(objMIAMIPLOT, objGWA)
						
						if(isValidScript) {
							#fnOpenPng.Miami(objGWA, "miami")					
							fnOpenPlot(objGWA, objMIAMIPLOT)
							MIAMIPLOT.run(objMIAMIPLOT, objGWA)
							fnClosePlot()
						}
					}
					if(strScriptFun == "MHPLOT") {
						objMHPLOT <- MHPLOT(strCommand)
						MHPLOT.GWADATA.valid(objMHPLOT, objGWA)
						
						if(isValidScript) {
							#fnOpenPng.MH(objGWA, "mh")		
							fnOpenPlot(objGWA, objMHPLOT)
							MHPLOT.run(objMHPLOT, objGWA)						
							fnClosePlot()
						}
					}
					if(strScriptFun == "ANNOTATE") {
						objANNOT 	<- ANNOTATE(strCommand)
						ANNOTATE.GWADATA.valid(objANNOT, objGWA)		
						objGWA 	<- ANNOTATE.run(objANNOT, objGWA)							
					}
					if(strScriptFun == "INDEP") {
						objINDEP 	<- INDEP(strCommand)
						INDEP.GWADATA.valid(objINDEP, objGWA)							
						#if(isValidScript) stop()
						lsOut 	<- INDEP.run(objINDEP, objGWA, objREPORT)
						
						objGWA		 	<- lsOut[[1]]
						objGWA.indep 	<- lsOut[[2]]
						objGWA.indep.x 	<- lsOut[[3]]
						objREPORT 		<- lsOut[[4]]
						rm(lsOut)
						if(isValidScript) {
							REPORT.write(objREPORT)
							strTagIndep <- ifelse(objINDEP@strTag != "", paste(".",objINDEP@strTag,sep=""), objINDEP@strTag)
							GWADATA.write(objGWA.indep, strSuffix = paste(strTagIndep,".indep",sep=""))
							GWADATA.write(objGWA.indep.x, strSuffix = paste(strTagIndep,".indepX",sep=""))
							rm(objGWA.indep)
							rm(objGWA.indep.x)
							rm(strTagIndep)
						}
					}
					if(strScriptFun == "CLUMP") {
						objCLUMP 	<- CLUMP(strCommand)
						CLUMP.GWADATA.valid(objCLUMP, objGWA)							
						#if(isValidScript) stop()
						if(isValidScript) {
							lsOut 	<- CLUMP.run(objCLUMP, objGWA, objREPORT)
							objGWA		 			<- lsOut[[1]]
							objGWA.clump 			<- lsOut[[2]]
							objGWA.clump.x 			<- lsOut[[3]]
							objGWA.clump.nomatch 	<- lsOut[[4]]
							objREPORT 				<- lsOut[[5]]
							rm(lsOut)
							REPORT.write(objREPORT)
							strTagClump <- ifelse(objCLUMP@strTag != "", paste(".",objCLUMP@strTag,sep=""), objCLUMP@strTag)
							GWADATA.write(objGWA.clump, strSuffix = paste(strTagClump,".clump",sep=""))
							GWADATA.write(objGWA.clump.x, strSuffix = paste(strTagClump,".clumpX",sep=""))
							if(nrow(objGWA.clump.nomatch@tblGWA)>0) {
								GWADATA.write(objGWA.clump.nomatch, strSuffix = paste(strTagClump,".clump_nomatch",sep=""))
							}
							rm(objGWA.clump)
							rm(objGWA.clump.x)
							rm(objGWA.clump.nomatch)
							rm(strTagClump)
						}
					}
					if(strScriptFun == "BONFERRONI") {
						objBF 	<- BONFERRONI(strCommand)
						lsOut 	<- BONFERRONI.run(objBF, objGWA, objREPORT,isValidScript)
						objGWA 		<- lsOut[[1]]
						objREPORT 	<- lsOut[[2]]
						rm(lsOut)
						if(isValidScript) {
							REPORT.write(objREPORT)
						}
					}
					if(strScriptFun == "FDR") {
						objFDR 	<- FDR(strCommand)				
						objGWA 	<- FDR.run(objFDR, objGWA)
					}
					if(strScriptFun == "JOINTTEST") {
						objJT 	<- JOINTTEST(strCommand)				
						objGWA 	<- JOINTTEST.run(objJT, objGWA)
					}
					if(strScriptFun == "CALCPDIFF") {
						objPD 	<- CALCPDIFF(strCommand)				
						lsOut 	<- CALCPDIFF.run(objPD, objGWA, objREPORT)
						objGWA 		<- lsOut[[1]]
						objREPORT 	<- lsOut[[2]]
						rm(lsOut)
						gc()
						if(isValidScript) REPORT.write(objREPORT)
					}
					if(strScriptFun == "CALCPHET") {
						objCP 	<- CALCPHET(strCommand)				
						objGWA 	<- CALCPHET.run(objCP, objGWA)
					}
					if(strScriptFun == "METAANALYSIS") {
						objMA 	<- METAANALYSIS(strCommand)				
						objGWA 	<- METAANALYSIS.run(objMA, objGWA)
					}
					if(strScriptFun == "MERGEEASYIN") {
						
						fileOutShortName = sub(".ecf","",fileECFName)
						if(objGWA@fileInMergeTag != "") 
							fileOutShortName = paste(fileOutShortName, objGWA@fileInMergeTag, sep=".")
						
						objME <- MERGEEASYIN(strCommand, fileOutShortName)
						MERGEEASYIN.GWADATA.valid(objME, objGWA)
						if(iGWA == 1) {
							### start first
							objGWA.merged <- MERGEEASYIN.start(objME, objGWADATA.default, objGWA)
						}
						if(iGWA > 1 & objGWA.merged@fileInMergeTag == objGWA@fileInMergeTag) {
							### add 
							objGWA.merged <- MERGEEASYIN.run(objME, objGWA.merged, objGWA)
						}
						if(iGWA > 1 & objGWA.merged@fileInMergeTag != objGWA@fileInMergeTag) {
							### start over
							container.GWADATA.out <- c(container.GWADATA.out, objGWA.merged)
							objGWA.merged <- MERGEEASYIN.start(objME, objGWADATA.default, objGWA)
						} 
						if(iGWA == nGWA) {
							### add last pair to output
							container.GWADATA.out <- c(container.GWADATA.out, objGWA.merged)
						}
					}
					# if(strScriptFun == "RENAMEMARKER") {
						
						# iRename = iRename + 1
						# #if(iGWA == 1) {
						# if(iRename>length(container.Rename)) {
							
							# objRENAMEMARKER <- RENAMEMARKER(strCommand)
							# objRENAMEMARKER <- RENAMEMARKER.read(objRENAMEMARKER, blnReadAll = isValidScript)
							# container.Rename[[iRename]] <- objRENAMEMARKER
						# }
							
						# objRENAMEMARKER <- container.Rename[[iRename]]
						# RENAMEMARKER.valid(objRENAMEMARKER)
						# RENAMEMARKER.GWADATA.valid(objRENAMEMARKER, objGWA)
						# ## 
										
						# lsOut <- RENAMEMARKER.run(objRENAMEMARKER, objGWA, objREPORT, isValidScript)

						# objGWA 		<- lsOut[[1]]
						# objREPORT 	<- lsOut[[2]]
						# rm(lsOut)
						# rm(objRENAMEMARKER)
					# }
					# if(strScriptFun == "AFCHECK") {
						# objAC <- AFCHECK(strCommand)
						# ## objAC inherits from classes SPLOT and from ADJUSTALLELES (+ MERGE)
						# if(objAC@blnAAMerge) {
							# iMerge = iMerge + 1
							# if(iMerge>length(container.Merge)) {
							# #if(iGWA == 1) {
								# objAC <- MERGE.init(objAC) 
							
								# if(isValidScript) objAC <- GWADATA.read(objAC)
								# else objAC <- GWADATA.read.10rows(objAC)
								
								# container.Merge[[iMerge]] <- objAC
							# }
							# objAC <- container.Merge[[iMerge]]
							# MERGE.GWADATA.valid(objAC, objGWA)
							
							# lsOut <- MERGE.run(objAC, objGWA, objREPORT, isValidScript) ### WILL JUST BE USED FOR THE PLOT!
							# objGWA <- lsOut[[1]]
							# objREPORT <- lsOut[[2]]
							# rm(lsOut)
							# if(isValidScript) REPORT.write(objREPORT)
						# } 
						# if(objAC@blnAdjAllele) {
							# ADJUSTALLELES.GWADATA.valid(objAC, objGWA)
							# lsOut <- ADJUSTALLELES.run(objAC, objGWA, objREPORT, isValidScript)
							# objGWA 			<- lsOut[[1]]
							# objREPORT 		<- lsOut[[2]]
							# # objGWA.miss 	<- lsOut[[3]]
							# # objGWA.invalid 	<- lsOut[[4]]
							# rm(lsOut)
						# }
						# SPLOT.GWADATA.valid(objAC, objGWA)
						
						# if(isValidScript) {
							# # if(objAC@strMode == "singleplot") {
								# # fnOpenPng(objGWA, "sp")
								# # SPLOT.run(objAC, objGWA)
								# # fnClosePlot()
							# # } 
							# if(objAC@strMode == "subplot") {
								# iMultiPlot = iMultiPlot + 1
								# if(iMultiPlot > length(container.MultiPlots)) {
									# #fnOpenMultiPng(fileOutBody, paste(iMultiPlot,".sp",sep=""), nGWA)
									# fileOutMultiPlot = fnOpenMultiPlot(fileOutBody, objAC, nGWA, afileOutMultiPlot)
									# afileOutMultiPlot = c(afileOutMultiPlot, fileOutMultiPlot)
									# SPLOT.run(objAC, objGWA)
									# container.MultiPlots[[iMultiPlot]]  <- objAC
								# } else {
									## fnAddPlot(iMultiPlot)
									# fnAddPlot(iMultiPlot + idxDeviceOffset)
									# objAC <- container.MultiPlots[[iMultiPlot]]
									# SPLOT.run(objAC, objGWA)
								# }
								# if(iGWA == nGWA) fnClosePlot()
							# }
							# lsOut <- AFCHECK.run(objAC, objGWA, objREPORT)
							# #objGWA.outlier <- lsOut[[1]]
							# objGWA 		<- lsOut[[1]]
							# objREPORT 	<- lsOut[[2]] 
							# rm(lsOut)
							# # if(objAC@blnWriteOutlier & dim(objGWA.outlier@tblGWA)[1] > 0) {
								# # GWADATA.write(objGWA.outlier, strSuffix = paste(".",objAC@strTag,".outlier",sep=""))
								# # rm(objGWA.outlier)
							# # }
							# REPORT.write(objREPORT)
							# #rm(objGWA.tmp)
							# rm(objAC)
						# }
						
					# }
					if(strScriptFun == "ADJUSTALLELES") {
						objAA <- ADJUSTALLELES(strCommand)
						
						if(objAA@blnAAMerge) {
							iMerge = iMerge + 1
							if(iMerge>length(container.Merge)) {
							#if(iGWA == 1) {
								objAA <- MERGE.init(objAA) 
							
								if(isValidScript) objAA <- GWADATA.read(objAA)
								else objAA <- GWADATA.read.10rows(objAA)
								
								container.Merge[[iMerge]] <- objAA
							}
							objAA <- container.Merge[[iMerge]]
							MERGE.GWADATA.valid(objAA, objGWA)
							lsOut <- MERGE.run(objAA, objGWA,objREPORT,isValidScript)
							objGWA 				<- lsOut[[1]]
							objREPORT 			<- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) REPORT.write(objREPORT)
						}
						ADJUSTALLELES.GWADATA.valid(objAA, objGWA)
						lsOut <- ADJUSTALLELES.run(objAA, objGWA, objREPORT, isValidScript)
						
						objGWA 			<- lsOut[[1]]
						objREPORT 		<- lsOut[[2]]
						if(isValidScript) REPORT.write(objREPORT)
						
						# objGWA.miss 	<- lsOut[[3]]
						# objGWA.invalid 	<- lsOut[[4]]
						rm(lsOut)
						rm(objAA)
					}
					if(strScriptFun == "MERGE") {
						iMerge = iMerge + 1
						if(iMerge>length(container.Merge)) {
						# if(iGWA == 1) {
							objMERGE <- MERGE(strCommand)
							if(isValidScript) objMERGE <- GWADATA.read(objMERGE)
							else objMERGE <- GWADATA.read.10rows(objMERGE)
							
							container.Merge[[iMerge]] <- objMERGE
						}
						objMERGE <- container.Merge[[iMerge]]
						if(objMERGE@fileRefTag == "1" | objMERGE@fileRefTag == objGWA@fileInTag) {
							MERGE.GWADATA.valid(objMERGE, objGWA)
							lsOut <- MERGE.run(objMERGE, objGWA, objREPORT, isValidScript)
							objGWA 		<- lsOut[[1]]
							objREPORT 	<- lsOut[[2]]
							rm(lsOut)
							if(isValidScript) REPORT.write(objREPORT)
						}
						rm(objMERGE)
					}
					if(strScriptFun == "WRITE") {
						objWRITE 	<- WRITE(strCommand)			
						if(isValidScript) {
							objREPORT <- WRITE.run(objWRITE, objGWA, objREPORT)
							REPORT.write(objREPORT)
						}
					}
					if(strScriptFun == "ADDCOL") {
						objADDCOL 	<- ADDCOL(strCommand)				
						objGWA 		<- ADDCOL.run(objADDCOL, objGWA)
					}
					if(strScriptFun == "EDITCOL") {
						objEDITCOL 	<- EDITCOL(strCommand)				
						EDITCOL.GWADATA.valid(objEDITCOL, objGWA)
						objGWA 		<- EDITCOL.run(objEDITCOL, objGWA)
					}
					if(strScriptFun == "STRSPLITCOL") {
						objSTRSPLITCOL 	<- STRSPLITCOL(strCommand)				
						STRSPLITCOL.GWADATA.valid(objSTRSPLITCOL, objGWA)
						objGWA 		<- STRSPLITCOL.run(objSTRSPLITCOL, objGWA)
					}
					if(strScriptFun == "EXTRACTSNPS") {
						objES 	<- EXTRACTSNPS(strCommand)
						EXTRACTSNPS.GWADATA.valid(objES, objGWA)
						lsOut 	<- EXTRACTSNPS.run(objES, objGWA, objREPORT)
						tblOut 		<- lsOut[[1]]
						objREPORT 	<- lsOut[[2]]
						rm(lsOut)
						if(isValidScript) {				
							REPORT.write(objREPORT)
							strTagExtract <- ifelse(objES@strTag != "", paste(".",objES@strTag,sep=""), objES@strTag)
							if(strsplit(objGWA@pathOut,"")[[1]][nchar(objGWA@pathOut)] != "/") pathOutExtract <- paste(objGWA@pathOut,"/",sep="")
							else pathOutExtract<-objGWA@pathOut
							fileExtractOut = paste(pathOutExtract,objGWA@fileInShortName,strTagExtract,".extracted",sep="")
							write.table(tblOut, fileExtractOut, row.names=F, quote=F, sep="\t", na="NA")
							rm(tblOut)
							rm(strTagExtract)
							rm(pathOutExtract)
							rm(fileExtractOut)
						}
					}
					# if(strScriptFun == "FLIPSTRAND") {
						# objFS 	<- FLIPSTRAND(strCommand)
						# FLIPSTRAND.GWADATA.valid(objFS, objGWA)
						# lsOut 	<- FLIPSTRAND.run(objFS, objGWA, objREPORT)
						# objGWA 		<- lsOut[[1]]
						# objREPORT 	<- lsOut[[2]]
						# rm(lsOut)
						# if(isValidScript) REPORT.write(objREPORT)
					# }
					if(strScriptFun == "EVALSTAT") {
						objEVALSTAT <- EVALSTAT(strCommand)
						objREPORT 	<- EVALSTAT.run(objEVALSTAT, objGWA, objREPORT)
						if(isValidScript) REPORT.write(objREPORT)
					}
					if(strScriptFun == "GETNUM") {
						objGN 		<- GETNUM(strCommand)
						objREPORT 	<- GETNUM.run(objGN, objGWA, objREPORT)
						if(isValidScript) REPORT.write(objREPORT)
					}
					if(strScriptFun == "CALCULATE") {
						objCALC 	<- CALCULATE(strCommand)				
						objREPORT 	<- CALCULATE.run(objCALC, objGWA, objREPORT)
						if(isValidScript) REPORT.write(objREPORT)
					}
					if(strScriptFun == "FILTER") {
						objFILTER 	<- FILTER(strCommand)
						lsOut 		<- FILTER.run(objFILTER, objGWA, objREPORT)
						objGWA 		<- lsOut[[1]]
						objREPORT 	<- lsOut[[2]]
						rm(lsOut)
						if(isValidScript) REPORT.write(objREPORT)
					}
					if(strScriptFun == "CLEAN") {
						objCLEAN 	<- CLEAN(strCommand)
						lsOut 		<- CLEAN.run(objCLEAN, objGWA, objREPORT, !isValidScript)
						objGWA 		<- lsOut[[1]]
						objREPORT 	<- lsOut[[2]]
						objGWA.Cleaned 	<- lsOut[[3]]
						rm(lsOut)
						if(isValidScript) {
							REPORT.write(objREPORT)
							if(objCLEAN@blnWriteCleaned&nrow(objGWA.Cleaned@tblGWA)>0) 
								GWADATA.write(objGWA.Cleaned, strSuffix = paste(".",objCLEAN@strCleanName,sep=""))
							rm(objGWA.Cleaned)
						}
					}
					if(strScriptFun == "CRITERION") {
						objCRIT 	<- CRITERION(strCommand)
						lsOut 	<- CRITERION.run(objCRIT, objGWA, objREPORT)
						
						objGWA.Crit <- lsOut[[1]]
						objREPORT <- lsOut[[2]]
						rm(lsOut)
						if(isValidScript) {
							REPORT.write(objREPORT)
							#GWADATA.write(objGWA.Crit, strSuffix = ".criterion")
							if(objCRIT@blnWriteEmpty | nrow(objGWA.Crit@tblGWA)>0) 
								GWADATA.write(objGWA.Crit, strSuffix = paste(".",sub("numSNP_","",objCRIT@strCritName),sep=""))
							rm(objGWA.Crit)
						}
					}
					if(strScriptFun == "GETCOLS") {
						objGC 	<- GETCOLS(strCommand)
						objGWA 	<- GETCOLS.run(objGC, objGWA)
					}
					if(strScriptFun == "RENAMECOL") {
						objRC 		<- RENAMECOL(strCommand)
						lsOut 		<- RENAMECOL.run(objRC, objGWA, objREPORT)
						objGWA 		<- lsOut[[1]]
						objREPORT 	<- lsOut[[2]]
						rm(lsOut)
						if(isValidScript) REPORT.write(objREPORT)
					}
					if(strScriptFun == "REMOVECOL") {
						objRC 	<- REMOVECOL(strCommand)
						objGWA 	<- REMOVECOL.run(objRC, objGWA)
					}
					# if(strScriptFun == "CLEANDUPLICATES") {
						# objCD 		<- CLEANDUPLICATES(strCommand)
						# CLEANDUPLICATES.GWADATA.valid(objCD, objGWA)
						# if(isValidScript){
							# lsOut 		<- CLEANDUPLICATES.run(objCD, objGWA, objREPORT)
							# objGWA 		<- lsOut[[1]]
							# objREPORT 	<- lsOut[[2]]
							# rm(lsOut)
							# REPORT.write(objREPORT)
						# }
					# }
					if(strScriptFun == "GC") {
						objGC 	<- GC(strCommand)
						GC.GWADATA.valid(objGC, objGWA)
						lsOut 	<- GC.run(objGC, objGWA, objREPORT)
						objGWA 		<- lsOut[[1]]
						objREPORT 	<- lsOut[[2]]
						rm(lsOut)
						if(isValidScript) REPORT.write(objREPORT)
					}
					if(strScriptFun == "QQPLOT") {
						
						objQQPLOT <- QQPLOT(strCommand)
						QQPLOT.GWADATA.valid(objQQPLOT, objGWA)
						
						if(isValidScript) {
							if(objQQPLOT@strMode == "singleplot") {
								fnOpenPlot(objGWA, objQQPLOT)
								QQPLOT.run(objQQPLOT, objGWA)
								fnClosePlot()
							} 
							if(objQQPLOT@strMode == "subplot") {
								iMultiPlot = iMultiPlot + 1
								if(iMultiPlot > length(container.MultiPlots)) {
									#fnOpenMultiPng(fileOutBody, paste(iMultiPlot,".sp",sep=""), nGWA)
									fileOutMultiPlot = fnOpenMultiPlot(fileOutBody, objQQPLOT, nGWA, afileOutMultiPlot)
									afileOutMultiPlot = c(afileOutMultiPlot, fileOutMultiPlot)
									QQPLOT.run(objQQPLOT, objGWA)
									container.MultiPlots[[iMultiPlot]]  <- objQQPLOT
								} else {
									#fnAddPlot(iMultiPlot)
									fnAddPlot(iMultiPlot + idxDeviceOffset)
									objQQPLOT <- container.MultiPlots[[iMultiPlot]]
									QQPLOT.run(objQQPLOT, objGWA)
								}
								if(iGWA == nGWA) fnClosePlot()
							}
						}
					}
					if(strScriptFun == "SPLOT") {
						
						objSPLOT <- SPLOT(strCommand)
						SPLOT.GWADATA.valid(objSPLOT, objGWA)
						
						if(isValidScript) {
							if(objSPLOT@strMode == "singleplot") {
								fnOpenPlot(objGWA, objSPLOT)
								SPLOT.run(objSPLOT, objGWA)
								fnClosePlot()
							} 
							if(objSPLOT@strMode == "subplot") {
								iMultiPlot = iMultiPlot + 1
								if(iMultiPlot > length(container.MultiPlots)) {
									#fnOpenMultiPng(fileOutBody, paste(iMultiPlot,".sp",sep=""), nGWA)
									fileOutMultiPlot = fnOpenMultiPlot(fileOutBody, objSPLOT, nGWA, afileOutMultiPlot)
									afileOutMultiPlot = c(afileOutMultiPlot, fileOutMultiPlot)
									SPLOT.run(objSPLOT, objGWA)
									container.MultiPlots[[iMultiPlot]]  <- objSPLOT
								} else {
									# fnAddPlot(iMultiPlot)
									fnAddPlot(iMultiPlot + idxDeviceOffset)
									objSPLOT <- container.MultiPlots[[iMultiPlot]]
									SPLOT.run(objSPLOT, objGWA)
								}
								if(iGWA == nGWA) fnClosePlot()
							}
						}
					}
					# if(strScriptFun == "PZPLOT") {
						# objPZPLOT <- PZPLOT(strCommand)
						# PZPLOT.GWADATA.valid(objPZPLOT, objGWA)
						
						# if(isValidScript) {
							# if(objPZPLOT@strMode == "singleplot") {
								# #fnOpenPng(objGWA, "pz")
								# #fnOpenPlot(objGWA, objPZPLOT, strSuffix = "pz")
								# fnOpenPlot(objGWA, objPZPLOT)
								# SPLOT.run(objPZPLOT, objGWA)
								# fnClosePlot()
							# } 
							# if(objPZPLOT@strMode == "subplot") {
								# iMultiPlot = iMultiPlot + 1
								# #if(iGWA == 1) {
								# if(iMultiPlot > length(container.MultiPlots)) {
									# ## set up new container
									# #fnOpenMultiPng(fileOutBody, paste(iMultiPlot,".pz",sep=""), nGWA)
									# fileOutMultiPlot = fnOpenMultiPlot(fileOutBody, objPZPLOT, nGWA, afileOutMultiPlot)
									# afileOutMultiPlot = c(afileOutMultiPlot, fileOutMultiPlot)
									# SPLOT.run(objPZPLOT, objGWA)
									# container.MultiPlots[[iMultiPlot]]  <- objPZPLOT
								# } else { 
									# ## use existing container
									## fnAddPlot(iMultiPlot)
									# fnAddPlot(iMultiPlot + idxDeviceOffset)
									# objPZPLOT <- container.MultiPlots[[iMultiPlot]]
									# SPLOT.run(objPZPLOT, objGWA)
								# }
								# if(iGWA == nGWA) fnClosePlot()
							# }
						# }
					# }
					if(strScriptFun == "RPLOT") {
							
						#if(iGWA == nGWA) {
						objRPLOT <- RPLOT(strCommand)
						RPLOT.REPORT.valid(objRPLOT, objREPORT)
						
						if(isValidScript) {
							iReportPlot = iReportPlot + 1
							container.ReportPlots[[iReportPlot]]  <- objRPLOT
							# fnOpenPlot(objREPORT, objRPLOT)
							# RPLOT.run(objRPLOT, objREPORT)
							# fnClosePlot()
						}
						#}
					}
					
					if(isGarbageCleaning) gc()
					
# if(isValidScript & iGWA==2) objGWA@tblGWA<-data.frame()					
					if(dim(objGWA@tblGWA)[1]==0) {
						if(isValidScript) {
						cat(paste("EASY WARNING:",strScriptFun,"\n After evaluation of \n",strScriptFun,"\n the GWA data set is empty!!!\n", sep=""))
						cat(paste("Skipping EasyX evaluation for \n",objGWA@fileInShortName,"!!!\n", sep=""))
						}
						break
					}
				}

				if(REPORT.getval(objREPORT, "numSNPsOut")=="NA") objREPORT	<-	REPORT.setval(objREPORT,"numSNPsOut",dim(objGWA@tblGWA)[1])
				if(isValidScript) REPORT.write(objREPORT)
				
				objREPORT 	<- 	REPORT.newrow(objREPORT)
				
				if(!isValidScript) {
					cat("OK")
				} 
			}
			if(strScriptFun == "MERGEEASYIN") container.GWADATA <- container.GWADATA.out
			else container.GWADATA <- list()
			
			## create collected RPLOTS:
			if(length(container.ReportPlots)>0) {
				for(iRP in 1:length(container.ReportPlots)) {
					objRP <- container.ReportPlots[[iRP]]
					fnOpenPlot(objREPORT, objRP)
					RPLOT.run(objRP, objREPORT)
					fnClosePlot()
				}
			}

			iCode = iCode + 1
		}
	} ## i
	
	########################
	################################################################################################################################################
	################################################################################################################################################
	
	cat("\n")
	cat("\n")
	cat("+++++\n")
	.timeStop <- Sys.time()
	cat(paste("Succesfully finished EasyX:", .timeStop,"\n"))
	cat(paste("Elapsed time:", .timeStop-.timeStart,"\n"))
	cat("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n")
	
	#graphics.off()
	#closeAllConnections()
	
	#	Copy log file to output folder
	if(paste(fileECF,".out",sep="") != paste(pathOut,"/",fileECFName,".out",sep="")) {
		file.copy(paste(fileECF,".out",sep=""), paste(pathOut,"/",fileECFName,".out",sep=""),overwrite=TRUE)
	}
	
	################################################################################################################################################
	################################################################################################################################################	
}
EasyX <- function(fileECF) { 
	
	### Wrapper for EasyX.run
	valout <-TRUE
	
	if(!file.exists(fileECF))
		stop(paste("EASY ERROR:\n ECF-file \n ",fileECF,"\n does not exist!!!\n", sep=""))
	
	#graphics.off()
	#closeAllConnections()
	options(show.error.messages = FALSE)
	options(warn=0)
	
	## get current device and connections
	idxDeviceIn <- dev.cur()
	tblConIn 	<- showConnections()
	
	fileLog <- paste(fileECF,".out",sep="")
	if(file.exists(fileLog)) file.remove(fileLog)
	
	#connection.Logfile <- file(fileLog,open='w')
	connection.Logfile <- file(fileLog,open='a')
	sink(connection.Logfile, split=TRUE)

	out <-  NULL
	out <- 	try(withCallingHandlers(
					EasyX.run(fileECF), 
					warning=function(w) { cat(w$message); cat("\n") })
				)
	#out <- try(withCallingHandlers(EasyX.run.warning(fileECF,fileOutLog), warning=function(w) {cat(w$message); cat("\n")}))
	#withCallingHandlers({out <- try(EasyX.run.warning(fileECF,fileOutLog))}, warning=function(w) {cat(w$message); cat("\n")})
	
	if(class(out) == "try-error") {
		cat(out)
		valout <- FALSE
	}
	cat("\n")
	options(show.error.messages = TRUE)
	
	## close all open devices except those that were open at the start
	idxDeviceOut <- dev.cur()
	while(idxDeviceOut!=idxDeviceIn) {
		dev.off()
		idxDeviceOut <- dev.cur()
	}
	## close all open connections except those that were open at the start
	sink()
	close(connection.Logfile)
	# tblConOut 	<- showConnections()
	# while(nrow(tblConOut)>nrow(tblConIn)) {
		# sink()
		# tblConOut <- showConnections()
	# }
	#graphics.off()
	#closeAllConnections()
	return(valout)
}

EasyStrata <- function(fileECF) { 
	return(EasyX(fileECF))
}
EasyQC <- function(fileECF) { 
	return(EasyX(fileECF))
}


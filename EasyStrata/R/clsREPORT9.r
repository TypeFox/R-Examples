setClass("REPORT",
	representation = representation(
						fileReportBody				=	"character",
						fileReport					=	"character",
						tblReport					=	"data.frame",
						iCurrentRow					=	"numeric",
						iCurrentCol					=	"numeric"
						),
	prototype = prototype(
						fileReportBody				=	"",
						fileReport					=	"",
						tblReport					=	data.frame(),
						iCurrentRow					=	1,
						iCurrentCol					=	1
						)
)

setGeneric("setREPORT", function(object) standardGeneric("setREPORT"))
setMethod("setREPORT", signature = (object = "REPORT"), function(object) {
	
	#if(file.exists(object@fileReport)) file.remove(object@fileReport)

	fileReport	<- paste(object@fileReportBody,".rep",sep="") 

	i = 1
	while(file.exists(fileReport)) {
		fileReport <- paste(object@fileReportBody,".",i,".rep",sep="") 
		i = i + 1
	}
	
	object@fileReport <- fileReport
	
	return(object)
})

#############################################################################################################################
REPORT.newrow <- function(objREPORT) {
	
	objREPORT@iCurrentRow <- objREPORT@iCurrentRow + 1
	
	aNewRow = rep("NA", ncol(objREPORT@tblReport))
	objREPORT@tblReport = rbind(objREPORT@tblReport, aNewRow)
	
	objREPORT@iCurrentCol = 1
	
	return(objREPORT)
}

##### Method to add values to Report
setGeneric("REPORT.addval", function(obj1,obj2,obj3, ...) standardGeneric("REPORT.addval"))
## ... from given Table (1 row)
setMethod("REPORT.addval", signature(obj1 = "REPORT", obj2 = "data.frame", obj3="missing"), function(obj1,obj2,obj3, ...) {
	
	objREPORT 	<- obj1
	tblAdd 		<- obj2
	
	if(ncol(tblAdd) > 0) {
		for(i in 1:ncol(tblAdd)) tblAdd[,i] <- as.character(tblAdd[,i])
		
		if(objREPORT@iCurrentRow == 1 & objREPORT@iCurrentCol == 1) {
			### Start new report
			objREPORT@tblReport 	<- tblAdd
			
		} else if(objREPORT@iCurrentRow == 1 & objREPORT@iCurrentCol > 1) {
			### Add to first row
			## Rename column name if it is already present in report table
			# # colNamesAdd <- names(tblAdd)
			# # isAlreadyInTable = colNamesAdd%in%names(objREPORT@tblReport)
			# # iCount = 1
			# # while(any(isAlreadyInTable)) {
				# # names(tblAdd)[isAlreadyInTable] = paste(colNamesAdd[isAlreadyInTable],".",iCount,sep="")
				# # isAlreadyInTable = names(tblAdd)%in%names(objREPORT@tblReport)
				# # iCount = iCount + 1
			# # }
			## Add to first row
			objREPORT@tblReport 	<- cbind(objREPORT@tblReport, tblAdd, stringsAsFactors = FALSE)
			
		} else {
			### Add to currentRow at position currentCol
			# iMatch = match(names(tblAdd), names(objREPORT@tblReport))
			# objREPORT@tblReport[objREPORT@iCurrentRow,iMatch] <- tblAdd[1,]
			
			iRow <- objREPORT@iCurrentRow
			aiCol <- objREPORT@iCurrentCol:(objREPORT@iCurrentCol+ncol(tblAdd)-1)
			while(ncol(objREPORT@tblReport)<aiCol[length(aiCol)]) {
				objREPORT@tblReport <- cbind(objREPORT@tblReport,matrix(NA,nrow(objREPORT@tblReport),ncol(tblAdd)))
				names(objREPORT@tblReport)[aiCol] <- names(tblAdd)
			}
			objREPORT@tblReport[iRow,aiCol] <- tblAdd[1,]
		}
		
		#write.table(objREPORT@tblReport,objREPORT@fileReport,row.names=F,quote=F,sep="\t")
		
		objREPORT@iCurrentCol <- objREPORT@iCurrentCol + ncol(tblAdd)
	}
	
	isDuplCol = duplicated(names(objREPORT@tblReport))
	if(any(isDuplCol)) {
		strDuplCol = names(objREPORT@tblReport)[which(isDuplCol)]
		stop(paste(" EASY ERROR:REPORT\n Created REPORT column \n",paste(strDuplCol,collapse=","), "\n is duplicated in the REPORT. \n PLease use a different name or (for MERGE,ADJUSTALLELES,AFCHECK,GC or RENAMEMARKER) influence the naming by setting --strTag!", sep=""))
	}
	
	return(objREPORT)
})
## from given single character

setMethod("REPORT.addval", signature(obj1 = "REPORT", obj2 = "character",obj3="character"), function(obj1,obj2,obj3, ...) {
	
	objREPORT 	<- obj1
	strValName 	<- obj2
	strVal	 	<- obj3
	
	if(objREPORT@iCurrentRow == 1 & objREPORT@iCurrentCol == 1) {
		### Start new report
		objREPORT@tblReport 	<- data.frame(NewCol=strVal, stringsAsFactors = FALSE)
		names(objREPORT@tblReport)[ncol(objREPORT@tblReport)] <- strValName
		
	} else if(objREPORT@iCurrentRow == 1 & objREPORT@iCurrentCol > 1) {
		### Add to first row
		# strValName2 <- strValName
		# isAlreadyInTable = strValName%in%names(objREPORT@tblReport)
		# iCount = 1
		# while(isAlreadyInTable) {
			# strValName = paste(strValName2,".",iCount,sep="")
			# isAlreadyInTable = strValName%in%names(objREPORT@tblReport)
			# iCount = iCount + 1
		# }

		objREPORT@tblReport 	<- cbind(objREPORT@tblReport, NewCol=strVal, stringsAsFactors = FALSE)
		names(objREPORT@tblReport)[ncol(objREPORT@tblReport)] <- strValName
		
	} else {
		### Add to currentRow
		# iMatch = match(strValName, names(objREPORT@tblReport))
		# objREPORT@tblReport[objREPORT@iCurrentRow,iMatch] <- strVal
		
		iRow <- objREPORT@iCurrentRow
		iCol <- objREPORT@iCurrentCol
		while(ncol(objREPORT@tblReport)<iCol) {
			objREPORT@tblReport <- cbind(objREPORT@tblReport,NewCol=rep(NA,nrow(objREPORT@tblReport)))
			names(objREPORT@tblReport)[ncol(objREPORT@tblReport)] <- strValName
		}
		objREPORT@tblReport[iRow,iCol] <- strVal
	}
	
	#write.table(objREPORT@tblReport,objREPORT@fileReport,row.names=F,quote=F,sep="\t")
	
	objREPORT@iCurrentCol <- objREPORT@iCurrentCol + 1
	
	isDuplCol = duplicated(names(objREPORT@tblReport))
	if(any(isDuplCol)) {
		strDuplCol = names(objREPORT@tblReport)[which(isDuplCol)]
		stop(paste(" EASY ERROR:REPORT\n Created REPORT column \n",paste(strDuplCol,collapse=","), "\n is duplicated in the REPORT. \n PLease use a different name or (for MERGE,ADJUSTALLELES,AFCHECK,GC or RENAMEMARKER) influence the naming by setting --strTag!", sep=""))
	}
	
	return(objREPORT)
	

})
## from given single numeric
setMethod("REPORT.addval", signature(obj1 = "REPORT", obj2 = "character",obj3="numeric"), function(obj1,obj2,obj3, ...) {

	objREPORT <- REPORT.addval(obj1,obj2,as.character(obj3))
	return(objREPORT)

})


##### Method to SET specific values to Report 
setGeneric("REPORT.setval", function(obj1,obj2,obj3, ...) standardGeneric("REPORT.setval"))
###
setMethod("REPORT.setval", signature(obj1 = "REPORT", obj2 = "character",obj3="character"), function(obj1,obj2,obj3, ...) {
	
	objREPORT 	<- obj1
	strValName 	<- obj2
	strVal	 	<- obj3
	
	iRow <- objREPORT@iCurrentRow
	iCol <- match(strValName, names(objREPORT@tblReport))
	
	if(!is.na(iCol) & length(iCol)==1) objREPORT@tblReport[iRow,iCol] <- strVal
	
	return(objREPORT)

})
## from given single numeric
setMethod("REPORT.setval", signature(obj1 = "REPORT", obj2 = "character",obj3="numeric"), function(obj1,obj2,obj3, ...) {

	objREPORT <- REPORT.setval(obj1,obj2,as.character(obj3))
	return(objREPORT)

})


REPORT.getval <- function(objREPORT, strColName) {		
	
	iRow <- objREPORT@iCurrentRow
	iCol <- match(strColName, names(objREPORT@tblReport))
	
	strValOut = objREPORT@tblReport[iRow, iCol]
	
	return(strValOut)
}

REPORT.resort <- function(objREPORT, strColName, iNewPos) {		
	
	tblIn <- tblOut <- objREPORT@tblReport
	iOldPos = match(strColName, names(tblIn))
	if(!is.na(iOldPos) & length(iOldPos) == 1) {
		aiSortOld = 1:ncol(tblIn)
		aiSortNew = aiSortOld[-iOldPos]
		if(iNewPos > 1 & iNewPos < length(aiSortNew)) {
			aiSortOut = c(aiSortNew[1:(iNewPos-1)], iOldPos, aiSortNew[(iNewPos):length(aiSortNew)])
		} else if(iNewPos == 1) {
			aiSortOut = c(iOldPos, aiSortNew)
		} else if(iNewPos == length(aiSortOld)) {
			aiSortOut = c(aiSortNew, iOldPos)
		} else 	aiSortOut = 1:ncol(tblIn)
		
		tblOut <- tblIn[,aiSortOut]
		objREPORT@tblReport <- tblOut
	}
	return(objREPORT)
}

REPORT.write <- function(objREPORT) {		
	
	write.table(objREPORT@tblReport,objREPORT@fileReport,row.names=F,quote=F,sep="\t")

}

REPORT <- function(fileOutBody){ 
	## Wrapper for class definition
	#REPORTout <- setREPORT(new("REPORT", fileReport = paste(objECF@pathOut,objECF@aSetTags[iSet],".EASY.REPORT.txt",sep="")))
	REPORTout <- setREPORT(new("REPORT", fileReportBody = fileOutBody))
	return(REPORTout)
	
}

# setValidity("ADDCOL", function(object){
	# print("ADDCOL-CHECK")
	
	
	
	# print(TRUE)
	# return(TRUE)
# })


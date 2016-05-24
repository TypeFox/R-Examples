setClass("EASYMERGE",
	representation = representation(
						### EASYIN PARAMS
						fileIn						=	"character",
						strMissing					= 	"character",
						strSeparator				= 	"character",
						acolIn						=	"character",
						acolInClasses				=	"character",
						acolNewName					=	"character",
						### Set within class
						aHeader						= 	"character",
						aClasses					= 	"character",
						aHeaderRead					= 	"character",
						aClassesRead				= 	"character",
						tblGWA						=	"data.frame",
						colMerge					=	"character",
						strSuffix					=	"character",
						iMergeID					=	"numeric"
						),
	prototype = prototype(
						### EASYIN PARAMS
						fileIn						=	"",
						strMissing					= 	"NA",
						strSeparator				= 	"WHITESPACE",
						acolIn						=	"",
						acolInClasses				=	"",
						acolNewName				=	"",
						#### Set withinh class
						aHeader						= 	"",
						aClasses					= 	"",
						aHeaderRead					= 	"",
						aClassesRead				= 	"",
						tblGWA						=	data.frame(),
						colMerge					=	"",
						strSuffix					=	"",
						iMergeID					=	0
						)
)
EASYMERGE.easymerge <- function(objEM, strConfigCommand, icount.GWADATA) {

	#aEqcSlotNamesIn = c("fileIn","fileInShortName", "fileInTag", "fileInType", "strMissing", "strSeparator", "acolIn", "acolInClasses", "pathOut")
	aEqcSlotNamesIn = c("fileIn", "strMissing", "strSeparator", "acolIn", "acolInClasses", "acolNewName", "colMerge", "strSuffix")
	#aEcfSlotNamesIn = c("arcdAddCol", "astrAddColNames")

	objEqcReader <- EqcReader(strConfigCommand,aEqcSlotNamesIn)
	
	for(i in 1:length(objEqcReader@lsEqcSlotsOut)) {
		tmpSlot <- names(objEqcReader@lsEqcSlotsOut)[i]
		tmpSlotVal <- objEqcReader@lsEqcSlotsOut[[i]]
		
		if(all(!is.na(tmpSlotVal))) slot(objEM, tmpSlot) <- tmpSlotVal
	}
	
	objEM@iMergeID <- icount.GWADATA
	
	return(objEM)
}

EASYMERGE.init <- function(object) {
		

	####################################################################
	##### Reset file separator
	if(object@strSeparator == "TAB") 		object@strSeparator <- "\t"
	if(object@strSeparator == "WHITESPACE") object@strSeparator <- ""
	if(object@strSeparator == "SPACE") 		object@strSeparator <- " "
	if(object@strSeparator == "COMMA") 		object@strSeparator <- ","
	
	##### Check File Separator
	if(all(object@strSeparator != c("\t", "", " ", ",")))
		stop("EASY ERROR:\n Wrong File separator defined.\n Please use TAB, WHITESPACE, SPACE or COMMA!")

	####################################################################
	##### Check availablity of GWA file
	if(!file.exists(object@fileIn)) {
		stop(paste("EASY ERROR:\n File \n",object@fileIn,"\n does not exist!!!\n", sep=""))
	}
	
	####################################################################
	##### Set object@aHeader

	strHeaderTmp = scan(file=object@fileIn,what="character",n=1,sep="\n",quiet=TRUE)

	if(object@strSeparator == "\t") {
		isTAB = grepl("\t", strHeaderTmp, fixed=T)
		if(!isTAB) 
			warning(paste("EASY WARNING:\n There is no TAB in the header of file \n",object@fileIn,"\n Please make sure that you have defined the correct delimiter!" ,sep="" ))
		aHeaderTmp <- strsplit(strHeaderTmp,"\t")[[1]]
		#object@ls_afileInHeaders[[i]][k] <- aHeaderTmp
	} else if(object@strSeparator == ",") {
		isCOMMA = grepl(",", strHeaderTmp, fixed=T)
		if(!isCOMMA) 
			warning(paste("EASY WARNING:\n There is no COMMA in the header of file \n",object@fileIn,"\n Please make sure that you have defined the correct delimiter!" ,sep="" ))
		aHeaderTmp <- strsplit(strHeaderTmp,",")[[1]]
		#object@ls_afileInHeaders[[i]][k] <- aHeaderTmp
	} else if(object@strSeparator == " ") {
		isSPACE = grepl(" ", strHeaderTmp, fixed=T)
		if(!isSPACE) 
			warning(paste("EASY WARNING:\n There is no SPACE in the header of file \n",object@fileIn,"\n Please make sure that you have defined the correct delimiter!" ,sep="" ))
		aHeaderTmp <- strsplit(strHeaderTmp," ")[[1]]
		#object@ls_afileInHeaders[[i]][k] <- aHeaderTmp
	} else {
		#### Whitespace consisting of spaces and tabs
		strHeaderTmp2 <- gsub("\t", " ", strHeaderTmp)
		aHeaderTmp <- strsplit(strHeaderTmp2," ")[[1]][strsplit(strHeaderTmp2," ")[[1]]!=""]
		
	}
	
	### Allow for P-value -> P.value
	
	aHeaderTmpOld = aHeaderTmp
	aHeaderTmp <- gsub("-",".",aHeaderTmp)
	if(any(grepl("-",aHeaderTmpOld))) 
		warning(paste("EASY WARNING:\n Columns \n",paste(aHeaderTmpOld[which(grepl("-",aHeaderTmpOld))],collapse=","),"\n contain a '-' that will be renamed to '.'! Therefore the new column names \n",paste(aHeaderTmp[which(grepl("-",aHeaderTmpOld))],collapse=","),"\n must be used throughout the ecf file!" ,sep="" ))
	
	object@aHeaderRead <- object@aHeader <- aHeaderTmp
	
	
	####################################################################	
	#### Set object@aClasses to enable fast reading of input file
	## 
	if(length(object@acolIn) != length(object@acolInClasses) & object@acolInClasses[1] != "")
		stop(paste("EASY ERROR:EASYMERGE\n Length of --acolIn differs from length of --acolInClasses for file\n",object@fileIn,"\n Please check DEFINE or EASYIN statements !!!", sep=""))
	
	if(!all(object@acolNewName == "")) {
		## acolNewName defined
		if(length(object@acolIn) != length(object@acolNewName))
			stop(paste("EASY ERROR:EASYMERGE\n Length of --acolIn differs from length of --acolNewName for file\n",object@fileIn,"\n Please check DEFINE or EASYIN statements !!!", sep=""))	
	}	
	
	
	aClassesTmp <- rep("NULL",length(aHeaderTmp))
	
	#if(all(object@acolIn == "")) {
	if(object@acolInClasses[1] == "") {
		## acolInClasses not defined
		## use best guess class from first 10 rows for all columns
		tbl_10rows <- read.table(object@fileIn, nrows = 10, header=T, sep = object@strSeparator, na.strings = object@strMissing, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "")
		aClasses_10rows <- sapply(tbl_10rows, class)
		#aClassesTmp[is.na(aiMatchHeader)] <- aClasses_10rows[is.na(aiMatchHeader)] ##best guess
		aClassesTmp <- aClasses_10rows
	} else {
		## acolIn defined for a subset of columns
		## only use defined columns
		
		#aiMatchColIn = match(object@acolIn, object@aHeader)
		aiMatchColIn = match(tolower(object@acolIn), tolower(object@aHeader))
		if(any(is.na(aiMatchColIn)))
			stop(paste("EASY ERROR:EASYMERGE\n Defined column \n",paste(object@acolIn[which(is.na(aiMatchColIn))],collapse=";")," not available in file\n",object@fileIn,"\n Please check !!!", sep=""))
		aClassesTmp[aiMatchColIn] <- object@acolInClasses
	}
	
	object@aClassesRead <- object@aClasses <- aClassesTmp
	
	#### Check class definitions
	
	isClassOk = object@aClasses%in%c("character","numeric","integer","double","logical","NULL")
	
	if(any(!isClassOk)) 
		stop(paste("EASY ERROR:EASYMERGE\n Class \n",paste(object@aClasses[which(!isClassOk)],collapse="\n")," not defined\n Please define class 'character','numeric','double','logical', 'integer' or 'NULL' for colums\n ",paste(object@aHeader[which(!isClassOk)],collapse="\n")," !!!", sep=""))

		
	return(object)
}

EASYMERGE.read <- function(object) {

	cat(paste(" + Reading ",object@fileIn, "... \n"))

	object@tblGWA 	<- tryCatch(
		read.table(object@fileIn, header=T, sep = object@strSeparator, na.strings = object@strMissing, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", colClasses = object@aClassesRead),
		error = function(err) {
			strError = err$message
			val=strsplit(strError,"'",fixed=T)[[1]][length(strsplit(strError,"'",fixed=T)[[1]])]
			g=scan(file = object@fileIn, what=character(0), n = -1,sep = "\n",quiet=TRUE)
			iRow = which(grepl(paste(val,object@strSeparator,sep=""),g,fixed=T) | grepl(paste(val,"\n",sep=""),g,fixed=T))[1]
			stop(paste(strError,"\n EASY ERROR:\n Cannot read '",val,"' from row '",iRow,"' !!!\n Please specifiy correct column class in --acolInClasses .\n ", sep=""))
		}
	)
	
	if(dim(object@tblGWA)[1]==0)
		stop(paste("EASY ERROR:EASYMERGE\n There are no rows available in \n",object@fileIn,"\n The file is empty!!!\n", sep=""))
	
	iRemoveHead = which(object@aClassesRead == "NULL")
	if(length(iRemoveHead)>0) {
		object@aHeader <- object@aHeaderRead[-iRemoveHead]
		object@aClasses <- object@aClassesRead[-iRemoveHead]
	}
	
	if(all(object@acolIn != "")) {
		#  Sort according to acolIn, case unsensitive!
		
		iMatchSort=match(tolower(object@acolIn),tolower(object@aHeader))
		# Resort:
		object@aHeader = object@aHeader[iMatchSort]
		object@aClasses = object@aClasses[iMatchSort]
		object@tblGWA = object@tblGWA[,iMatchSort]
		# Rename:
		object@aHeader <- object@acolIn
		names(object@tblGWA) <- object@acolIn
		
		if(length(object@acolNewName)==length(object@acolIn)) {
			## acolNewName defined
			object@aHeader <- object@acolNewName
			names(object@tblGWA) <- object@acolNewName
		}
	
	}
	
	return(object)
}


EASYMERGE.read.10rows <- function(object) {

	#object@tblGWA 	<- read.table(object@fileIn, nrows = 10, header=T, sep = object@strSeparator, na.strings = object@strMissing, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", colClasses = object@aClasses)	
	
	object@tblGWA 	<- tryCatch(
		read.table(object@fileIn, nrows = 10, header=T, sep = object@strSeparator, na.strings = object@strMissing, stringsAsFactors=FALSE, strip.white = TRUE, comment.char = "", colClasses = object@aClasses),
		error = function(err) {
			strError = err$message
			val=strsplit(strError,"'",fixed=T)[[1]][length(strsplit(strError,"'",fixed=T)[[1]])]
			g=scan(file = object@fileIn, what=character(0), n = -1,sep = "\n",quiet=TRUE)
			iRow = which(grepl(paste(val,object@strSeparator,sep=""),g,fixed=T) | grepl(paste(val,"\n",sep=""),g,fixed=T))[1]			
			stop(paste(strError,"\n EASY ERROR:\n Cannot read '",val,"' from row '",iRow,"' !!!\n Please specifiy correct column class in --acolInClasses .\n ", sep=""))
		}
	)
	
	if(dim(object@tblGWA)[1]==0)
		stop(paste("EASY ERROR:EASYMERGE\n There are no rows available in \n ",object@fileIn,"\n The file is empty!!!\n", sep=""))
	
	iRemoveHead = which(object@aClassesRead == "NULL")
	if(length(iRemoveHead)>0) {
		object@aHeader <- object@aHeaderRead[-iRemoveHead]
		object@aClasses <- object@aClassesRead[-iRemoveHead]
	}
	
	if(all(object@acolIn != "")) {
		#  Sort according to acolIn, case unsensitive!
		
		iMatchSort=match(tolower(object@acolIn),tolower(object@aHeader))
		# Resort:
		object@aHeader = object@aHeader[iMatchSort]
		object@aClasses = object@aClasses[iMatchSort]
		object@tblGWA = object@tblGWA[,iMatchSort]
		# Rename:
		object@aHeader <- object@acolIn
		names(object@tblGWA) <- object@acolIn
		
		if(length(object@acolNewName)==length(object@acolIn)) {
			## acolNewName defined
			object@aHeader <- object@acolNewName
			names(object@tblGWA) <- object@acolNewName
		}
	
	}
	
	return(object)
}

EASYMERGE.GWADATA.valid <- function(objEASYMERGE, objGWA) {
	
	isAv <- objEASYMERGE@colMerge %in% objGWA@aHeader
	if(!isAv)
		stop(paste(" EASY ERROR:EASYMERGE\n Defined column colMerge \n",objEASYMERGE@colMerge, "\n is not available in GWA data-set \n",objGWA@fileIn,"\n PLease specify correct column name.", sep=""))

	### Check duplicate colnames in merged data set

						
	aHeadIn <- objGWA@aHeader
	iMarkerIn = match(objEASYMERGE@colMerge, aHeadIn)
	aHeadIn[-iMarkerIn] <- paste(aHeadIn[-iMarkerIn],"",sep="")
	
	aHeadRef <- objEASYMERGE@aHeader
	iMarkerRef = match(objEASYMERGE@colMerge, aHeadRef)
	aHeadRef[-iMarkerRef] <- paste(aHeadRef[-iMarkerRef],objEASYMERGE@strSuffix,sep="")
	
	aHeadIntersect = intersect(aHeadIn[-iMarkerIn],aHeadRef[-iMarkerRef])
	if(length(aHeadIntersect)>0)
		stop(paste(" EASY ERROR:EASYMERGE\n Column \n",paste(aHeadIntersect,collapse="\n"), 
					"\n will be duplicated after merging data to \n",objGWA@fileIn,
					"\n PLease specify correct suffixes strSuffix or rename columns in one file.", sep="")
			)
}

EASYMERGE.run <- function(objEASYMERGE, objGWA) {

	 objGWA.merged <- GWADATA.merge(objGWA, objEASYMERGE, 
						strSuffix.In = "", 
						strSuffix.Add = objEASYMERGE@strSuffix, 
						blnAll.In = TRUE, 
						blnAll.Add = FALSE, 
						strBy.In = objEASYMERGE@colMerge, 
						strBy.Add = objEASYMERGE@colMerge
					) 
		
	return(objGWA.merged)
}

################################################################################################################################
################################################################################################################################
###### Wrapper for class setting
################################################################################################################################
##### Wrapper for constructing the object WITH validity checks
#EASYMERGE <- function(fileIn, fileInTag, strMissing, strSeparator, colMarker, acolPrimaryKey, aHeader, aClasses){ 
EASYMERGE <- function(){ 
	## Wrapper/constructor for class definition
	EASYMERGEout <- new("EASYMERGE")
	return(EASYMERGEout)

}

################################################################################################################################
################################################################################################################################

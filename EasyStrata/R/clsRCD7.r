setClass("RCD",
	representation = representation(
						rcdIn			=	"character",
						valOut			=	"NULL"
						),
	prototype = prototype(
						rcdIn			=	"",
						valOut			=	NULL
						)
	#contains = c("ADDCOL", "EVALSTAT") # inherits from ADDCOL class
	#contains = c("ECF5","ADDCOL") # inherits from ADDCOL class
	### R checks for validity of ADDCOL, when validObject(ECF) is called!
)
################################################################################################################################
################################################################################################################################
###### Setting of the ECF object
################################################################################################################################
##### Constructor for the object

setGeneric("setRCD", function(object) standardGeneric("setRCD"))
setMethod("setRCD", signature = (object = "RCD"), function(object) {
	
	object@rcdIn <- object@rcdIn
	
	return(object)
})

RCD.eval <- function(objRCD, objGWA){
   
	
	for(iCol in 1:ncol(objGWA@tblGWA)) {
		strColName = names(objGWA@tblGWA)[iCol]
		
		if(grepl(strColName,objRCD@rcdIn)) {
			#assign(strColName,as.numeric(as.character(tblIn[,iColIn])))
			
			assign(strColName,objGWA@tblGWA[,iCol])
			expressionCast <- parse(text=paste("class(",strColName,")<-'",objGWA@aClasses[iCol],"'",sep=""))
			eval(expressionCast)
		}
	}
	expressionIn=parse(text=objRCD@rcdIn)	
	
	valOut <- tryCatch(	
			eval(expressionIn),
			error = function(err) stop(paste("EASY ERROR:\n Failed to evaluate rcd \n",objRCD@rcdIn,"\nwith file\n",objGWA@fileIn,"\nPlease check Column Names used!!!\n",err, sep=""))
	)
	
    return(valOut)
}


RCD.eval.report <- function(objRCD, objREPORT){
   
	
	for(iCol in 1:ncol(objREPORT@tblReport)) {
		strColName = names(objREPORT@tblReport)[iCol]
		
		if(grepl(strColName,objRCD@rcdIn)) {
			#assign(strColName,as.numeric(as.character(tblIn[,iColIn])))
			
			assign(strColName,objREPORT@tblReport[,iCol])
			
			suppressWarnings({
			expressionCast <- parse(text=paste("class(",strColName,")<-'numeric'",sep=""))
			eval(expressionCast)
			})

			if(eval(parse(text=paste("all(is.na(",strColName,"))",sep="")))) {
				
				assign(strColName,objREPORT@tblReport[,iCol])
				expressionCast <- parse(text=paste("class(",strColName,")<-'character'",sep=""))
				# expressionCast <- parse(text=paste("class(",strColName,")<-'numeric'",sep=""))
				eval(expressionCast)
				
				## Recheck, still all NA?:
				if(eval(parse(text=paste("all(is.na(",strColName,"))",sep="")))) {
				
					assign(strColName,objREPORT@tblReport[,iCol])
					suppressWarnings({
					expressionCast <- parse(text=paste("class(",strColName,")<-'numeric'",sep=""))
					eval(expressionCast)
					})
					### this suppresses error when plotting solely NA
				}
			}
			
		}
	}
	expressionIn=parse(text=objRCD@rcdIn)	
	
	valOut <- tryCatch(	
			eval(expressionIn),
			error = function(err) stop(paste("EASY ERROR:\n Failed to evaluate rcd \n",objRCD@rcdIn,"\nwith report\n",objREPORT@fileReport,"\nPlease check Column Names used!!!\n",err, sep=""))
	)
	
    return(valOut)
}

#validRCD <- function(object, aHeader, aClasses, strFileName){
# validRCD <- function(objRCD, objGWA){
   
	# for(iCol in 1:ncol(objGWA@tblGWA)) {
		# strColName = names(objGWA@tblGWA)[iCol]
		# if(grepl(strColName,objRCD@rcdIn)) {
			# #assign(strColName,as.numeric(as.character(tblIn[,iColIn])))
			
			# assign(strColName,NULL)
			# expressionCast <- parse(text=paste("class(",strColName,")<-'",objGWA@aClasses[iCol],"'",sep=""))
			# eval(expressionCast)
		# }
	# }
	# expressionIn=parse(text=objRCD@rcdIn)	
	
	# tryCatch(	
			# suppressWarnings(eval(expressionIn)),
			# error = function(err) stop(paste("EASY ERROR:\n Failed to evaluate rcd \n ",objRCD@rcdIn,"\n with file\n ",objGWA@fileIn,"\n Please check Column Names and Classes used!!!\n",err, sep=""))
	# )
	
    # return(TRUE)
	
# }

################################################################################################################################
################################################################################################################################
###### Wrapper for class setting
################################################################################################################################
##### Wrapper for constructing the object WITH validity checks
RCD <- function(rcdIn){
	## Wrapper for class definition
	RCDout <- setRCD(new("RCD", rcdIn = rcdIn))
	#RCDout.valid <- validRCD(RCDout)
	return(RCDout)
	#return(RCDout)
	
}
################################################################################################################################
################################################################################################################################

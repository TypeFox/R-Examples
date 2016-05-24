setClass("EqcReader",
	representation = representation(
						strEqcCommand			=	"character",
						aEqcSlotNamesIn			=	"character",
						lsEqcSlotsOut			=	"list"
						),
	prototype = prototype(
						strEqcCommand			=	"",
						aEqcSlotNamesIn			=	"",
						lsEqcSlotsOut			=	list()
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

setGeneric("setEqcReader", function(object) standardGeneric("setEqcReader"))
setMethod("setEqcReader", signature = (object = "EqcReader"), function(object) {
	
	# if(!file.exists(object@fileECF)) 
		# stop(paste("EASY ERROR:\necf-File",object@fileECF,"\ndoes not exist!", sep=""))
		
	# ##### Read EQC data from config file
	# print("Reading EasyQC config-data...")
	# aEcfRowIn	<-	scan(file = object@fileECF, what=character(0), n = -1, comment.char = "#",  strip.white = TRUE,sep = "\n")
	
	#object@aEcfRowIn
	#object@aEqcSlotNamesIn
	
	strEqcCommand		<-	object@strEqcCommand
	aEqcSlotNamesIn 	<- 	object@aEqcSlotNamesIn
	lsEqcSlotsOut 		<- list()
	
	astrEqcCommand = strsplit(strEqcCommand,"--")[[1]]

	if(toupper(astrEqcCommand[1]) != astrEqcCommand[1]) 
		stop(paste(" EASY ERROR:\n ",astrEqcCommand[1]," is not a EasyQC function: \n Please correct, uncomment or remove this row."))
		
	
	if(length(astrEqcCommand) >1) {
		for(i in 2:length(astrEqcCommand)) {
			
			strCommandPart = astrEqcCommand[i]
			
			aCommandSplit = strsplit(strCommandPart," ")[[1]]
			Vname	<- aCommandSplit[1]
			###
			# Update, remove leading spaces
			if(any(aCommandSplit == "")) aCommandSplit = aCommandSplit[-which(aCommandSplit == "")]
			###
			Vvalue	<- ifelse(length(aCommandSplit)>1,paste(aCommandSplit[2:length(aCommandSplit)],collapse=" "),"NA")
						
			if(Vname%in%aEqcSlotNamesIn) {
				if(	substring(Vname,1,3) == "str" | 
					substring(Vname,1,3) == "col" | 
					substring(Vname,1,3) == "rcd" | 
					substring(Vname,1,4) == "path" | 
					substring(Vname,1,4) == "file" ) {
					
					### Property of type string
					Vvalue[Vvalue == "NA"] <- NA
					#slot(object, Vname) <- Vvalue
					eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=Vvalue",sep="")))
					
					
				} else if(	substring(Vname,1,4) == "astr" | 
							substring(Vname,1,4) == "acol" | 
							substring(Vname,1,4) == "arcd" |
							substring(Vname,1,5) == "afile" ) {
					
					### Property of type string array
					#astrTemp = strsplit(Vvalue,"--")[[1]]
					astrTemp = strsplit(Vvalue,";")[[1]]
					astrTemp[astrTemp == "NA"] <- NA
					#slot(object, Vname) <- astrTemp
					eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=astrTemp",sep="")))
					
				} else if(	substring(Vname,1,3) == "num" ) {
					
					### Property of type numeric
					#ls_ECF[[iProp]] <- as.numeric(Vvalue)
					# slot(object, Vname) <- tryCatch(
						# if(Vvalue == "NA") NA
						# else as.numeric(Vvalue), 
						# warning = function(e) {
							# stop(paste(e,Vname,"=",Vvalue,"not defined:\n  Please define numeric value for", Vname))
						# }
					# )
					VvalueTmp <- tryCatch(
						if(Vvalue == "NA") NA
						else as.numeric(Vvalue), 
						error = function(err) {
							stop(paste(err,"\n EASY ERROR:\n ",Vname," = ",Vvalue," not defined:\n Please define numeric value for ", Vname,sep=""))
						},
						warning = function(err) {
							stop(paste(err,"\n EASY ERROR:\n ",Vname," = ",Vvalue," not defined:\n Please define numeric value for ", Vname,sep=""))
						}
					)
					eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=VvalueTmp",sep="")))
					
				} else if(	substring(Vname,1,4) == "anum" ) {
					
					### Property of type string array
					#astrTemp = strsplit(Vvalue,"--")[[1]]
					astrTemp = strsplit(Vvalue,";")[[1]]
					astrTemp[astrTemp == "NA"] <- NA
					VvalueTmp <- tryCatch(
						as.numeric(astrTemp),
						error = function(err) {
							stop(paste(err,"\n EASY ERROR:\n ",Vname," = ",Vvalue," not defined:\n Please define numeric value for ", Vname,sep=""))
						},
						warning = function(err) {
							stop(paste(err,"\n EASY ERROR:\n ",Vname," = ",Vvalue," not defined:\n Please define numeric value for ", Vname,sep=""))
						}
					)
					eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=VvalueTmp",sep="")))
					
				# } else if(	substring(Vname,1,5) == "aastr" |
							# substring(Vname,1,5) == "aacol" | 
							# substring(Vname,1,5) == "aarcd" ) {
					
					# lsstrTemp = as.list(strsplit(Vvalue,";")[[1]])
					# lsstrTemp2 = lapply(lsstrTemp,function(x) strsplit(x,"--")[[1]])
					# lsstrTemp3 = lapply(lsstrTemp2,function(x) x[x == "NA"] <- NA)
					# #slot(object, Vname) <- lsstrTemp3
					# eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=lsstrTemp3",sep="")))
				
				} else if( substring(Vname,1,3) == "bln" ) {
					
					### Property of type boolean
					if(!(Vvalue == "NA" | Vvalue == "0" | Vvalue == "1")) {
						stop(paste(" EASY ERROR:\n ",Vname," = ",Vvalue," not defined:\n Please use logical value (0,1) or NA for ", Vname, sep=""))
					}
					
					#slot(object, Vname) <- ifelse(Vvalue == "1", TRUE, ifelse(Vvalue == "0", FALSE, NA))
					VvalueTmp <- ifelse(Vvalue == "1", TRUE, ifelse(Vvalue == "0", FALSE, NA))
					eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=VvalueTmp",sep="")))
					
				} else if(	substring(Vname,1,4) == "abln" ) {
					
					### Property of type boolean array					
					#astrTemp = strsplit(Vvalue,"--")[[1]]
					astrTemp = strsplit(Vvalue,";")[[1]]
					
					isOk = astrTemp == "NA" | astrTemp == "0" | astrTemp == "1"
					if(!all(isOk)) {
						stop(paste(" EASY ERROR:\n ",Vname," = ",paste(astrTemp[!isOk], collapse = ",")," not defined:\n Please use logical value (0,1) or NA for each element of ", Vname, sep=""))
					}

					#slot(object, Vname) <- ifelse(astrTemp == "1", TRUE, ifelse(astrTemp == "0", FALSE, NA))
					VvalueTmp <- ifelse(astrTemp == "1", TRUE, ifelse(astrTemp == "0", FALSE, NA))
					eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=VvalueTmp",sep="")))
				}	
				# } else if(	substring(Vname,1,5) == "aabln" ) {
					
					# lsstrTemp = as.list(strsplit(Vvalue,";")[[1]])
					# lsstrTemp2 <- lapply(lsstrTemp,function(x) strsplit(x,"--")[[1]])
					# lapply(lsstrTemp2,
						# function(x) {
							# isOk = x == "NA" | x == "0" | x == "1"
							# if(!all(isOk)) {
										# stop(paste(Vname,"=",paste(x[!isOk], collapse = ","),"not defined:\n  Please define logical value (0,1) or NA for each element of", Vname))
									# }
							# }
					# )
					# #slot(object, Vname) <- lapply(lsstrTemp2,function(x) ifelse(x == "1", TRUE, ifelse(x == "0", FALSE, NA)))
					# VvalueTmp <- lapply(lsstrTemp2,function(x) ifelse(x == "1", TRUE, ifelse(x == "0", FALSE, NA)))
					# eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=VvalueTmp",sep="")))
				# } else if( toupper(Vname) == Vname ) {
					
					# ### Property of type EASY-FUNCTIONS
					# ###	0 - Do not use
					# ### 1 - Use but supress output
					# ### 2 - Use and output as txt (only for Functions that return objGWA)
					# ### 3 - Use and output as gz (only for Functions that return objGWA)
					
					# if(Vvalue == "NA") Vvalue = "0"
					
					# if(!(Vvalue == "0" | Vvalue == "1" | Vvalue == "2" | Vvalue == "3")) {
						# stop(paste(Vname,"=",Vvalue,"not defined:\n  Please define Easy Function output value (0,1,2,3) for Easy Function ", Vname))
					# }
					
					# #slot(object, Vname) <- ifelse(Vvalue == "1", TRUE, ifelse(Vvalue == "0", FALSE, NA))
					# #VvalueTmp <- ifelse(Vvalue == "1", TRUE, ifelse(Vvalue == "0", FALSE, NA))
					# #VvalueTmp <- ifelse(Vvalue == "1", TRUE, ifelse(Vvalue == "0", FALSE, NA))
					# eval(parse(text=paste("lsEqcSlotsOut$",Vname,"=Vvalue",sep="")))
					
				# }
			} else {
				stop(paste(" EASY ERROR:\n ",Vname," is not a EasyQC parameter for function",astrEqcCommand[1],":\n Please correct, uncomment or remove this parameter row.\n Allowed parameters are [",paste(aEqcSlotNamesIn,collapse=","),"]!!!"))
			}
		}
	}
	object@lsEqcSlotsOut <- lsEqcSlotsOut
	
	return(object)
})

################################################################################################################################
################################################################################################################################
###### Wrapper for class setting
################################################################################################################################
##### Wrapper for constructing the object WITH validity checks
EqcReader <- function(strEqcCommand,aEqcSlotNamesIn){
	## Wrapper for class definition
	EqcReader <- setEqcReader(new("EqcReader", strEqcCommand = strEqcCommand, aEqcSlotNamesIn = aEqcSlotNamesIn))
	return(EqcReader)
	# validECF(ECFout)
	# return(ECFout)
	## Identical:
	# ECFin <- new("ECF4", fileECF = fileECFIn) 
	# ECFout <- setECF4(ECFin)
	# return(ECFout)
}
################################################################################################################################
################################################################################################################################

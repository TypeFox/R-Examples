# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 10 dec 2009
# licence	: GPL

# Author(s) : Juhui Wang, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 240                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2011-11-10 11:16:2#$: date of the last spread

#---------------------------------------------------------------------------

##########################################################################
## mtkParsor class,  definition and methods

DEBUG.mtkParsor=FALSE


#' The  mtkParsor class parses the XML file and extracts the information about the factors
#'  and processes used in a sensitivity analysis session.
#' @slot xmlPath the XML file's path and name
#' @exportClass mtkParsor
#' @title The mtkParsor class

setClass(Class="mtkParsor",
		representation=representation(xmlPath="character"),
		validity = function(object){
			if(!file.exists(object@xmlPath)) stop("The file ",object@xmlPath," does not exist!\n")
		} 
)

###########################################################################
## Methods definition
###########################################################################

###########################################################################

#' The constructor
#' @param xmlPath the XML file's path and name.
#' @return an object of class \code{\linkS4class{mtkParsor}}
#' @export mtkParsor

mtkParsor= function(xmlPath) {
	res <- new("mtkParsor", xmlPath=xmlPath)
	return(res)
}



###########################################################################


#' Sets the xml File and  tests if it's openable.
#' @param this an object of class \code{\linkS4class{mtkParsor}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod setXMLFilePath

setMethod(f="setXMLFilePath", signature=c(this="mtkParsor",xmlPath="character"),
		definition=function(this,xmlPath) {
			
			nameThis <- deparse(substitute(this))
			
			if(file.exists(xmlPath))
				this@xmlPath <- xmlPath
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		}) 

###########################################################################

#' Parses the xml file and creates the factors and processes from the XML file.
#' @param this the underlying object of class \code{\linkS4class{mtkParsor}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
#' @title The run method

setMethod(f="run", signature=c(this="mtkParsor", context="mtkExpWorkflow"),
		definition=function(this, context){
			
			nameContext <- deparse(substitute(context))
			
			
			# Parsor the file this@xmlPath, create and modify the slots of the context
			# according to the data extracted from the xmlFile
			
			
			doc <- xmlTreeParse(this@xmlPath,ignoreBlanks=TRUE, trim=TRUE, useInternalNodes=TRUE) 
			
			facteurs <- getNodeSet(xmlRoot(doc), "//mxd:factor")
			
			factors <- list()
			
			for (facteur in facteurs ) {
				
				
				if(DEBUG.mtkParsor)show(facteur)
				
				
				facteur.id<-xmlGetAttr(facteur,"id",converter=str_trim) ## get la valeur de l'attribut "id", si "id" n'existe pas, retourne NULL
				facteur.name <- xmlGetAttr(facteur,"name",converter=str_trim)
				facteur.unit <- xmlGetAttr(facteur,"unit",converter=str_trim)
				if(is.null(facteur.id)) facteur.id <- 'unknown'
				if(is.null(facteur.name)) facteur.name <- 'unknown'
				if(is.null(facteur.unit)) facteur.unit <- ''
				
				
				facteur.domaine <- facteur[["domain"]] ## le noeud "domaine" du facteur, si le noeud "domaine" n'existe pas retourne NULL
				domaine.nominalValue <- xmlGetAttr(facteur.domaine,"nominalValue",converter=str_trim)
				domaine.distributionName <- xmlGetAttr(facteur.domaine,"distributionName",converter=str_trim)
				domaine.valueType <- xmlGetAttr(facteur.domaine,"valueType",converter=str_trim) # str_trim: a function from the package stringr
				
				if(is.null(domaine.nominalValue)) domaine.nominalValue <- 0
				if(is.null(domaine.distributionName)) domaine.distributionName <- 'unknown'
				if(is.null(domaine.valueType)) domaine.valueType <- 'float'
				
#				dom <- mtkDomain(distributionName=domaine.distributionName,
#						nominalValType=xsOff(domaine.valueType),
#						nominalValString=domaine.nominalValue)
				
				## TODO il faut peut etre gerer le cas ou domaine.valueType ne contient pas de : (pas de namespace)
				domValType <- str_sub(domaine.valueType, str_locate(domaine.valueType,":")[1]+1 )
				domNominalVal <- switch(domValType,
						"character" = as.character(domaine.nominalValue),
						"string" = as.character(value),
						"double" = as.double(domaine.nominalValue),
						"float" = as.double(domaine.nominalValue),
						"logical" = as.logical (domaine.nominalValue),
						"integer" = as.integer (value)
				)
				dom <- mtkDomain( distributionName=domaine.distributionName, domNominalVal)				
				
				
				##if(DEBUG.mtkParsor) cat("domaine:*", domaine.distributionName, domaine.nominalValue,domaine.valueType, "*\n")
				if(DEBUG.mtkParsor) cat("domaine:*", domaine.distributionName, domNominalVal, domValType,"*\n")
				
				levels <- list()
				weights <- list()
				
				domaine.levels <- facteur.domaine["level"]
				if(length(domaine.levels)!=0)
					for (level in domaine.levels) {
						if(DEBUG.mtkParsor)show(level) # TODO show plante
						value <- xmlGetAttr(level,"value",converter=str_trim)
						weight <- as.numeric(xmlGetAttr(level,"weight"),converter=str_trim)
						
						levels <- c(levels,value)
						weights <- c(weights, weight)
						
					}
				
				if(length(levels) != 0) {
					mtkL <- mtkLevels('categorical', levels, weights)
					setLevels(dom,mtkL)
				}
				domaine.distributionParameters <- facteur.domaine["distributionParameter"]
				parameters<-list()
				
				if(length(domaine.distributionParameters)!=0)
					for (parameter in domaine.distributionParameters) {
						if(DEBUG.mtkParsor)show(parameter)
						name <- xmlGetAttr(parameter,"name",converter=str_trim)
						value <- xmlGetAttr(parameter,"value",converter=str_trim)
						valueType <- xmlGetAttr(parameter,"valueType",converter=str_trim)
						if(is.null(name)) name <- 'unknown'
						if(is.null(valueType)) valueType <- 'float'
						
#						parameters <- c(parameters, mtkParameter(valName=name,type=xsOff(valueType),valString=value))
						## HR refactoring syntaxe mtkValue, mtkParameter et mtkDomain
						## TODO il faut peut etre gerer le cas ou domaine.valueType ne contient pas de : (pas de namespace)
					##	browser()

						goodValueType <- str_sub(valueType, str_locate(valueType,":")[1]+1 )
						goodValue <- switch(goodValueType,
								"character" = as.character(value),
								"string" = as.character(value),
								"double" = as.double(value),
								"float" = as.double(value),
								"logical" = as.logical (value),
								"integer" = as.integer (value)
						)
						
                                                ## HM, le 12/10/2012
						## assign(name,goodValue, pos=1)
                                                zzz <- mtkParameter()
                                                zzz@name <- name
                                                zzz@val <- goodValue
                                                zzz@type <- typeof(goodValue)
						
						if(DEBUG.mtkParsor) cat("parameter:*", name,  value, goodValueType, "*\n")
						## parameters <- c(parameters, mtkParameter(name))
						parameters <- c(parameters, zzz)
						
					}
				
				if(length(parameters) != 0)setDistributionParameters(dom,parameters)
				
				
				facteur.features <- facteur["feature"]
				features <- list()
				
				if(length(facteur.features)!=0)
					for (feature in facteur.features) {
						if(DEBUG.mtkParsor)show(feature)
						name <- xmlGetAttr(feature,"name",converter=str_trim)
						value <- xmlGetAttr(feature,"value",converter=str_trim)
						valueType <- xmlGetAttr(feature,"valueType", converter=str_trim)
						if(is.null(name)) name <- 'unknown'
						if(is.null(valueType)) valueType <- 'float'
						## HR refactoring syntaxe mtkValue, mtkParameter et mtkDomain
						goodValueType <- str_sub(valueType, str_locate(valueType,":")[1]+1 )
						goodValue <- switch(goodValueType,
								"character" = as.character(value),
								"string" = as.character(value),
								"double" = as.double(value),
								"float" = as.double(value),
								"logical" = as.logical (value),
								"integer" = as.integer (value)
						)
						
                                                ## HM, le 12/10/2012
						## assign(name,goodValue, pos=1)
                                                zzz <- mtkFeature()
                                                zzz@name <- name
                                                zzz@val <- goodValue
                                                zzz@type <- typeof(goodValue)
						
						## features <- c(features, mtkFeature(name))
						features <- c(features, zzz)
					}
				fac<-mtkFactor(name=facteur.name, id=facteur.id, unit=facteur.unit,domain=dom)
				if(length(features) != 0) setFeatures(fac,features)
				
				factors <- c(factors,fac)
				
				
			}
			
			context@expFactors <- mtkExpFactors(expFactorsList=factors)
			
			
			
			if(DEBUG.mtkParsor)show(context@expFactors)
			
			
			
			processes <- getNodeSet(xmlRoot(doc), "//mxd:process")
			
			
			
			for (process in processes ) {
				if(DEBUG.mtkParsor)show(process)
				stage <- xmlGetAttr(process,"stage")
				call  <- xmlValue(process[["call"]])
				
				protocolSeparated <- strsplit(call, "://")
				protocol <- protocolSeparated[[1]][1]
				siteSeparated <- strsplit(protocolSeparated[[1]][2],"/")
				site <- siteSeparated[[1]][1]
				service <- siteSeparated[[1]][2]
				
				
				parameters=process[["parameters"]]["parameter"]
				
				p <- vector(mode="raw", length=0)
				if(length(parameters)!=0)
					for (parameter in parameters) {
						if(DEBUG.mtkParsor)show(parameter)
						name <- xmlGetAttr(parameter,"name",converter=str_trim)
						value <- xmlGetAttr(parameter,"value",converter=str_trim)
						valueType <- xmlGetAttr(parameter,"valueType",converter=str_trim)
						if(is.null(name)) name <- 'unknown'
						if(is.null(valueType)) valueType <- 'float'
						#browser()
						## HR refactoring syntaxe mtkValue, mtkParameter et mtkDomain
						goodValueType <- str_sub(valueType, str_locate(valueType,":")[1]+1 )
						goodValue <- switch(goodValueType,
								"character" = as.character(value),
								"string" = as.character(value),
								"double" = as.double(value),
								"float" = as.double(value),
								"logical" = as.logical (value),
								"integer" = as.integer (value)
						
						)
						
                                                ## HM, le 12/10/2012
						## assign(name,goodValue, pos=1)
                                                zzz <- mtkParameter()
                                                zzz@name <- name
                                                zzz@val <- goodValue
                                                zzz@type <- typeof(goodValue)

						## p <- c(p,mtkParameter(name))
                                                p <- c(p,zzz)
						
					}
				if(stage == "design") cService <- "mtkDesigner"
				if(stage == "evaluate") cService <- "mtkEvaluator"
				if(stage == "analyze") cService <- "mtkAnalyser"
				obj<-new(cService,name=stage, protocol=protocol, site=site, service=service,
						parameters=p,ready=TRUE, state=FALSE, result=NULL)
				setProcess(context,obj,stage)
				
			}
			
			assign(nameContext,context, envir=parent.frame())
			
			return(invisible())
		})


#' Extracts the sub-string B from a string of pattern A:B such xs:integer.
#' @param str a string of pattern A:B such xs:integer
#' @return the sub-string B of str
#' @title The xsOff function

xsOff<-function(str){
	tmp <- strsplit(str, ":")
	return(tmp[[1]][2])
}



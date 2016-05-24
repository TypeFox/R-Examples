# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 8 feb 2011
# licence	: GPL

# Author(s) : Juhui WANG, based a version from Hervé Richard INRA MIA BioSP
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 200                 $: revision number of the last spread
#' $Author:: hrichard         $: author of the last spread
#' $Date:: 2011-10-19 17:20:4#$: date of the last spread

#-----------------------------------------------------------------------
#' Nota :
# HR fusion de rev 25 avec rev 29
#=======================================================================
#' The  mtkDomain class
#' @slot distributionName a string representing the  distribution law
#' @slot nominalValue mtkValue corresponding to nominal value
#' @slot levelList list of mtkLevel objects
#' @slot distributionParameters list of mtkParameter
#' @title The mtkDomain class
#' @exportClass mtkDomain
#====================================================================
# TODO :
# 1 - on devrait manipuler des vecteurs plutot que des listes
# pour s'assurer que l'objet est homogene
# 2 - prevoir surcharge de l'operateur $ pour extraire levelList et
# distributionParameters de l'objet mtkDomain
#====================================================================


setClass("mtkDomain",
		representation=representation
				(
				distributionName="character",
				nominalValue="ANY",
				levels="mtkLevels",
				distributionParameters="list"
		)
)
#===================================================================
# initialize definition

setMethod(f="initialize",
		signature=c("mtkDomain"),
		definition=function(.Object, distributionName="unknown", domainNominalValue=0, distributionParameters=list())
		{
			
			.Object@distributionName <- distributionName
			
			if (is.na(match("mtkValue",(class(domainNominalValue))))){
				nomVal <- mtkValue(name='nominalValue', type = typeof(domainNominalValue), val=domainNominalValue)
			}
			else nomVal<- domainNominalValue
			
			.Object@nominalValue <- nomVal
      		if(distributionName == 'discrete'){
				if(length(distributionParameters)==0) levels <- mtkLevels()
				if(length(distributionParameters)> 0){
					if(exists("mtkLevels", where=distributionParameters)) levels <- distributionParameters['mtkLevels']
					t <- 'categorical'
					if(exists("type", where=distributionParameters)) 
						t <- unlist(distributionParameters['type'])
					l <- list()
					if(exists("levels", where=distributionParameters)) 
						l <- unlist(distributionParameters['levels'])
					w <- numeric(0)
					if(exists("weights", where=distributionParameters))
						 w <- unlist(distributionParameters['weights'])
		
					levels <- mtkLevels(t, l, w)
				}
				.Object@levels <- levels
			}else {
			if(length(distributionParameters) > 0 && class(distributionParameters[[1]])!="mtkParameter")
					distributionParameters <- make.mtkParameterList(distributionParameters)
			
			.Object@distributionParameters=distributionParameters
			}
			return(.Object)
		}
)
#===================================================================
#' The constructor method
#' @param distributionName the distribution name.
#' @param domainNominalValue a mtkValue value or arg to build mtkValue (see mtkValue Constructor)
#' @export mtkDomain

mtkDomain = function(distributionName="unknown", domainNominalValue=0,distributionParameters=list())
{
	
	new(Class="mtkDomain", 
			distributionName=distributionName,
			domainNominalValue=domainNominalValue,
      		distributionParameters=distributionParameters
	)
}
#===================================================================
# Method definition
#' Return the distribution's name .
#' @param this an object of class \code{\linkS4class{mtkDomain}}
#' @title The getDistributionName method
#' @exportMethod getDistributionName
#

setMethod(f="getDistributionName",
		signature=c(this="mtkDomain"),
		definition=function(this="mtkDomain"){

			return (this@distributionName)
		}
)
#
#
#===================================================================
# Method definition
#' The getNominalValueType method
#' @param this the underlying object of class \code{\linkS4class{mtkDomain}}
#' @exportMethod getNominalValueType
#' @title The getNominalValueType method

setMethod(f="getNominalValueType",
		signature=c("mtkDomain"), definition=function(this) {
			return(this@nominalValue@type)
		}
)

#===================================================================
# Method definition
#' The getNominalValueType method
#' @param this the underlying object of class \code{\linkS4class{mtkDomain}}
#' @exportMethod getNominalValueType
#' @title The getNominalValueType method

setMethod(f="getNominalValue",
		signature=c("mtkDomain"), definition=function(this) {
			return(this@nominalValue@val)
		}
)

#===================================================================
# Method definition
#' The getLevels method
#' @param this the underlying object of class \code{\linkS4class{mtkDomain}}
#' @return TRUE or FALSE
#' @exportMethod getLevels
#
# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getLevels", function(this) {standardGeneric("getLevels")})

setMethod(f="getLevels",
		signature=c(this="mtkDomain"),
		definition=function(this="mtkDomain"){
			return (getLevels(this@levels))
		}
)

#===================================================================
# Method definition
#' The getLevels method
#' @param this the underlying object of class \code{\linkS4class{mtkDomain}}
#' @return TRUE or FALSE
#' @exportMethod getLevels
#
# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getLevels", function(this) {standardGeneric("getLevels")})

setMethod(f="getWeights",
		signature=c(this="mtkDomain"),
		definition=function(this="mtkDomain"){
			return (getWeights(this@levels))
		}
)

#===================================================================
# Method definition
#' The getLevels method
#' @param this the underlying object of class \code{\linkS4class{mtkDomain}}
#' @return TRUE or FALSE
#' @exportMethod getLevels
#
# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getLevels", function(this) {standardGeneric("getLevels")})

setMethod(f="getDiscreteDistributionType",
		signature=c(this="mtkDomain"),
		definition=function(this="mtkDomain"){
			return (getType(this@levels))
		}
)
#===================================================================
# Method definition
#' The getDistributionParameters method
#' @param this the underlying object of class \code{\linkS4class{mtkDomain}}
#' @exportMethod getDistributionParameters

#
setMethod(f="getDistributionParameters",
		signature=c(this="mtkDomain"),
		definition=function(this){
			distribParms <- list()
			if(length(this@distributionParameters)>0){
				distribParms <- lapply(this@distributionParameters,getValue)
				n <- lapply(this@distributionParameters,getName)
				names(distribParms) <- n
			}
			return(distribParms)
		}
)
#===================================================================
# Method definition
#' The setLevels method
#' @param this the underlying object of class \code{\linkS4class{mtkDomain}}
#' @param levels a list of mtkLevel objects
#' @exportMethod setLevels

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("setLevels",function(this, levels) {standardGeneric("setLevels")})
setMethod(f="setLevels",
		signature=c(this="mtkDomain",levels="list"),
		definition=function(this, levels){
			nameThis <- deparse(substitute(this))
			if(length(levels)==0) levels <- mtkLevels()
			if(length(levels)> 0){
					t <- 'categorical'
					if(exists("type", where=levels)) 
						t <- unlist(levels['type'])
					l <- list()
					if(exists("levels", where=levels)) 
						l <- unlist(levels['levels'])
					w <- numeric(0)
					if(exists("weights", where=levels))
						w <- unlist(levels['weights'])
					
					levels <- mtkLevels(t, l, w)
				}
				
				
			this@levels <- levels
			
			assign(nameThis, this, envir <- parent.frame())
  			return(invisible())
		}
)

setMethod(f="setLevels",
		signature=c(this="mtkDomain",levels="mtkLevels"),
		definition=function(this, levels){
			nameThis <- deparse(substitute(this))
			this@levels <- levels
			assign(nameThis, this, envir <- parent.frame())
			return(invisible())
		}
)
#===================================================================
# Method definition
#' The setDistributionParameters method
#' @param this the underlying object of class \code{\linkS4class{mtkDomain}}
#' @param list a list of mtkParameter objects
#' @exportMethod setDistributionParameters

setMethod(f="setDistributionParameters",
		signature=c(this="mtkDomain",aDistParamList="list"),
		definition=function(this,aDistParamList ){
			localThis <- deparse(substitute(this))
			if(this@distributionName == 'discrete'){
				if(length(aDistParamList)==0) levels <- mtkLevels()
				if(length(aDistParamList)> 0){
					if(exists("mtkLevels", where=aDistParamList)) levels <- aDistParamList['mtkLevels']
					t <- 'categorical'
					if(exists("type", where=aDistParamList)) 
						t <- unlist(aDistParamList['type'])
					l <- list()
					if(exists("levels", where=aDistParamList)) 
						l <- unlist(aDistParamList['levels'])
					w <- numeric(0)
					if(exists("weights", where=aDistParamList))
						w <- unlist(aDistParamList['weights'])
					
					levels <- mtkLevels(t, l, w)
				}
				this@levels <- levels
			}else {
				if(length(aDistParamList) >0 && class(aDistParamList[[1]]) != "mtkParameter")
					aDistParamList <- make.mtkParameterList(aDistParamList)
			
					this@distributionParameters <- aDistParamList
			}
			assign(localThis, this, envir <- parent.frame())
			return(invisible())
			
		}
	
)


setMethod(f="show", signature=c(object="mtkDomain"),
		definition=function(object){
			print(object)
			invisible()
		}
)


setMethod(f="print", signature=c(x="mtkDomain"), 
		definition=function(x,...){
			## general information 
			cat("------------------------------\n")
			cat('An object of class "mtkDomain" \n')
			cat("------------------------------\n")
			
			## information on the distribution
			cat( 'Distribution : ', getDistributionName(x),  '\n')
			cat( 'Nominal value : ', getNominalValue(x), '\n')
			cat( 'Nominal Value Type : ', getNominalValueType(x), '\n')
			if(getDistributionName(x)=='discrete'){
				cat('Type : ', getDiscreteDistributionType(x), '\n')
				## information on the levels
				cat('Levels : ', unlist(getLevels(x)), '\n')
				cat('Weights : ', unlist(getWeights(x)), '\n')
				## information on the features
			} else {	
				nbp <- length(x@distributionParameters)
				if(nbp > 0){
					cat("Distribution parameters :",nbp," \n")
					print(unlist(getDistributionParameters(x)),...)
				}
			}
			cat("------------------------------\n")
			invisible()
			
		}
)

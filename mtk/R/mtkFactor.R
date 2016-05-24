# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 
# licence	: GPL

# Author(s) : Juhui WANG, 22 Avril 2014
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 264                 $: revision number of the last spread
#' $Author:: hmonod           $: author of the last spread
#' $Date:: 2011-12-15 18:46:2#$: date of the last spread

#-----------------------------------------------------------------------
#' Nota :
#' The  mtkFactor class
#SetLevelList
# Attributs :
#' The mtkFactor class
#' @slot name string, name of the mtkFactor
#' @slot id string, single Id of the mtkFactor
#' @slot unit string, optional
#' @slot type string, a string to specify the data type of the factor's levels
#' @slot domain mtkDomain object for specify the mtkFactor domain
#' @slot featureList list of mtkFeature objects for specific characteristic of a mtkFactor
#' @exportClass mtkFactor
# =======================================================================
setClass(Class="mtkFactor",
		representation=representation(
				name="character",
				id="character",
				unit="character",
				type="character",
				domain="mtkDomain",
				featureList="list"
		),
		prototype=prototype(name="unkown", id="unkown", unit="", type="numeric", domain=mtkDomain(), featureList=list())
)
#===================================================================
# initialize definition
# call with new()
# WARNING : args function need default value!
#
#' The initialize method
#' @param .Object the underlying object
#' @param name a string giving the factor name
#' @param id a string giving the single id of the factor in the code
#' @param unit a string giving the measurement units of the factor levels
#' @paeam type a string giving the data type of the factor's levels
#' @param domain an mtkDomain object
#' @return an object of class mtkFactor

setMethod(
		f="initialize",
		signature=c("mtkFactor"),
		definition=function(.Object, name="unkown", id="unkown", unit="", type="numeric", domain=mtkDomain(), featureList=list())
		{
			.Object@name <- name
			.Object@id <- id
			.Object@type <- type
			.Object@unit <- unit
			.Object@domain <- domain
			if(length(featureList)>0 && class(featureList[[1]]) != 'mtkFeature')
				featureList<- make.mtkFeatureList(featureList)
			.Object@featureList <- featureList
			##validObject(object)#call to validity method
			return(.Object)
		}
)
#===================================================================
#' The constructor
#' @param name a string to name the factor.
#' @param id a string giving the single id of the factor in the code
#' @param unit a string giving the measurement units of the factor levels
#' @param type a string giving the data type associated with the factor
#' @param domain an mtkDomain object
#' @return an object of class \code{\linkS4class{mtkFactor}}
#' @export mtkFactor

mtkFactor = function(name="unkown", id="unkown", unit="", type="numeric", domain=mtkDomain(), featureList=list())
{
	new(Class="mtkFactor",
			name=name,
			id=id,
			unit=unit,
			type = type,
			domain=domain,
			featureList = featureList
			)
}

#===================================================================
# Method definition
#' The setFeatures method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @param aFList a list of mtkFeature objects
#' @exportMethod setFeatures
#
# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("setFeatures",function(this, aFList) {standardGeneric("setFeatures")})
setMethod(f="setFeatures",
		signature=c(this="mtkFactor",aFList="list"),
		definition=function(this, aFList){
     	 localThis <- deparse(substitute(this))
     	 
	  	if(length(aFList)>0 && class(aFList[[1]]) != 'mtkFeature')
			  aFList<- make.mtkFeatureList(aFList)
	  	this@featureList <- aFList
		
      	assign(localThis, this, envir <- parent.frame())
     	return(invisible())

		}
)



################# made by Juhui WANG ########################################################
#' sets the data type for the factor.
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @param type a string indicating the data type associated with the factor.
#' @return invisble()
#' @exportMethod setType
#' @title The setType method

setMethod(f="setType", signature=c(this="mtkFactor", type="character"),
		definition=function(this, type ) {
			nameThis <- deparse(substitute(this))
			this@type <-  type
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})


################## made by Juhui WANG ########################################################
#' gets the type associated with the factor.
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a string indicating the data type associated with the factor.
#' @exportMethod getType
#' @title The getType method

setMethod(f="getType",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
			return (this@type)
		}
)
#

#===================================================================
# Method definition
#' The getName method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return the name of the factor as a string.
#' @exportMethod getName
#
# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getName", function(this) {standardGeneric("getName")})
# Argument
# aMtkFactor : mtkfactor object unknown
# Return a string : the name of the mtkFactor
#
setMethod(f="getName",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getName call from mtkFactor class : argument this not valid, must be a mtkFactor object")
# 			}
			return (this@name)
		}
)

setMethod(f="setName", signature=c(this="mtkFactor", name="character"),
		definition=function(this, name) {
			nameThis <- deparse(substitute(this))
			this@name <-  name
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
# TODO la meme chose avec id : getFactorId
#===================================================================
# Method definition
#' The getDomain method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a mtkDomain object
#' @exportMethod getDomain

#
# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getDomain", function(this) {standardGeneric("getDomain")})
# Argument
# aMtkFactor : mtkfactor object
# Return a mtkDomain object, i.e. the domain of the mtkFactor
#
setMethod(f="getDomain",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getDomain call : argument this not valid, must be a mtkFactor object")
# 			}
			return (this@domain)
		}
)

setMethod(f="setDomain",
		signature=c(this="mtkFactor",domain="mtkDomain"),
		definition=function(this, domain){
			localThis <- deparse(substitute(this))
			
			if(class(domain) != 'mtkDomain')
				stop ("domain must be an object of the class mtkDomain")
			this@domain <- domain
			
			assign(localThis, this, envir <- parent.frame())
			return(invisible())
			
		}
)
#===================================================================
# Method definition
#' The getDistributionName method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a string
#' @exportMethod getDistributionName

#
# move in mtkAllGenerics.R, uncomment otherwise :
# (same signature as mtkFactor)
#setGeneric ("getDistributionName", function(this) {standardGeneric("getDistributionName")})
# Argument
# amtkFactor : mtkFactor object
# Return a string : the name of the mtkDomain
#
setMethod(f="getDistributionName",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getDistributionName call from mtkFactor : argument this not valid, must be a mtkFactor object")
# 			}
			dn <- getDistributionName(this@domain)
			return (dn)
		}
)

#===================================================================
# Method definition
#' The getDistributionName method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a string
#' @exportMethod getDistributionName

#
# move in mtkAllGenerics.R, uncomment otherwise :
# (same signature as mtkFactor)
#setGeneric ("getDistributionName", function(this) {standardGeneric("getDistributionName")})
# Argument
# amtkFactor : mtkFactor object
# Return a string : the name of the mtkDomain
#
setMethod(f="getDistributionNominalValue",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getDistributionName call from mtkFactor : argument this not valid, must be a mtkFactor object")
# 			}
			dn <- getNominalValue(this@domain)
			return (dn)
		}
)

#===================================================================
# Method definition
#' The getDistributionName method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a string
#' @exportMethod getDistributionName

#
# move in mtkAllGenerics.R, uncomment otherwise :
# (same signature as mtkFactor)
#setGeneric ("getDistributionName", function(this) {standardGeneric("getDistributionName")})
# Argument
# amtkFactor : mtkFactor object
# Return a string : the name of the mtkDomain
#
setMethod(f="getDistributionNominalValueType",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getDistributionName call from mtkFactor : argument this not valid, must be a mtkFactor object")
# 			}
			dn <- getNominalValueType(this@domain)
			return (dn)
		}
)

#===================================================================
# Method definition
#' The getDistributionName method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a string
#' @exportMethod getDistributionName

#
# move in mtkAllGenerics.R, uncomment otherwise :
# (same signature as mtkFactor)
#setGeneric ("getDistributionName", function(this) {standardGeneric("getDistributionName")})
# Argument
# amtkFactor : mtkFactor object
# Return a string : the name of the mtkDomain
#
setMethod(f="getDiscreteDistributionType",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getDistributionName call from mtkFactor : argument this not valid, must be a mtkFactor object")
# 			}
			dn <- getDiscreteDistributionType(this@domain)
			return (dn)
		}
)

#===================================================================
# Method definition
#' The getDistributionName method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a string
#' @exportMethod getDistributionName

#
# move in mtkAllGenerics.R, uncomment otherwise :
# (same signature as mtkFactor)
#setGeneric ("getDistributionName", function(this) {standardGeneric("getDistributionName")})
# Argument
# amtkFactor : mtkFactor object
# Return a string : the name of the mtkDomain
#
setMethod(f="getDiscreteDistributionLevels",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getDistributionName call from mtkFactor : argument this not valid, must be a mtkFactor object")
# 			}
			dn <- getLevels(this@domain)
			return (dn)
		}
)

#===================================================================
# Method definition
#' The getDistributionName method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a string
#' @exportMethod getDistributionName

#
# move in mtkAllGenerics.R, uncomment otherwise :
# (same signature as mtkFactor)
#setGeneric ("getDistributionName", function(this) {standardGeneric("getDistributionName")})
# Argument
# amtkFactor : mtkFactor object
# Return a string : the name of the mtkDomain
#
setMethod(f="getDiscreteDistributionWeights",
		signature=c(this="mtkFactor"),
		definition=function(this="mtkFactor"){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getDistributionName call from mtkFactor : argument this not valid, must be a mtkFactor object")
# 			}
			dn <- getWeights(this@domain)
			return (dn)
		}
)

#
#===================================================================
# Method definition
#' The getDistributionParameters method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a named list
#' @exportMethod getDistributionParameters
#
# move in mtkAllGenerics.R, uncomment otherwise :
# setGeneric ("getDistributionParameters", function(this) {standardGeneric("getDistributionParameters")})

setMethod(f="getDistributionParameters",
		signature=c(this="mtkFactor"),
		definition=function(this){
# 			if (!is.mtkFactor(this)) {
# 				stop ("Stop in getDistributionParameters call from mtkFactor: argument This not valid, must be a mtkFactor object")
# 			}
			# HR 01/07/11 correction du return.  HR 11/08/11 mauvaise idée (retour à la syntaxe v148)
			
			return (getDistributionParameters(this@domain))
		}
)

#===================================================================
# Method definition
#' The getDistributionParameters method
#' @param this the underlying object of class \code{\linkS4class{mtkFactor}}
#' @return a named list
#' @exportMethod getDistributionParameters
#
# move in mtkAllGenerics.R, uncomment otherwise :
# setGeneric ("getDistributionParameters", function(this) {standardGeneric("getDistributionParameters")})

setMethod(f="getFeatures",
		signature=c(this="mtkFactor"),
		definition=function(this){
			f <- list()
			if(length(this@featureList)>0){
				f <- lapply(this@featureList,getValue)
				n <- lapply(this@featureList,getName)
				names(f) <- n
			}
			return(f)
		}
)

setMethod(f="getMTKFeatures",
		signature=c(this="mtkFactor"),
		definition=function(this){
			
			return(this@featureList)
		}
)
#    
#)
#===================================================================
### HM, 5/6/11: debut
## HR 09/11 transfo en methode et renommage
## HM 2011-11-03 conservation du renommage, retour a une fonction
#===================================================================
# Method definition
#' The make.mtkFactor function creates a factor from simple R objects 
#' @param name a character string, giving the name of the factor
#' @param id an optional character string, giving the code name of the factor
#' @param unit an optional character string, giving the units of the factor values (to be used in graphics, for example)
#' @param type a string to tell the data type associated with the factor.
#' @param nominal an optional nominal value, character or numeric depending on the distribution
#' @param distribName a character string, giving the R or mtk name of the distribution
#' @param distribPara an optional named list of distribution parameters
#' @param levels an optional character or numeric vector of the factor levels
#' @param weights an optional vector of levels weights
#' @param features optional list of factor features
#' @return an object of class mtkfactor
#' @examples make.mtkFactor("A", distribName="unif", distribPara=list(min=0,max=1))

make.mtkFactor <- function(name="unkown", id="unkown", unit="", type="", nominal=NA, distribName='unknown', distribPara=list(), features=list()) {
  
    mNominal <- mtkValue('nominalvalue', type = type, val = nominal)
	factDomain <-mtkDomain(distributionName=distribName, domainNominalValue=mNominal, distributionParameters=distribPara)
 
    facteur <- mtkFactor(name = name, id = id, unit = unit, type = type, domain = factDomain, featureList = features)
  
 	 return(facteur)
}
#===================================================================
### HM, 5/6/11: fin
#===================================================================

#===================================================================
### HM, 5/11/11
### show method
#===================================================================

#' The show method for an object of class \code{mtkFactor}
#' @param object an object of class \code{mtkFactor}
#' @value nil
setMethod(f="show", signature=c(object="mtkFactor"),
		definition=function(object){
			print(object)
			invisible()
		}
)


#' The show method for an object of class \code{mtkFactor}
#' @param object an object of class \code{mtkFactor}
#' @value nil
setMethod(f="print", signature=c(x="mtkFactor"), 
				definition=function(x,...){
					## general information 
					cat("------------------------------\n")
					cat('An object of class "mtkFactor" \n')
					cat("------------------------------\n")
					cat( 'Name : ', x@name, '\n')
					cat( 'Id : ', x@id, '\n')
					cat( 'Unit : ', x@unit, '\n')
					## information on the distribution
					cat( 'Nominal value : ', getDistributionNominalValue(x), '\n')
					cat('Distribution : ', getDistributionName(x),  '\n')
					if(getDistributionName(x)=='discrete'){
						cat('Type : ', getDiscreteDistributionType(x), '\n')
						## information on the levels
						cat('Levels : ', unlist(getDiscreteDistributionLevels(x)), '\n')
						cat('Weights : ', unlist(getDiscreteDistributionWeights(x)), '\n')
						## information on the features
					} else {	
						nbp <- length(getDistributionParameters(x))
						if(nbp > 0){
							cat("Distribution parameters :",nbp," \n")
							print(unlist(getDistributionParameters(x)),...)
							
						}
					}
					nbf <- length(x@featureList)
					
					if(nbf > 0){
						cat("Features :",nbf," \n")
						print(unlist(getFeatures(x)),...)
					}
					
					cat("------------------------------\n")
					invisible()
          }
          )

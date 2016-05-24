# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 14 jan 2010
# licence	: GPL

# Author(s) : Hervé Richard INRA MIA BioSP
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 228                 $: revision number of the last spread
#' $Author:: hmonod           $: author of the last spread
#' $Date:: 2011-11-07 16:24:2#$: date of the last spread

#-----------------------------------------------------------------------
#
#=======================================================================
# TODO mettre un attribut nom mtkExpFactors ?
# (dans ce cas faire suivre dans le schema)
#===================================================================
# class Declaration
#
mtkExpFactors.DEBUG=FALSE

#' The mtkExpFactors class
#' @slot expFactorsList list of mtkFactor objects
#' @exportClass mtkExpFactors

setClass(Class="mtkExpFactors",
		representation=representation(expFactorsList="list"),
		prototype=prototype(expFactorsList=list())
)

#===================================================================
# initialize definition
# call with new()
# WARNING : args function need default value!
#
#' The initialize method
#' @param .Object the underlying object
#' @param expFactorsList a list of mtkFactor objects
#' @return an object of class mtkExpFactors
# @exportMethod initialize

setMethod(
		f="initialize",
		signature=c("mtkExpFactors"),
		definition=function(.Object, expFactorsList=list())
		{
			.Object@expFactorsList <- expFactorsList
                        
                        ## to ensure that .Object@expFactorsList is correctly named
                        ## VERY IMPORTANT for extractions to work correctly
            fnames <- sapply(.Object@expFactorsList,getName)
            names(.Object@expFactorsList) <- fnames

			##validObject(object)#call to validity method
			return(.Object)
		}
)
#===================================================================
#' The constructor method
#' @param expFactorsList a list of mtkFactor objects
#' @return an object of class mtkExpFactors
#' @export mtkExpFactors
mtkExpFactors <- function(expFactorsList=list())
{
	new(Class="mtkExpFactors", expFactorsList=expFactorsList)
}

# Method definition

#' Assigns a list of mtkFactor objects to the object mtkExpFactors
#' @title The setFactors method
#' @param this an object of class \code{\linkS4class{mtkExpFactors}}
#' @param aFactList a list of mtkFactor objects
#' @exportMethod setFactors

# Return Nil
# Arguments :
# aMtkExpFactors : a valid instance of mtkExpFactors
# aListMtkFactor : a valid list of mtkFactor elements
#
# move in mtkAllGenerics.R, uncomment otherwise
# setGeneric ("setFactors", function(this, aFactList) {standardGeneric("setFactors")})
#
setMethod(f="setFactors",
		signature=c(this="mtkExpFactors",aFactList="list"),
		definition=function(this, aFactList){
			
			if(length(aFactList) >0 && class(aFactList[[1]]) == 'mtkFactor')
			{
				localThis <- deparse(substitute(this))
				this@expFactorsList <- aFactList
				assign(localThis, this, envir <- parent.frame())
				return(invisible())
			}
			else
			{
				stop("Stop in setFactors call : argument aFactList not a list or not a mtkFactor list")
			}
		}
)
# TODO le test sur le premier element de la liste me parait insuffisant
#===================================================================
# Method definition
# TODO addMtkFactor
# Add a factor to the list of the mtkExpFactors
# Return Nil
# Arguments :
# aMtkExpFactors : a valid instance of mtkExpFactors
# aMtkFactor : a valid instance of mtkFactor
#
# move in mtkAllGenerics.R, uncomment otherwise
# setGeneric ("addMtkFactor", function(this, p) {standardGeneric("addMtkFactor")})
#
#===================================================================
# Method definition
# TODO deleteMtkFactor
#
#
# move in mtkAllGenerics.R, uncomment otherwise
# setGeneric ("deleteMtkFactor", function(this, p) {standardGeneric("deleteMtkFactor")})

#===================================================================
# Method definition

#' Gets a list of mtkExpFactors names.
#' @param this the underlying object of class \code{\linkS4class{mtkExpFactors}}
#' @return a list of string incluing all the names of the factors.
#' @title The getNames method
#' @exportMethod getNames
#'
# Argument
# aMtkExpFactors a mtkExpFactors object
# Return
# a list of string including all the names of the mtkExpFactors
#
# move in mtkAllGenerics.R, uncomment otherwise
# setGeneric ("getNames", function(this) {standardGeneric("getNames")})
#
setMethod(f="getNames",
		signature=c(this="mtkExpFactors"),
		definition=function(this="mtkExpFactors"){
			
			mtkFactorsNames <- sapply(this@expFactorsList,getName)
			
			return (mtkFactorsNames)
		}
)


#===================================================================
# Method definition
setMethod(f="getFactors",
		signature=c(this="mtkExpFactors"),
		definition=function(this="mtkExpFactors"){
		return (this@expFactorsList)
		}
)

#===================================================================
# Method definition
setMethod(f="getFactorNames",
		signature=c(this="mtkExpFactors"),
		definition=function(this="mtkExpFactors"){
			
			fNames <- sapply(this@expFactorsList,getName)
			
			return (fNames)
		}
)



#' Gets a list of mtkExpFactors names.
#' @param this the underlying object of class \code{\linkS4class{mtkExpFactors}}
#' @return a list of string incluing all the distributionnames of the domains of mtkExpFactors.
#' @exportMethod getDistributionNames
#' @title The getDistributionNames method
# Argument
# aMtkExpFactors a mtkExpFactors object
# Return
# a list of string incluing all the distributionnames of the domains of mtkFactor objects
#
# move in mtkAllGenerics.R, uncomment otherwise
# setGeneric ("getDistributionNames", function(this) {standardGeneric("getDistributionNames")})
#
setMethod(f="getDistributionNames",
		signature=c(this="mtkExpFactors"),
		definition=function(this="mtkExpFactors"){
			
			distributionNames <- sapply(this@expFactorsList,getDistributionName)

			return (distributionNames)
		}
)

#' Gets a list of mtkExpFactors names.
#' @param this the underlying object of class \code{\linkS4class{mtkExpFactors}}
#' @return a list of string incluing all the distributionnames of the domains of mtkExpFactors.
#' @exportMethod getDistributionNames
#' @title The getDistributionNames method
# Argument
# aMtkExpFactors a mtkExpFactors object
# Return
# a list of string incluing all the distributionnames of the domains of mtkFactor objects
#
# move in mtkAllGenerics.R, uncomment otherwise
# setGeneric ("getDistributionNames", function(this) {standardGeneric("getDistributionNames")})
#
setMethod(f="getDistributionNominalValues",
		signature=c(this="mtkExpFactors"),
		definition=function(this="mtkExpFactors"){
			
			nValues <- sapply(this@expFactorsList,getDistributionNominalValue)
			
			return (nValues)
		}
)

setMethod(f="getDistributionNominalValueTypes",
		signature=c(this="mtkExpFactors"),
		definition=function(this="mtkExpFactors"){
			
			nTypes <- sapply(this@expFactorsList,getDistributionNominalValueType)
			
			return (nTypes)
		}
)
#
#===================================================================
# Method definition
#' The getDistributionParameters method
#' @param this an object of class \code{\linkS4class{mtkExpFactors}}
#' @return a list of lists of mtkParameter objects, with one sublist per factor
#' @exportMethod getDistributionParameters

# retourne une liste de listes, exemple:
#		list( x1=list(min=-pi,max=pi), x2=list(min=-pi,max=pi), x3=list(min=-pi,max=pi))

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getDistributionParameters", function(this) {standardGeneric("getDistributionParameters")})

setMethod(f="getDistributionParameters",
		signature=c(this="mtkExpFactors"),
		definition=function(this){
			
			distribParam <- lapply(this@expFactorsList,getDistributionParameters)
			return(distribParam)
			
		}
)


#
#===================================================================
# Method definition
#' The getDistributionParameters method
#' @param this an object of class \code{\linkS4class{mtkExpFactors}}
#' @return a list of lists of mtkParameter objects, with one sublist per factor
#' @exportMethod getDistributionParameters

# retourne une liste de listes, exemple:
#		list( x1=list(min=-pi,max=pi), x2=list(min=-pi,max=pi), x3=list(min=-pi,max=pi))

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getDistributionParameters", function(this) {standardGeneric("getDistributionParameters")})

setMethod(f="getFactorFeatures",
		signature=c(this="mtkExpFactors"),
		definition=function(this){
			fFeatures <- lapply(this@expFactorsList,getFeatures)
			return(fFeatures)
			
		}
)

#' The show method for an object of class \code{mtkFactor}
#' @param object an object of class \code{mtkFactor}
#' @value nil
setMethod(f="show", signature=c(object="mtkExpFactors"),
		definition=function(object){
			print(object)
			invisible()
		}
)


#' The show method for an object of class \code{mtkFactor}
#' @param object an object of class \code{mtkFactor}
#' @value nil
setMethod(f="print", signature=c(x="mtkExpFactors"), 
		definition=function(x,...){
			## general information 
			cat("------------------------------\n")
			cat('An object of class "mtkExpFactors" \n')
			cat("------------------------------\n")
			cat( 'Number of factors : ', length(x@expFactorsList), '\n')
			if(length(x@expFactorsList)>0)
				f <- lapply(x@expFactorsList,print)
			
			cat("------------------------------\n")
			invisible()
		}
)

#===================================================================
# Method definition
#' The mtkReadFactors method create a list of factors from a csv formated file (see specifications)  
#' @param file string for the name of the file to read
#' @param path string for the path of this file (default=getwd())
#' @return a list of mtkfactor for the slot expFactorsList
#' @exportMethod mtkReadFactors

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("mtkReadFactors", function(this, path) {standardGeneric("mtkReadFactors")})

#setMethod(f="mtkReadFactors", 
#          signature = c(file="character", path="character"),
#          definition = function(file, path)
#          {
#                        # is file existing ? warning : pathfile is platform dependant
#            
#            theFile <- paste(path, file, sep=.Platform$file.sep)
#            
#            if(!file.exists(theFile))
#              {
#                stop("stop in mtkReadFactors : file or path does not exist ! ")
#                }
#                
#            # global var
#            nbCol <- 8 # number of colomns of csv file
#            mySep1 <- ";"
#            mySep2 <- ","
#              
#            # how many lines ?
#            nbFact <-length(count.fields(theFile, sep = "\n"))#,comment.char = "#" ))
#            nbLines <- length(count.fields(theFile, sep = "\n", comment.char=""))
#            
#            factorList <- list()
#            
#            # parsing factors 
#            for (i in 1:nbLines - 1)  
#              {
#                # for each factor (i.e. line treatment )
#                runningFact <- scan(file=theFile, what = "",sep="\n",skip=i,nlines=1)
#                if (!grepl("^#", runningFact)) # skip comment (lines beginning by #)
#                {
#                  t1<-gsub("\t", "", runningFact)
#                  t2 <-strsplit(t1,mySep1)
#                  
#                  # factor name
#### HM, 7/11/11 apparently str_replace_all is obsolete, replaced by str_replace
#                  factName <- str_replace(sapply(t2,"[",1)," ","")
#                  # factor Id 
#                  factId <- str_replace(sapply(t2,"[",2)," ","")
#                   
#                  if (!is.na(charmatch(factId,table="",nomatch=NA)))
#                    factId <- factName
#                    
#                  # factor unit                  
#                  factUnit <-str_replace(sapply(t2,"[",3)," ","")
#                  if (!is.na(charmatch(factUnit,table="",nomatch=NA)))
#                    factUnit <- "UNIT"
#                  
#                  # factor domain
#                  # domain nominal value
#
#                  factNominalValue <- str_replace(sapply(t2,"[",4)," ","")
#                  if (!is.na (as.double(factNominalValue))) {
#                    factNominalValue <- as.double(factNominalValue)
#                    goodFactNominalVal <-mtkValue(factNominalValue, val=factNominalValue)
#                  }
#                  else { # nominalValue is a string
#                    nv <- "nominalValue"
#                    
#                    goodFactNominalVal <-mtkValue(nv,val=factNominalValue )
#                  }
#                  
#                 
#                  # domain distribution name
#                  # be careful : distribution must be a R fuction
#                  factDistributionName <- str_replace(sapply(t2,"[",5)," ","")
#                  
#                  # domain distribution parameters 
#                  # param1 = val1 , param2 = val2 , param3 = val3 , ... , paramN = valN
#                 factDistributionParam <- str_replace(sapply(t2,"[",6)," ","")
#                  splitFactDistributionParam<-strsplit(factDistributionParam,",")[1][[1]]
#                  nbParam <-length(str_locate_all(splitFactDistributionParam,"="))
#                  # for each param, name parameter is before =, value is after
#                  paramList <- list()
#                  
#                  if (nbParam > 0) {
#                    for (j in 1:nbParam)
#                    {
#                      paramTmp <- splitFactDistributionParam[j]
#                      nameParamTmp <- str_split(paramTmp,"=")[1][[1]][1]# before = 
#                      
#                      valParamTmp <- as.double(str_split(paramTmp,"=")[1][[1]][2]) # after =
#                      
#                      # 3 types : numeric, character, boolean
#                      
#                      valParamTmp <- str_split(paramTmp,"=")[1][[1]][2] # after =
#                      
#                      if(is.na(as.double(valParamTmp))) {
#                        trueLogical <-c("TRUE","True","true")
#                        falseLogical <-c("FALSE","False","false")
#                        
#                        if ((sum(valParamTmp == trueLogical)) > 0) {
#                          valParamTmp <-TRUE
#                        }
#                        else
#                          if ((sum(valParamTmp == falseLogical)) > 0) {
#                            valParamTmp <- FALSE
#                          }
#                          else { # string
#                            valParamTmp <- as.character(valParamTmp)
#                          }
#                      }
#                      else { # numeric
#                        valParamTmp <- as.double(valParamTmp)
#                      }
#                      
#                     
#                      mtkPtemp <- mtkParameter(name=nameParamTmp,val=valParamTmp)
#                      paramList <-c(paramList, mtkPtemp)
#                    }                  
#                  }
#                  
#                  # domain levels
#                  # level1 : weight1 ; level2 : weight2 ; level3 : weight3 ; ... ; levelN : weightN
#                  
#                  factLevels <- str_replace(sapply(t2,"[",7)," ","")
#                  splitFactLevels <- strsplit(factLevels,",")[1][[1]]
#                  nbLevels <-length(str_locate_all(splitFactLevels,":"))
#                  
#                  levelList <- list()
#				  weightList <- list()
#                  
#                  if (nbLevels > 0) {
#                    for (j in 1:nbLevels)
#                    {
#                      levelTmp <- splitFactLevels[j]
#                      valLevelTmp <- str_split(levelTmp,":")[1][[1]][1]# before : 
#                      if (is.na(as.double(str_split(levelTmp,":")[1][[1]][2]))) {
#                        weigthLevelTmp <-1/nbLevels
#                      }
#                      else {
#                        weigthLevelTmp <- as.double(str_split(levelTmp,":")[1][[1]][2]) # after ":"
#                        if (weigthLevelTmp == 1) {
#                          weigthLevelTmp <- 1/nbLevels
#                        }
#                      }
#                      
#                      levelList <-c(levelList, valLevelTmp)
#					  weightList <- c(weightList, weigthLevelTmp)
#                      }
#                  }
#				  levelList <- mtkLevels('categorical', levels=levelList, weights=weightList)
#                                    
#                  # features (same treatment as parameters)
#                  # feature1 = val1 , feature2 = val2 , feature3 = val3 , ... , featureN = valN
#                  userFeaturesList <- str_replace(sapply(t2,"[",8)," ","")
#                  splitUserFeaturesList<-strsplit(userFeaturesList,",")[1][[1]]
#                  nbFeatures <-length(str_locate_all(splitUserFeaturesList,"="))
#                  # for each feature, name feature is before =, value is after
#                  featureList <- list()
#                  
#                  if (nbFeatures > 0) {
#                    for (j in 1:nbFeatures)
#                    {
#                      featureTmp <- splitUserFeaturesList[j]
#                      nameFeatureTmp <- str_split(featureTmp,"=")[1][[1]][1]# before = 
#                      # 3 types : numeric, character, boolean
#                      
#                      valFeatureTmp <- str_split(featureTmp,"=")[1][[1]][2] # after =
#                      
#                      if(is.na(as.double(valFeatureTmp))) {
#                        trueLogical <-c("TRUE","True","true")
#                        falseLogical <-c("FALSE","False","false")
#                        if ((sum(valFeatureTmp == trueLogical)) > 0) {
#                          valFeatureTmp <-TRUE
#                        }
#                        else
#                          if ((sum(valFeatureTmp == falseLogical)) > 0) {
#                            valFeatureTmp <- FALSE
#                          }
#                          else { # string
#                            valFeatureTmp <- as.character(valFeatureTmp)
#                          }
#                      }
#                      else { # numeric
#                        valFeatureTmp <- as.double(valFeatureTmp)
#                      }
#                      
#                      mtkFtemp <- mtkFeature(nameFeatureTmp,valFeatureTmp )
#                      featureList <-c(featureList, mtkFtemp)
#                    }
#                  }
#                  
#                  # creating factor : 
#                  ## first, the domain
#                  ## slot distributionName a string representing the  distribution law
#                  ## slot nominalValue mtkValue corresponding to nominal value
#                  ## slot levelList list of mtkLevel objects
#                  ## slot distributionParameters list of mtkParameter
#                  runningFactDomain <-mtkDomain(distributionName=factDistributionName,domainNominalValue=goodFactNominalVal)
#                  if (length(paramList) > 0)
#                  {
#                    setDistributionParameters(runningFactDomain,paramList)
#                  }
#                  if (length(levelList) > 0)
#                  {
#                    setLevels(this=runningFactDomain,levels=levelList)
#                  }
#                  
#                  # second, factor's attributes
#                  ## slot name string, name of the mtkFactor
#                  ## slot id string, single Id of the mtkFactor
#                  ## slot unit string, optional
#                  ## slot domain mtkDomain object for specify the mtkFactor domain
#                  
#                  goodFactor <- mtkFactor(name=factName, id=factId, unit=factUnit, domain=runningFactDomain )
#                 
#                  ## slot featureList, list of mtkFeature objects for specific characteristic of a mtkFactor
#                 if (length(featureList) > 0)
#                  {
#                    setFeatures(this=goodFactor,aFList=featureList) 
#                  }
#                  # Add new factor
#                  #print(goodFactor)
#                  #scan()
#                  factorList <- c(factorList,goodFactor)
#                  
#                }# end of skip if comment
#                
#                #browser()
#              }# end of creating one factor
#              
#              return(mtkExpFactors(expFactorsList=factorList))
#              
#        
#
#            }
### HM, 05/06/2011 : Mis en commentaires 
###,
###          where = topenv(parent.frame()),
###          valueClass = NULL, sealed = FALSE)
#)



#' Extract or Replace Parts of an object of class \code{mtkExpFactors}
#' @param x an object of class \code{mtkExpFactors}
#' @param i an index specifying the element to extract
#' @return an object of class \code{mtkFactor}
setMethod(f="[[", signature=c(x="mtkExpFactors", i="ANY"),
          definition=function(x,i){
            out <- (x@expFactorsList)[[i]]
            out
          })

#' Extract or Replace Parts of an object of class \code{mtkExpFactors}
#' @param x an object of class \code{mtkExpFactors}
#' @param i an index specifying elements to extract or replace
#' @return an object of class \code{mtkFactor}
setMethod(f="[", signature=c(x="mtkExpFactors", i="ANY"),
          definition=function(x,i){
            x@expFactorsList <- (x@expFactorsList)[i]
            x
          })

#' Extract or Replace Parts of an object of class \code{mtkExpFactors}
#' @param x an object of class \code{mtkExpFactors}
#' @return an object of class \code{mtkExpFactors}
setMethod(f="$", signature=c(x="mtkExpFactors"),
          definition=function(x,name){
            out <- (x@expFactorsList)[[name]]
            out
          })

##' Extract or Replace Parts of an object of class \code{mtkExpFactors}
##' @param x an object of class \code{mtkExpFactors}
##' @param i an index specifying elements to extract or replace
##' @param value an object of class \code{mtkFactor}
##' @return an object of class \code{mtkExpFactors}
#setMethod(f="[<-", signature=c(x="mtkExpFactors", i="numeric", j="missing", value="mtkFactor"),
#          definition=function(x,i,value){
#            ## a loop is necessary to keep the factor names ok
#            li <- length(i)
#            for(k in seq(li)){
#              value.k <- value
#              value.k@name <- names(x@expFactorsList)[i[k]]
#              value.k@id <- names(x@expFactorsList)[i[k]]
#              x@expFactorsList[i[k]] <- value.k
#            }
#            x
#          })

#' Extract or Replace Parts of an object of class \code{mtkExpFactors}
#' @param x an object of class \code{mtkExpFactors}
#' @param i an index specifying elements to extract or replace
#' @param value an object of class \code{mtkFactor}
#' @return an object of class \code{mtkExpFactors}
#setMethod(f="[<-", signature=c(x="mtkExpFactors", i="character", j="missing", value="mtkFactor"),
#          definition=function(x,i,value){
#            ## a loop is necessary to keep the factor names ok
#            li <- length(i)
#            for(k in seq(li)){
#              value.k <- value
#              value.k@name <- i[k]
#              value.k@id <- i[k]
#              exists <- which(getNames(x) == i[k])
#              if(length(exists) > 0){
#                x@expFactorsList[exists] <- value.k
#              }
#              else{
#                xnames <- names(x@expFactorsList)
#                x@expFactorsList <- c(x@expFactorsList, value.k)
#                names(x@expFactorsList) <- c(xnames, i[k])
#              }
#            }
#            x
#          })
#' Extract or Replace Parts of an object of class \code{mtkExpFactors}
#' @param x an object of class \code{mtkExpFactors}
#' @param value an object of class \code{mtkFactor}
#' @return an object of class \code{mtkExpFactors}
#setMethod("$<-", signature(x = "mtkExpFactors", value = "mtkFactor"),
#    function (x, name, value) 
#    {
#            ## a loop is necessary to keep the factor names ok
#            li <- length(name)
#            for(k in seq(li)){
#              value.k <- value
#              value.k@name <- name[k]
#              value.k@id <- name[k]
#              exists <- which(getNames(x) == name[k])
#              if(length(exists) > 0){
#                x@expFactorsList[exists] <- value.k
#              }
#              else{
#                xnames <- names(x@expFactorsList)
#                x@expFactorsList <- c(x@expFactorsList, value.k)
#                names(x@expFactorsList) <- c(xnames, name[k])
#              }
#            }
#            x
#    }
#)
##================================================================================
## HM, 05/11/11 (END)
## extract and replace method and other basic methods for the R user 
##================================================================================


# Mexico Toolkit
#
# version 	: 0.03
# date		: 30 nov 2009
# MAJ   	: 05 June 2011
# licence	: GPL

# Author(s) : Hervé Monod, Hervé Richard, Juhui Wang, MIA INRA
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 291                  $: revision number of the last spread
#' $Author:: jwang             $: author of the last spread
#' $Date:: 2012-06-14 11:23:30#$: date of the last spread

#-----------------------------------------------------------------------
#' Nota :
#
### R package : mexico toolkit
### file : Mexico Toolkit AllGenerics methods
### describe all generic methods used in Mexico Toolkit Package
###
#=======================================================================
## General methods
#==================================
### run
##
#' Launchs a task defined in a process or workflow.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkParsor}}], [\code{\linkS4class{mtkExpWorkflow}}], [\code{\linkS4class{mtkProcess}}] and its sub-classes .
#' Examples of "run" are:
#' \itemize{
#' \item{run(this, context)}{this is an object of class [\code{\linkS4class{mtkDesigner}}], and context is an object of class [\code{\linkS4class{mtkExpWorkflow}}].}
#' \item{run(this, context)}{this is an object of class [\code{\linkS4class{mtkEvaluator}}], and context is an object of class [\code{\linkS4class{mtkExpWorkflow}}].}
#' }
#' 
#' @param this an object corresponding to the task to launch
#' @param context an object specifying the context which manages the task
#' @title The generic function run
setGeneric ("run", function(this, context) {standardGeneric("run")})

###
### report

#' Reports the results produced by a process or workflow.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}], [\code{\linkS4class{mtkExpWorkflow}}].
#' @param this a process that produces  data to report
#' @title The generic function report

setGeneric ("report", function(this) {standardGeneric("report")})



### A supprimer car redondant avec make.mtkParameterList
#### mtkParameters
##' mtkParameters method to create a list of mtkParameters
##' @title The generic function mtkParameters
##' @param listParam list of parameters to convert in list of mtkParameters
#setGeneric ("mtkParameters", function(listParam) {standardGeneric("mtkParameters")})

###
### make.mtkFactor
##' make.mtkFactor method to create a factor interactively
##' @title The method make.mtkFactor
##' @param name a character string, giving the name of the factor
##' @param id an optional character string, giving the code name of the factor
##' @param unit an optional character string, giving the units of the factor values (to be used in graphics, for example)
##' @param nominal an optional nominal value, character or numeric depending on the distribution
##' @param distribName a character string, giving the R or mtk name of the distribution
##' @param distribPara an optional named list of distribution parameters
##' @param levels an optional character or numeric vector of the factor levels
##' @param weights an optional vector of levels weights
##' @param features optional list of factor features

##setGeneric ("make.mtkFactor", function(name, id, unit, nominal, distribName, distribPara, levels, weights, features) {standardGeneric("make.mtkFactor")})


### make.mtkLevelList
##' make.mtkLevelList method to create a list of mtkLevel
##' @title The generic function make.mtkLevelList
##' @param x a numeric or character vector 
##' @param weights an optional numeric vector

## setGeneric ("make.mtkLevelList", function(x, weights) {standardGeneric("make.mtkLevelList")})


### make.mtkParameterList
##' make.mtkParameterList method to create a list of mtkLevel
##' @title The generic function make.mtkParameterList
##' @param list a numeric or character vector 

## setGeneric ("make.mtkParameterList", function(list) {standardGeneric("make.mtkParameterList")})

###
### make.mtkFeatureList
##' make.mtkFeatureList method to create a list of mtkFeature
##' @title The generic function make.mtkFeatureList
##' @param list a numeric or character vector 

## setGeneric ("make.mtkFeatureList", function(x) {standardGeneric("make.mtkFeatureList")})

########## HM, 2011-11-03 (FIN) ##########


### commenter (en attente de suppression) de la déclaration de mtkFeatures car redondant avec make.mtkFeatureList
#### mtkFeatures
##' mtkFeatures method to create a list of mtkFeatures
##' @title The generic function mtkFeatures
##' @param listParam list of parameters to convert in list of mtkFeatures
#
#setGeneric ("mtkFeatures", function(listParam) {standardGeneric("mtkFeatures")})

###
# mtkReadFactors
#' mtkReadFactors method create a list of factors from a csv formated file (see specifications) 
#' @title Generic method to read factors from a csv-like file
#' @param file string for the name of the file to read
#' @param path string for the path of this file (default=getwd())
#' @return a list of mtkfactor for the slot expFactorsList

setGeneric ("mtkReadFactors", function(file, path) {standardGeneric("mtkReadFactors")})

###

#==================================
## get.XXX methods
#==================================
### getName

#' Gets the name of an object.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkValue}}], [\code{\linkS4class{mtkFactor}}].
#' @title The generic function getName
#' @param this an object with a slot "name"
#' @return the "name"  of the underlying object

setGeneric ("getName", function(this) {standardGeneric("getName")} )

###
### getValue
#' Gives the typed value of an object managing a triple (name, valueType, value).
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkValue}}].
#' @title The generic function getValue
#' @param this the underlying object.
#' @return the value with the relevant type
setGeneric ("getValue", function(this) {standardGeneric("getValue")} )

###
### getValueType
#' Gives the type of a mtkValue object with the triple (name, valueType, value).
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkValue}}].
#' @title The generic function getValueType
#' @param this an object of class mtkValue
#' @return the type of the value represented by \code{this}
setGeneric ("getType", function(this) {standardGeneric("getType")} )

###
### getNominalValueType
#' Gives the type of the nominal value of an object related to a probability distribution. 
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkDomain}}].
#' @title The generic function getNominalValueType
#' @param this an object related to a probability distribution
#' @return the value type of the nominalValue of the object
setGeneric ("getNominalValue", function(this) {standardGeneric("getNominalValue")} )
setGeneric ("getDistributionNominalValue", function(this) {standardGeneric("getDistributionNominalValue")} )
setGeneric ("getDistributionNominalValues", function(this) {standardGeneric("getDistributionNominalValues")} )

setGeneric ("getNominalValueType", function(this) {standardGeneric("getNominalValueType")} )
setGeneric ("getDistributionNominalValueType", function(this) {standardGeneric("getDistributionNominalValueType")} )
setGeneric ("getDistributionNominalValueTypes", function(this) {standardGeneric("getDistributionNominalValueTypes")} )
###
### getNames
#' Gives the names of the components of an object 
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkExpFactors}}].
#' @param this an object with a slot 'name'
#' @return a list of strings corresponding to the names of the object.
#' @title The generic function getNames

setGeneric ("getNames", function(this) {standardGeneric("getNames")})
setGeneric ("getFactorNames", function(this) {standardGeneric("getFactorNames")})
###
# getDistributionName
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkDomain}}], [\code{\linkS4class{mtkFactor}}].
#' @param this an object with a slot 'distributionName' or containing an object with a slot 'distributionName'
#' @return a string giving the probability distribution name of an object
#' @title The generic function getDistributionName.

setGeneric ("getDistributionName", function(this) {standardGeneric("getDistributionName")})
setGeneric ("getDiscreteDistributionType", function(this) {standardGeneric("getDiscreteDistributionType")})
###
# getDistributionNames
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkExpFactors}}].
#' @param this an object containing objects with a slot 'distributionName'
#' @return a list of strings the probability distribution names of an object
#' @title The generic function getDistributionNames.

setGeneric ("getDistributionNames", function(this) {standardGeneric("getDistributionNames")})

###
### getLevels
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkDomain}}].
#' @param this an object with a slot 'mtkLevelList'
#' @return a list of mtkLevel objects.
#' @title The generic function getLevels.

setGeneric ("getLevels", function(this) {standardGeneric("getLevels")})
setGeneric ("getDiscreteDistributionLevels", function(this) {standardGeneric("getDiscreteDistributionLevels")})
###
### getWeights
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkDomain}}].
#' @param this an object with a slot 'mtkLevelList'
#' @return a list of mtkLevel objects.
#' @title The generic function getLevels.

setGeneric ("getWeights", function(this) {standardGeneric("getWeights")})
setGeneric ("getDiscreteDistributionWeights", function(this) {standardGeneric("getDiscreteDistributionWeights")})
###
### getDistributionParameters

#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkDomain}}], [\code{\linkS4class{mtkFactor}}], [\code{\linkS4class{mtkExpFactors}}].
#' @param this an object related to probability distributions
#' @return a list of mtkParameter objects or a list of lists of mtkParameter objects
#' @title The generic function getDistributionParameters.

setGeneric ("getDistributionParameters", function(this) {standardGeneric("getDistributionParameters")})

###
###getFeatures

#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkFactor}}].
#' @param this an object with a slot 'feature'
#' @return a list of mtkFeature objects
#' @title The generic function getFeatures.

setGeneric ("getFeatures", function(this) {standardGeneric("getFeatures")})

setGeneric ("getMTKFeatures", function(this) {standardGeneric("getMTKFeatures")})

setGeneric ("getFactorFeatures", function(this) {standardGeneric("getFactorFeatures")})
###

### getDomain

#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkFactor}}].
#' @param this an object with a 'domain' slot
#' @return an object of class mtkDomain
#' @title The generic function getDomain.

setGeneric ("getDomain", function(this) {standardGeneric("getDomain")})
setGeneric ("setDomain", function(this, domain) {standardGeneric("setDomain")})
###
### getParameters
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}] and its inherited classes.
#' @param this an object with a slot 'parameter'
#' @return a named list
#' @title The generic function getParameters.
#' @note The vector of parameters are  converted as  a named list:
#' \code{(name of parameter1 = value of parameter1, name of parameter2 = value of parameter 2, ...)}.

setGeneric ("getParameters",         function(this)        {standardGeneric("getParameters")} )

###
### getResult
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}] and its inherited classes.
#' @param this an object of class [\code{\linkS4class{mtkProcess}}.
#' @return an object of class [\code{\linkS4class{mtkResult}}]
#' @title The generic function getResult.

setGeneric ("getResult",         function(this)        {standardGeneric("getResult") })

###
### getData
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}] and its inherited classes.
#' @param this an object of class [\code{\linkS4class{mtkProcess}}.
#' @return a data frame corresponding to the results produced by the process.
#' @title The generic function getData.

setGeneric ("getData",         function(this)        {standardGeneric("getData")} )


###
### extractData
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkExpWorkflow}}] and its inherited classes.
#' @param this an object of class [\code{\linkS4class{mtkExpWorkflow}}.
#' @param name a vector of string from "design", "evaluate", or "analyze"  indicating the results to extract.
#' @return a data frame corresponding to the results produced by the process.
#' @title The generic function extractData.

setGeneric ("extractData",         function(this, name)        {standardGeneric("extractData")})


###
### getProcess
#' gets search a process of the workflow by name.
#' This function is implemented in the following classe:
#' [\code{\linkS4class{mtkExpWorkflow}}].
#' @param this an object of class [\code{\linkS4class{mtkExpWorkflow}}].
#' @param name  a string from "design", "evaluate", or "analyze"  indicating the  process to get.
#' @title The generic function getProcess.

setGeneric ("getProcess",         function(this, name)        {standardGeneric("getProcess")})


###
#==================================
### set.XXX methods
#==================================
###
### setProcess
#' Change a process of the workflow to a new one.
#' This function is implemented in the following classe:
#' [\code{\linkS4class{mtkExpWorkflow}}].
#' @param this an object of class [\code{\linkS4class{mtkExpWorkflow}}].
#' @param p an object of class [\code{\linkS4class{mtkProcess}}].
#' @param name a string from (design, evaluate, analyze) to associate with the underlying process.
#' @title The generic function setProcess.


setGeneric ("setProcess", function(this, p, name) {standardGeneric("setProcess")})

### setXMLFilePath
#' Changes the XML file to parse.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkParsor}}]
#' @param this the parsor.
#' @param xmlPath the XML file with path.
#' @title The generic function setXMLFilePath.

setGeneric ("setXMLFilePath", function (this, xmlPath) {standardGeneric("setXMLFilePath")})
###
### setReady
#' Makes the process ready to run.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}], [\code{\linkS4class{mtkExpWorkflow}}]
#' @param this a process or workflow.
#' @param switch a logical to tell if the process is ready to run.
#' @title The generic function setReady.
setGeneric ("setReady", function(this, switch) {standardGeneric("setReady")} )

### setState
#' Marks the process as aleady run.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}].
#' @param this a process.
#' @param switch a logical to tell if the process is finished.
#' @title The generic function setState.
setGeneric ("setState", function(this, state) {standardGeneric("setState")} )

###
### setParameters
#' Assigns a new vector of parameters to the process.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}].
#' @param this the underlying process.
#' @param f a vector of objects from class [\code{\linkS4class{mtkParameter}}]
#' @title The generic function getParameters.

setGeneric ("setParameters", function(this, f) {standardGeneric("setParameters")} )

###
### serializeOn
#' Serialize the process.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}] and its sub-classes.
#' @param this the underlying process.
#' @param f a vector of objects from class [\code{\linkS4class{mtkParameter}}]
#' @title The generic function serializeOn.

setGeneric ("serializeOn", function(this) {standardGeneric("serializeOn")} )

###
### setName
#' Assigns  a processing step (design, evaluate, anslysis) to a process.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}].
#' @param this the underlying process.
#' @param name a string from (design, evaluate, analyze).
#' @title The generic function setName.
setGeneric ("setName", function(this, name) {standardGeneric("setName")} )
###
### setLevels
#' Assigns a new slot 'levelList' to an object related to a discrete probability distribution
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkDomain}}].
#' @param this an object with a slot 'levelList'
#' @title The generic function setlevelList.

setGeneric (name="setLevels", def=function(this, levels) {standardGeneric("setLevels")})
### setWeights
#' Assigns a new slot 'levelList' to an object related to a discrete probability distribution
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkDomain}}].
#' @param this an object with a slot 'levelList'
#' @title The generic function setlevelList.

setGeneric (name="setWeights", def=function(this, weights) {standardGeneric("setWeights")})

###
### setDistributionParameters

#' Assigns a new slot 'distributionParameters' to an object related to a probability distribution
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkDomain}}].
#' @param this an object with a slot 'distributionParameters'
#' @param aDistParamList a list of mtkParameter objects
#' @title The generic function setDistributionParameters.

setGeneric ("setDistributionParameters",   function(this, aDistParamList) {standardGeneric("setDistributionParameters")})

###
### setFeatures matchType

#' Assigns a new slot 'feature' to an object related to factors
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkFactor}}].
#' @param this an object with a slot 'feature'
#' @param aFList a list of mtkFeature objects
#' @title The generic function setFeatures
#signature=c(this="mtkFactor",aFList="list"),

setGeneric (name="setFeatures", def=function(this, aFList) {standardGeneric("setFeatures")})
###
### setFactors
#' Assigns a list of mtkFactor objects to the object mtkExpFactors
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkExpFactors}}].
#' @param this an object containing information about factors
#' @param aFactList a list of mtkFactor objects
#' @title The generic function setFactors.

setGeneric ("setFactors", function(this, aFactList) {standardGeneric("setFactors")})
setGeneric ("getFactors", function(this) {standardGeneric("getFactors")})
###
### setType
#' Assigns  a data type to a factor.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkFactor}}].
#' @param this an object of class [\code{\linkS4class{mtkFactor}}]
#' @param type a string corresponding to a data type.
#' @title The generic function setType.

setGeneric ("setType", function(this, type) {standardGeneric("setType")})
setGeneric ("setValue", function(this, val) {standardGeneric("setValue")})
###
#==================================
### add.XXX or delete.XXX methods
#==================================
###
### addProcess
#' Adds a process to the workflow.
#' This function is implemented in the following classe:
#' [\code{\linkS4class{mtkExpWorkflow}}].
#' @param this an object of class [\code{\linkS4class{mtkExpWorkflow}}].
#' @param p an object of class [\code{\linkS4class{mtkProcess}}].
#' @param name a string from (design, evaluate, analyze) to associate with the underlying process.
#' @title The generic function addProcess.

setGeneric ("addProcess", function(this, p, name) {standardGeneric("addProcess")})

###
### deleteProcess
#' Delete a process from the workflow.
#' This function is implemented in the following classe:
#' [\code{\linkS4class{mtkExpWorkflow}}].
#' @param this an object of class [\code{\linkS4class{mtkExpWorkflow}}].
#' @param name a string from "design", "evaluate" or "analyze", indicating which process to remove.
#' @title The generic function deleteProcess.

setGeneric ("deleteProcess", function(this, name) {standardGeneric("deleteProcess")} )

###
### reevaluate
#' evaluate the context of the processes of the workflow to know if we can re-use the results produced before 
#' or we  have to re-run them.
#' This function is implemented in the following classe:
#' [\code{\linkS4class{mtkExpWorkflow}}] and its sub-classes.
#' @param this an object of class [\code{\linkS4class{mtkExpWorkflow}}].
#' @param name a string from "design", "evaluate" or "analyze", indicating frow which process we
#' do the re-evaluate.
#' @title The generic function reevaluate.

setGeneric ("reevaluate", function(this, name) {standardGeneric("reevaluate")} )

###
#==================================
### is.XXX methods
#==================================

### is.ready
#' Return TRUE if a task (process or workflow) is ready to run.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}] and [\code{\linkS4class{mtkExpWorkflow}}]
#' @param this an object related to a process or a workflow
#' @return TRUE or FALSE.
#' @title The generic function is.ready.
setGeneric ("is.ready", function(this) {standardGeneric("is.ready")} )

###
### is.finished
#' Return TRUE if a task is finished run.
#' This function is implemented in the following classes:
#' [\code{\linkS4class{mtkProcess}}]
#' @param this an object related to a process
#' @return TRUE or FALSE.
#' @title The generic function is.finished.
setGeneric ("is.finished", function(this) {standardGeneric("is.finished")} )
###
###

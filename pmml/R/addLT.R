# PMML: Predictive Model Markup Language
#
# Copyright (c) 2009-2013, some parts by Togaware Pty Ltd and other by Zementis, Inc. 
#
# This file is part of the PMML package for R.
#
# The PMML package is free software: you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 2 of 
# the License, or (at your option) any later version.
#
# The PMML package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please see the
# GNU General Public License for details (http://www.gnu.org/licenses/).
######################################################################################
#
# Author: Tridivesh Jena
# Date: Aug 2013
#
#
# Possible future functions:
# Node -> InternalNode
# nodeStr <- toString.XMLNode(Node)
# internalDocument <- xmlTreeParse(nodeStr,asText=TRUE,useInternalNodes=TRUE)
# internalNode <- getNodeSet(internalDocument,"/")[[1]]
#
# InternalNode -> Node
# internalNodeStr <- toString.XMLNode(internalNode)
# internalNodeDocument <- xmlTreeParse(internalNodeStr,asText=TRUE,useInternalNodes=FALSE)
# node <- internalNodeDocument$doc$children[[1]]
# XMLInternalDocument -> XMLNode
# xmlTreeParse(toString.XMLNode(modelInternalDocument),useInternalNodes=FALSE)$doc$children[[1]]
#
#######################################################################
#
# Brief summary of useful commands from the XML package...most are needed due to the 
# unecessarily large variety of XML storage formats
#
# getNodeSet(XMLInternalDocument) -> list(XMLInternalNode)
# xmlTreeParse(file, useInternalNodes=FALSE) -> XMLDocument
# xmlTreeParse(file, useInternalNodes=TRUE) -> XMLInternalDocument 
# from 2 steps above: XMLDocument$doc$children -> XMLNode
# string <- toString.XMLNode(XMLNode)  or XMLInternalNode
# xmlTreeParse(string,asText=TRUE,useInternalNodes=TRUE/FALSE) same as before
# addChildren(XMLNode,XMLNode) -> new XMLNode
# addChildren(XMLInternalNode, XMLInternalNode, at=integer) -> same previous XMLInternalNode, modified
# xmlParent(XMLInternalNode) -> XMLInternalNode
# pmml is XMLNode, saveXML can save XMLInternalNode as XMLNode
#
# getNodeSet(XMLInternalDocument,xpath,namespace)
# if namespace is p=http://www.dmg.org/PMML-4_1  then Model is p:Model
# Eg: predicted MiningField from PMML model is 
# ans <- getNodeSet(pmml,"/p:PMML/*/p:MiningSchema/child::*[attribute::name='predictedScore']",
#                   c(p="http://www.dmg.org/PMML-4_1"))
# Note we have to start from "/"   not "." 
# Once we extract an element from PMML using getNodeSet, namespace is "" and the next application
# doesnt require a namespace 
#

addLT <- function(xmlmodel=NULL, transforms=NULL, namespace="4_2",...) {
warning("Deprecated; use pmml() with 'transforms' argument instead")
#to avoid malloc error
flush.console()
 if(!is.null(namespace))
 {
  namespace <- .getNamespace(namespace)
 } else
 {
  namespace <- .getNamespace("4_2")
 }

# we expect the input to be always a XML Node
# If a file is to be used, the function 'fileToXMLNode' should be used to convert the file to a XMLNode object

  modelstring <- toString.XMLNode(xmlmodel)
  modelInternalDocument <- xmlTreeParse(modelstring,asText=TRUE,useInternalNodes=TRUE)

# if no transforms, return xml model
# if input was filename, read it in, convert to xml, return 
  if(is.null(transforms))
  {
   return(xmlmodel)  
  }


# transforms is not null
# first read in the transform
  transformstring <- toString.XMLNode(transforms)
  transformInternalDocument <- xmlTreeParse(transformstring,asText=TRUE,useInternalNodes=TRUE)

# Get DerivedField nodes from input transforms node
  dfNodeSet <- getNodeSet(transformInternalDocument,"/child::*/DerivedField")


# DataDictionary and MiningSchema to be modified later
  dd <- getNodeSet(modelInternalDocument,"/p:PMML/p:DataDictionary",c(p=namespace))[[1]]
  ms <- getNodeSet(modelInternalDocument,"/p:PMML/child::*/p:MiningSchema",c(p=namespace))[[1]]
  mfn <- getNodeSet(modelInternalDocument,"/p:PMML/child::*/p:MiningSchema/p:MiningField/@name",c(p=namespace))
  dfn <- getNodeSet(modelInternalDocument,"/p:PMML/p:DataDictionary/p:DataField/@name",c(p=namespace))
# mfn is the same as formulaFields but I kept it here to be symmetric. We can delete and replace its occurances
# if memory becomes an issue

  formulaFields <- getNodeSet(modelInternalDocument,"/p:PMML/child::*/p:MiningSchema/p:MiningField/@name",c(p=namespace))
  derivedFields <- getNodeSet(transformInternalDocument,"/LocalTransformations/DerivedField/@name")
  originalFields <- getNodeSet(transformInternalDocument,"/LocalTransformations/DerivedField/descendant::*/@field")

  formulaString <- lapply(formulaFields,FUN=toString)
  derivedString <- lapply(derivedFields,FUN=toString)
  originalString <- lapply(originalFields,FUN=toString)

  od <- setdiff(originalString,derivedString)
  fd <- setdiff(formulaString,derivedString)
  nf <- union(od,fd)

  delete <- which(!(formulaFields %in% nf))
  if(length(delete) > 0)
  {
   for(i in 1:length(delete))
   {
    ms <- removeChildren(ms,delete[i]-i+1) 
   }
  }

  delete <- which(!(dfn %in% nf))
  if(length(delete) > 0)
  {
   for(i in 1:length(delete))
   {
    dd <- removeChildren(dd,delete[i]-i+1)
   } 
  }

  continuous <- c("Discretize","NormContinuous")
  discrete <- c("NormDiscrete","FieldColumnPair")

  for(i in which(!(nf %in% formulaFields)))
  {
    dfPath<-paste("/LocalTransformations/DerivedField[descendant::*[@field='",nf[i],"']]/child::*",sep="")
    ddfName <- getNodeSet(transformInternalDocument,dfPath)
    xformType <- lapply(ddfName,xmlName)
 
    if(any(xformType %in% continuous))
    {
      newmf <- newXMLNode("MiningField",attrs=c(name=paste(nf[[i]],sep=""),usageType="active",optype="continuous"))
      newdf <- newXMLNode("DataField",attrs=c(name=paste(nf[[i]],sep=""),optype="continuous", dataType="double"))
      ms <- addChildren(ms,newmf)
      dd <- addChildren(dd,newdf)
    } else 
    {
      newmf <- newXMLNode("MiningField",attrs=c(name=paste(nf[[i]],sep=""),usageType="active",optype="categorical"))
      newdf <- newXMLNode("DataField",attrs=c(name=paste(nf[[i]],sep=""),optype="categorical",dataType="string"))
      ms <- addChildren(ms,newmf)
      dd <- addChildren(dd,newdf)
    }
   }


# Does input PMML have LocalTransformations already?
  ltNodeSet <- getNodeSet(modelInternalDocument,"/p:PMML/child::*/p:LocalTransformations",c(p=namespace))

# if not, just add the whole input LocalTransformations element
# else add as last children in LT
  position <- 0
  if(is.null(ltNodeSet) || length(ltNodeSet)==0)
  {
# go through the expected sequence of elements; find the last element from the sequence
# (parent model)->Extension->MiningSchema->Output->ModelStats->ModelExplanation->Targets
    ms <- getNodeSet(modelInternalDocument,"/p:PMML/*/p:MiningSchema",c(p=namespace))[[1]]
# modelElement = getNodeSet(modelInternalDocument,"/p:PMML/*[child::p:MiningSchema]",c(p=namespace))[[1]]

    if(is.null(ms) || length(ms)==0)
      stop("MiningSchema element is required.")

    mdl<-getNodeSet(modelInternalDocument,"/p:PMML/*[child::p:MiningSchema]",c(p=namespace))[[1]]
    modelElementList <- toString.XMLNode(getNodeSet(mdl,"./*"))

    if(grepl("Extension",modelElementList))
      position <- position + 1

    # can skip this, but not...just to be nicer looking!
    if(grepl("MiningSchema",modelElementList))
      position <- position + 1

    if(grepl("Output",modelElementList))
      position <- position + 1

    if(grepl("ModelStats",modelElementList))
      position <- position + 1

    if(grepl("ModelExplanation",modelElementList))
      position <- position + 1

    if(grepl("Targets",modelElementList))
      position <- position + 1

    ltNode <- getNodeSet(transformInternalDocument,"/*")[[1]]
    addChildren(mdl,ltNode,at=position)
  } else
  {
    addChildren(ltNodeSet[[1]],dfNodeSet,at=0)
  }

  ddfPath <- paste("/p:PMML/p:DataDictionary/p:DataField",sep="")
  numdf <- getNodeSet(modelInternalDocument,ddfPath,c(p=namespace))
  if(length(numdf) != 0)
  {
    dd <- addAttributes(dd,numberOfFields=length(numdf))
  }
  return(modelInternalDocument) 
}


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
# Date: May 2014
#
#
#######################################################################

addMSAttributes <- function(xmlmodel=NULL, attributes=NULL, namespace="4_2", ...) {

#to avoid malloc error
 flush.console()
# if(!is.null(namespace))
# {
  namespace <- .getNamespace(namespace)
# } else
# {
#  namespace <- .getNamespace("4_2")
# }

 if(!is.data.frame(attributes))
 {
  print("Please give attribute information in data frame format.")
 }

# we expect the input to be always a XML Node
# If a file is to be used, the function 'fileToXMLNode' should be used to convert the file to a XMLNode object

  modelstring <- toString.XMLNode(xmlmodel)
  modelInternalDocument <- xmlTreeParse(modelstring,asText=TRUE,useInternalNodes=TRUE)

# DataDictionary and MiningSchema to be modified later
  ms <- getNodeSet(modelInternalDocument,"/p:PMML/child::*/p:MiningSchema",c(p=namespace))[[1]]
  formulaFields <- getNodeSet(modelInternalDocument,"/p:PMML/child::*/p:MiningSchema/p:MiningField/@name",c(p=namespace))
  formulaString <- lapply(formulaFields,FUN=toString)

# remove levels from data frame (strings stored as factors) so that new values can be added 
  for(k in 1:ncol(attributes)){attributes[[k]]<-as.character(attributes[[k]])}
# or use attributes <- read.table(filename,stringsAsFactors=FALSE)  
  fieldNames <- colnames(attributes)
  internalAttr <- attributes
  if(FALSE %in% is.element(fieldNames,formulaFields))
  {
    print(cat(fieldNames[which(!(is.element(fieldNames,formulaFields)))],"not in model",sep=""))
    for(k in 1:length(fieldNames[which(!(is.element(fieldNames,formulaFields)))]))
    {
      internalAttr[,fieldNames[which(!(is.element(fieldNames,formulaFields)))][k]]<-NULL
    }
  }

  for(i in 1:ncol(internalAttr))
  {
    mf <- getNodeSet(ms,paste0("p:MiningField[@name='",colnames(internalAttr)[i],"']"),c(p=namespace))[[1]]
    s = structure(names=as.character(rownames(internalAttr)[which(!is.na(internalAttr[,i]))]),internalAttr[which(!is.na(internalAttr[,i])),i])
#    for(j in 1:nrows(internalAttr))
#    {
    addAttributes(mf,.attrs=s,append=TRUE)

#      addAttributes(mf,rownames(internalAttr)[j]=internalAttr[j,i],append=FALSE)
#    }
  }
 
  return(modelInternalDocument) 
}


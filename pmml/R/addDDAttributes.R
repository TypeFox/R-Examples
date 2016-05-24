# PMML: Predictive Model Markup Language
#
# Copyright (c) 2009-2015, some parts by Togaware Pty Ltd and other by Zementis, Inc. 
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

addDDAttributes <- function(xmlmodel=NULL, attributes=NULL, field=NULL, namespace="4_2", ...) {

#to avoid malloc error
 flush.console()
 namespace <- .getNamespace(namespace)

 if(!is.data.frame(attributes) && !is.vector(attributes))
 {
  print("Please give attribute information in data frame format or as a vector.")
 }

# we expect the input to be always a XML Node
# If a file is to be used, the function 'fileToXMLNode' should be used to convert the file to a XMLNode object

  modelstring <- toString.XMLNode(xmlmodel)
  modelInternalDocument <- xmlTreeParse(modelstring,asText=TRUE,useInternalNodes=TRUE)
# DataDictionary and MiningSchema to be modified later
  dd <- getNodeSet(modelInternalDocument,"/p:PMML/p:DataDictionary",c(p=namespace))[[1]]
  formulaFields <- getNodeSet(modelInternalDocument,"/p:PMML/p:DataDictionary/p:DataField/@name",c(p=namespace))
  formulaString <- lapply(formulaFields,FUN=toString)

# convert field information from vector to data frame format
  if(is.vector(attributes))
  {
   attributes <- data.frame(attributes)
   colnames(attributes) <- field 
   fieldNames <- field
  } else
   fieldNames <- colnames(attributes)

# remove levels from data frame (strings stored as factors) so that new values can be added 
  attributes[] <- lapply(attributes,as.character) 
# or use attributes <- read.table(filename,stringsAsFactors=FALSE)

  if(!all(is.element(fieldNames,formulaFields)))
  {
    print("WARNING: the following field additions will be ignored")
    cat(fieldNames[which(!(is.element(fieldNames,formulaFields)))],"not in model\n")
    attributes[,c(which(!(is.element(fieldNames,formulaFields))))] <- NA
  }

  for(i in 1:ncol(attributes))
  {
   if(!all(is.na(attributes[,i])))
   {
    df <- getNodeSet(dd,paste0("p:DataField[@name='",colnames(attributes)[i],"']"),c(p=namespace))[[1]]
    s = structure(names=as.character(rownames(attributes)[which(!is.na(attributes[,i]))]),attributes[which(!is.na(attributes[,i])),i])
    addAttributes(df,.attrs=s,append=TRUE)
   }
  }
 
  return(modelInternalDocument) 
}


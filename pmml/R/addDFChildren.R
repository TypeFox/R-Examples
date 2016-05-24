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

addDFChildren <- function(xmlmodel=NULL, field=NULL, intervals=NULL, values=NULL, namespace="4_2", ...) {

#to avoid malloc error
 flush.console()
 namespace <- .getNamespace(namespace)

# we expect the input to be always a XML Node
# If a file is to be used, the function 'fileToXMLNode' should be used to convert the file to a XMLNode object

  modelstring <- toString.XMLNode(xmlmodel)
  modelInternalDocument <- xmlTreeParse(modelstring,asText=TRUE,useInternalNodes=TRUE)
# DataDictionary and MiningSchema to be modified later
  dd <- getNodeSet(modelInternalDocument,"/p:PMML/p:DataDictionary",c(p=namespace))[[1]]
  formulaFields <- getNodeSet(modelInternalDocument,"/p:PMML/p:DataDictionary/p:DataField/@name",c(p=namespace))
  formulaString <- lapply(formulaFields,FUN=toString)

# convert field information from vector to data frame format
  fieldName <- field

# or use attributes <- read.table(filename,stringsAsFactors=FALSE)

  if(!is.element(fieldName,formulaFields))
  {
    print("WARNING: the following field additions will be ignored")
    stop(cat(fieldName,"not in model\n"))
  }

  df <- getNodeSet(dd,paste0("p:DataField[@name='",fieldName,"']"),c(p=namespace))[[1]]
  addChildren(df,kids=intervals,append=TRUE,supressNamespaceWarning=TRUE)
  addChildren(df,kids=values,append=TRUE,supressNamespaceWarning=TRUE)
  return(modelInternalDocument)
}


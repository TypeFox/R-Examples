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

#' Add Output nodes to a PMML object.
#' 
#' @param xmlmodel The PMML model to which the OutputField elements are to be added
#' @param outputNodes The Output nodes to be added. These may be created using the 
#'		'makeOutputNodes' helper function 
#' @param at Given an Output element, the 1 based index after which the given Output 
#'		child element should be inserted at 
#' @param xformText Post-processing information to be included in the OutputField element.
#'		This expression will be processed by the functionToPMML function
#' @param nodeName The name of the element to be added
#' @param attributes The attributes to be added 
#' @param whichOutput The name of the OutputField element
#' @param namespace The namespace of the PMML model 
#' 
#' @details
#'   This function is meant to add any post-processing information to an existing model via
#' the OutputField element. One can also use this to tell the PMML model to output other values
#' not automatically added to the model output.
#'   The first method is to use the 'makeOutputNodes' helper function to make a list of output 
#' elements to be added. 'whichOutput' lets the function know which of the Output elements we want to 
#' work with; there may be more than one in a multiple model file. One can then add those elements there,
#' at the desired index given by the 'at' parameter; the elements are inserted after the OutputField 
#' element at the 'at' index.
#'   This function can also be used with the 'nodeName' and 'attributes' to add the list of attributes to 
#' an OutputField element with name 'nodeName'
#' element using the 'xmlmodel', 'outputNodes' and 'at' parameters. 
#'   Finally, one can use this to add the transformation expression given by the 'xformText' parameter
#' to the node with name 'nodeName'. The string given via 'xformText' is converted to an XML expression similarly
#' to the functionToPMML function. 
#' 
#' @return Output node with the OutputField elements inserted.
#' 
#' @author Tridivesh Jena
#' 
#' @examples
#' # Load the standard iris dataset
#' data(iris)
#'
#' # Create a linear model and convert it to PMML
#' mod <- lm(Sepal.Length~.,iris)
#' pmod <- pmml(mod)
#' 
#' # Create additional output nodes
#' onodes0<-makeOutputNodes(name=list("OutputField","OutputField"),
#'                          attributes=list(list(name="dbl",
#'                          optype="continuous"),NULL),
#'                          expression=list("ln(x)","ln(x/(1-x))"))
#' onodes2<-makeOutputNodes(name=list("OutputField","OutputField"),
#'                          attributes=list(list(name="F1",
#'                          dataType="double",optype="continuous"),
#'                          list(name="F2")))
#' 
#' # Create new pmml objects with the output nodes appended
#' addOutputField(xmlmodel=pmod, outputNodes=onodes2, at="End", 
#'                xformText=NULL, nodeName=NULL, attributes=NULL,
#'                whichOutput=1)
#' pmod2<-addOutputField(xmlmodel=pmod, outputNodes=onodes0, at="End", 
#'                        xformText=NULL, nodeName=NULL, 
#'                        attributes=NULL,whichOutput=1)
#' 
#' # Create nodes with attributes and transformations
#' addOutputField(xmlmodel=pmod2, outputNodes=onodes2,at=2)
#' addOutputField(xmlmodel=pmod2, xformText=list("exp(x) && !x"), 
#'                nodeName="Predicted_Sepal.Length")
#' 
#' att <- list(datype="dbl",optpe="dsc")
#' addOutputField(xmlmodel=pmod2, nodeName="Predicted_Sepal.Length", 
#'                attributes=att)


addOutputField <- function(xmlmodel=NULL, outputNodes=NULL, at="End", xformText=NULL, nodeName=NULL,
                           attributes=NULL, whichOutput=1, namespace="4_2") {

#to avoid malloc error
 flush.console()
 namespace <- .getNamespace(namespace)

# we expect the input to be always a XML Node
# If a file is to be used, the function 'fileToXMLNode' should be used to convert the file to a XMLNode object

  modelstring <- toString.XMLNode(xmlmodel)
  modelInternalDocument <- xmlTreeParse(modelstring,asText=TRUE,useInternalNodes=TRUE)

  if(!is.null(outputNodes)) {
    warning("OutputField given, childNode and attributes ignored.")
    outNode <- getNodeSet(modelInternalDocument,paste0("/p:PMML/descendant::p:Output[",whichOutput,"]"), c(p=namespace))
    if(length(outNode) == 0)
      warning(paste0("No OutputField node with the name ",nodeName," found")) 
    if(at == "End") {
      addChildren(outNode[[1]],kids=list(outputNodes),append=TRUE,supressNamespaceWarning=TRUE)
    } else {
    # inserted after 'at' OutputField element 
      addChildren(outNode[[1]],kids=list(outputNodes),at=at,append=TRUE,supressNamespaceWarning=TRUE)
    }
  } 
  if(!is.null(attributes))
  {
    warning("attributes given, outputNode and childNode will be ignored")
    outputNode <- getNodeSet(modelInternalDocument,paste0("/p:PMML/descendant::p:Output/p:OutputField[@name='",nodeName,"']"),
                            c(p=namespace))
    for(i in 1:length(attributes)) {
      if(length(outputNode) == 1) {
        addAttributes(outputNode[[1]],.attrs=attributes[i])
      } else {
        addAttributes(outputNode[[1]],.attrs=attributes[[i]])
      }
    }
  }
  if(!is.null(xformText))
  {
    warning("childNode will be added. outputNode and attributes will be ignored.") 
    outputNode <- getNodeSet(modelInternalDocument,paste0("/p:PMML/descendant::p:Output/p:OutputField[@name='",nodeName,"']"),
                            c(p=namespace))
    for(i in 1:length(nodeName))  
      addChildren(outputNode[[i]],kids=list(.pmmlU(xformText[[i]])),append=TRUE,supressNamespaceWarning=TRUE)
  }
 
  return(modelInternalDocument)
}



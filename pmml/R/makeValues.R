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

#' Create Values element, most likely to add to a DataDictionary element 
#' 
#' @param value The 'value' attribute of each 'Value' element to be created in order.      
#' @param displayValue The 'displayValue' attribute of each 'Value' element to be created in order.
#' @param property The 'property' attribute of each 'Value' element to be created in order. 
#' @param namespace The namespace of the PMML model 
#' 
#' @details
#' The 'makeValues' function is used the same way as the 'makeIntervals' function. If certain attributes for an
#'  element should not be included, they should be input in the list as NULL.
#' 
#' @return PMML Values elements.
#' 
#' @author Tridivesh Jena
#' 
#' @examples
#' # define 3 values, none with a 'displayValue' attribute and 1 value 
#' # defined as 'invalid'. The 2nd one is 'valid' by default.
#' mv <- makeValues(list(1.1,2.2,3.3),list(NULL,NULL,NULL),
#'                  list("valid",NULL,"invalid"))
#' @seealso \code{\link{makeIntervals}} to make Interval child elements, \code{\link{addDFChildren}}
#' to add these xml fragments to the DataDictionary PMML element.

makeValues <- function(value=NULL,displayValue=NULL,property=NULL,namespace="4_2")
{
 namespace <- .getNamespace(namespace)
# if(is.vector(value))
 if((length(value) != length(displayValue)) || (length(value) != length(property)) 
     || (length(displayValue) != length(property)))
  stop("all parameters must have same length.")

 #property <- mapply(function(b){if(is.na(b)) return(NULL) else return(b)},property)
 nodes <- vector(mode="list",length=length(value))
 for(i in 1:length(value))
 {
  nodes[[i]] <- newXMLNode("Value",attrs=c(value=value[[i]],displayValue=displayValue[[i]],property=property[[i]]))
 }
 return(nodes) 
}


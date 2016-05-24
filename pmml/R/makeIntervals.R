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

#' Create Interval elements, most likely to add to a DataDictionary element 
#' 
#' @param closure The 'closure' attribute of each 'Interval' element to be created in order.
#' @param leftMargin The 'leftMargin' attribute of each 'Interval' element to be created in order. 
#' @param rightMargin The 'rightMargin' attribute of each 'Interval' element to be created in order. 
#' @param namespace The namespace of the PMML model 
#' 
#' @details
#'   The 'Interval' element allows 3 attributes, all of which may be defined in the 'makeIntervals' 
#' function. The value of these attributes should be provided as a list. Thus the elements of the 
#' 'leftMargin' for example define the value of that attribute for each 'Interval' element in order.
#' 
#' @return PMML Intervals elements.
#' 
#' @author Tridivesh Jena
#' 
#' @examples
#' # make 3 Interval elements
#' # we define the 3 Intervals as ,1]  (1,2)  and [2, 
#' mi<-makeIntervals(list("openClosed","openOpen","closedOpen"),
#'                    list(NULL,1,2),list(1,2,NULL))
#'
#' @seealso \code{\link{makeValues}} to make Values child elements, \code{\link{addDFChildren}}
#' to add these xml fragments to the DataDictionary PMML element.

makeIntervals <- function(closure=NULL,leftMargin=NULL,rightMargin=NULL, namespace="4_2")
{
 namespace <- .getNamespace(namespace)

# if(is.vector(closure))
 if((length(closure) != length(leftMargin)) || (length(closure) != length(rightMargin)) 
     || (length(leftMargin) != length(rightMargin)))
  stop("all parameters must have same length.")

 #rightMargin <- mapply(function(b){if(is.na(b)) return(NULL) else return(b)},rightMargin)
 nodes <- vector(mode="list",length=length(closure))
 for(i in 1:length(closure))
 {
  nodes[[i]] <- newXMLNode("Interval",attrs=c(closure=closure[[i]],leftMargin=leftMargin[[i]],rightMargin=rightMargin[[i]])) 
 }
 return(nodes) 
}



# PMML (Predictive Model Markup Language) Transformations 
#
# Copyright (c) 2015 Zementis, Inc.
#
# This file is part of the pmmlTransformations package. 
#
# The pmmlTransformations package is free: you can redistribute it and/or 
# modify it under the terms of the GNU General Public License as published 
# by the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# The pmmlTransformations package is distributed in the hope that it will 
# be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please see the
# GNU General Public License for details (http://www.gnu.org/licenses/).
############################################################################
#
# Author: Dmitriy Bolotov
#
#---------------------------------------------------------------------------


#' Add a function transformation to a WrapData object.
#' 
#' @param boxdata wrapper object obtained by using the WrapData function on raw data
#' @param origFieldName string specifying name(s) of the original data field(s) being used in the transformation
#' @param newFieldName name of the new field created by the transformation
#' @param newFieldDataType data type of the new field created by the transformation
#' @param formulaText string expression specifying the transformation
#' @param mapMissingTo value to be given to the transformed variable if the value of any input variable is missing
#'
#' @details 
#' 
#' Calculate the expression provided 
#' in formulaText for every row in the \code{boxdata$data} 
#' data frame. The \code{formulaText} argument must represent 
#' a valid R expression, and any functions used in 
#' \code{formulaText} must be defined in the current 
#' environment.
#' 
#' The name of the new field is optional (a default name is provided), but an error 
#' will be thrown if attempting to create a field with a name that already exists in 
#' the WrapData object.
#' 
#'
#' @return R object containing the raw data, the transformed data and data statistics. 
#' The \code{data} data frame will contain a new \code{newFieldName} column, and 
#' \code{fieldData} will contain a new \code{newFieldName} row.
#' 
#' @author Dmitriy Bolotov
#' 
#' @seealso \code{\link{WrapData}}
#' 
#' @examples
#' # Load the standard iris dataset
#' data(iris)
#' 
#' # Wrap the data
#' irisBox <- WrapData(iris)
#' 
#' # Perform a transform on the Sepal.Length field: 
#' # the value is squared and then divided by 100
#' irisBox <- FunctionXform(irisBox,origFieldName="Sepal.Length",
#'                          newFieldName="Sepal.Length.Transformed",
#'                          formulaText="(Sepal.Length^2)/100")
#' 
#' # Combine two fields to create another new feature:                      
#' irisBox <- FunctionXform(irisBox,
#'                          origFieldName="Sepal.Width, Petal.Width",
#'                          newFieldName="Width.Sum",
#'                          formulaText="Sepal.Width + Sepal.Length")
#'                          
#' # Create linear model using the derived features
#' fit <- lm(Petal.Length ~ 
#'          Sepal.Length.Transformed + Width.Sum, data=irisBox$data)
#' 
#' # Create pmml from the fit
#' library(pmml)
#' fit_pmml <- pmml(fit, transform=irisBox)

FunctionXform <- function (boxdata,origFieldName,newFieldName="newField",
                           newFieldDataType="numeric",formulaText,mapMissingTo=NA) {

  boxdata$data$newFieldName <- NA
  
  ## This loop makes it possible to apply an if-else formula to the new data column.
  for (n in 1:length(boxdata$data$newFieldName)) {
    boxrow <- boxdata$data[n,]
    boxdata$data$newFieldName[n] <- eval(parse(text=formulaText),boxrow)
  }
  
  names(boxdata$data)[names(boxdata$data)=="newFieldName"] <- newFieldName

  #new column for formula; only create if doesn't already exist; this is unnecessary if functionXform is already added by WrapData()
  if (!("functionXform" %in% colnames(boxdata$fieldData))) {
    boxdata$fieldData$functionXform <- "NA"
  }
    
  #make new row with "NA" entries
  temprow <- matrix(c(rep.int("NA",length(boxdata$fieldData))),nrow=1,ncol=length(boxdata$fieldData))
  newrow <- data.frame(temprow)
  colnames(newrow) <- colnames(boxdata$fieldData)
  boxdata$fieldData <- rbind(boxdata$fieldData,newrow)
  
  #add data to new row
  row.names(boxdata$fieldData)[nrow(boxdata$fieldData)] <- newFieldName
  
  levels(boxdata$fieldData$type)[2] <- "derived" #must create factor level first
  boxdata$fieldData[newFieldName,"type"] <- "derived"
  boxdata$fieldData[newFieldName,"dataType"] <- newFieldDataType #this could be string
  
  
  boxdata$fieldData[newFieldName,"functionXform"] <- formulaText
  
  #if origFieldName contains multiple fields, these will be combined into one string 
  boxdata$fieldData[newFieldName,"origFieldName"] <- paste(origFieldName,collapse=",")
  
  if(!is.na(mapMissingTo))
  {
    boxdata$fieldData[newFieldName,"missingValue"] <- mapMissingTo
  }
  
  return(boxdata)
}

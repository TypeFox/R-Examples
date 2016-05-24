# PMML (Predictive Model Markup Language) Transformations 
#
# Copyright (c) 2015 Zementis, Inc.
#
# This file is part of the pmmlTransformations package 
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
# Author: Tridivesh Jena
#
#---------------------------------------------------------------------------

RenameVar <- function(boxdata,xformInfo=NA,...)
{

	i <- NULL
        j <- NULL
	colnm <- NULL 
 
	boxData <- Initialize(boxdata)
  
	if(is.na(xformInfo))
	{
	  warning("No field name to rename found")
	  return(boxdata)
	} else
	{
	  # for each argument given
	  coln <- as.character(xformInfo)
	  # split to find initial and final names 	
	  if(grepl("[^-]->",coln))
	  {
	    st <- strsplit(coln,"->")
	  } else
	  {
	   st <- strsplit(coln,"-->")
	  }
          if(!is.na(st[[1]][2]))
          {
            derivedFieldName <- st[[1]][2]
          }
	  colnm <- st[[1]][1]
	  if(grepl("column",colnm,ignore.case=TRUE))
	  {
	    colnm <- gsub("column","",colnm,ignore.case=TRUE)
	  }
	  if(grepl("^[-,_]",colnm))
	  {
	    colnm <- gsub("^[-,_]*","",colnm)
	  }

          if(is.na(st[[1]][2]))
          {
            derivedFieldName <- paste("derived_",row.names(boxData$fieldData)[coln2],sep="")
          }

          # if column number, find the appropriate field 
	  if(suppressWarnings(!is.na(as.numeric(colnm))))
	  {
	    coln2 <- as.numeric(colnm)
	    dataType <- boxData$fieldData[names(boxData$data)[coln2],"dataType"]
	    if(dataType == "numeric")
	    {
	      row.names(boxData$fieldData)[coln2] <- derivedFieldName
	      names(boxData$data)[coln2] <- derivedFieldName
#new
	      if(!is.null(boxData$matrixData))
	      {
               names(boxData$matrixData)[coln2] <- derivedFieldName
	      }
            }
         } else
	 {
	   i <- which(names(boxData$data) == colnm)
	   if(is.null(i))
	   {
	     j <- which(names(boxData$data) == colnm)
	   }

	   if(is.null(i) && is.null(j))
	   {
	     stop("field name not found.")
	   }
           if(is.null(j))
           {
            row.names(boxData$fieldData)[i] <- derivedFieldName
            names(boxData$data)[i] <- derivedFieldName
#new
            if(!is.null(boxData$matrixData))
            {
             names(boxData$matrixData)[i] <- derivedFieldName
	    }
	   } else
	   {
            row.names(boxData$fieldData)[j] <- derivedFieldName
            names(boxData$data)[j] <- derivedFieldName
#new
            if(!is.null(boxData$matrixData))
            {
	     names(boxData$matrixData)[j] <- derivedFieldName
	    }
	   }
	 }
       }
	
       return(boxData)
}

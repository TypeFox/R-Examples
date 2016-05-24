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


WrapData <-function(indata,useMatrix=FALSE)
{
	dataBox <- NULL
	fieldNames <- NULL
	nrows <- NULL
	ncols <- NULL
	type <- NULL
	dataType <- NULL
	origFieldName <- NULL

	if(useMatrix)
 	{
	 dataBox$matrixData <- as.matrix(indata)
	} else
	{
	  dataBox$matrixData <- NULL
	}
        indatafrm <- data.frame(indata)
	dataBox$data <- indatafrm
	dataBox$nrows <- nrow(indatafrm)
	dataBox$ncols <- ncol(indatafrm)

#new
        if(is.matrix(indatafrm))
        {
          if(!is.numeric(indatafrm))
          {
            stop("Non-numeric matrices not yet supported for transformations")
          }
        }


#new 
        fieldNames <- names(indatafrm)

	for(i in 1:dataBox$ncols)
	{
		origFieldName <- NA
		type[i] <- "original"

		if(is.numeric(indata[,i]))
		{
		  dataType[i] <- "numeric"
		} else
		{
		  dataType[i] <- "factor"
		}

	}

    # mark all original field names by type=original and origFieldName=NA 
#     df<-data.frame(type,dataType,origFieldName,row.names=fieldNames)

    #add rest of fields with NA:
    sampleMin <- NA
    sampleMax <- NA
    xformedMin <- NA
    xformedMax <- NA
    centers <- NA
    scales <- NA
    fieldsMap <- NA
    transform <- NA
    default <- NA
    missingValue <- NA
    functionXform <- NA


    df <- data.frame(type,dataType,origFieldName,sampleMin,sampleMax,xformedMin,xformedMax,
                     centers,scales,fieldsMap,transform,default,missingValue,functionXform,row.names=fieldNames)

    dataBox$fieldData <- df 

	return(dataBox)
}

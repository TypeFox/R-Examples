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
# Author: Tridivesh Jena
#
#---------------------------------------------------------------------------

Initialize <- function(inbox)
{
        xformedMin <- NULL
        xformedMax <- NULL
        sampleMin <- NULL
        sampleMax <- NULL
        centers <- NULL
        scales <- NULL
        fieldsMap <- NULL
        transform <- NULL
      	default <- NULL
      	missingValue <- NULL
        functionXform <- NULL

#         print('start of Initialize')
        
        if(is.null(inbox$fieldData[1,"xformedMax"]))
        {
                for(i in 1:nrow(inbox$fieldData))
                {
                        xformedMin <- c(xformedMin,NA)
                        xformedMax <- c(xformedMax,NA)
                        sampleMin <- c(sampleMin,NA)
                        sampleMax <- c(sampleMax,NA)
                }
                inbox$fieldData <- cbind(inbox$fieldData,sampleMin)
                inbox$fieldData <- cbind(inbox$fieldData,sampleMax)
                inbox$fieldData <- cbind(inbox$fieldData,xformedMin)
                inbox$fieldData <- cbind(inbox$fieldData,xformedMax)
        }
        


        if(is.null(inbox$fieldData[1,"centers"]))
        {
                for(i in 1:nrow(inbox$fieldData))
                {
                        centers <- c(centers,NA)
                        scales <- c(scales,NA)
                }
                inbox$fieldData <- cbind(inbox$fieldData,centers)
                inbox$fieldData <- cbind(inbox$fieldData,scales)
        }

        if(is.null(inbox$fieldData[1,"fieldsMap"]))
        {
                for(i in 1:nrow(inbox$fieldData))
                {
                        fieldsMap <- c(fieldsMap,NA)
                }
                inbox$fieldData <- cbind(inbox$fieldData,fieldsMap)
        }

        if(is.null(inbox$fieldData[1,"transform"]))
        {
                for(i in 1:nrow(inbox$fieldData))
                {
                        transform <- c(transform,NA)
                }
                inbox$fieldData <- cbind(inbox$fieldData,transform)
        }

        if(is.null(inbox$fieldData[1,"default"]))
        {
                for(i in 1:nrow(inbox$fieldData))
                {
                        default <- c(default,NA)
                        missingValue <- c(missingValue,NA)
                }
                inbox$fieldData <- cbind(inbox$fieldData,default)
                inbox$fieldData <- cbind(inbox$fieldData,missingValue)
        }
  
        if(is.null(inbox$fieldData[1,"functionXform"]))
        {
                for(i in 1:nrow(inbox$fieldData))
                {
                      functionXform <- c(functionXform,NA)
                }
                inbox$fieldData <- cbind(inbox$fieldData,functionXform)
        }
        
        return(inbox)
}

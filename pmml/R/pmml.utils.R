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

.removeAsFactor <- function(fieldName)
{
    if(length(grep("as\\.factor\\(",fieldName)) == 1)
    {
        fieldName <- gsub("as.factor\\((\\w*)\\)","\\1", fieldName, perl=TRUE)
    } 
    return (fieldName)
}

.getNamespace <- function(x)
{
  return(paste("http://www.dmg.org/PMML-",x,sep=""))
}


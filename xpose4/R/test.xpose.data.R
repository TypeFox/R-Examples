# Xpose 4
# An R-based population pharmacokinetic/
# pharmacodynamic model building aid for NONMEM.
# Copyright (C) 1998-2004 E. Niclas Jonsson and Mats Karlsson.
# Copyright (C) 2005-2008 Andrew C. Hooker, Justin J. Wilkins, 
# Mats O. Karlsson and E. Niclas Jonsson.
# Copyright (C) 2009-2010 Andrew C. Hooker, Mats O. Karlsson and 
# E. Niclas Jonsson.

# This file is a part of Xpose 4.
# Xpose 4 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with this program.  A copy can be cound in the R installation
# directory under \share\licenses. If not, see http://www.gnu.org/licenses/.

test.xpose.data <- function(object) {

#  if(is.null(object@Data)) return("The object contains no Data")

  if(!is.null(SData(object))) {
    sdata <- SData(object)
    data  <- Data(object)

    if((ncol(data)+1) != ncol(sdata)) {
      return("The Data and the SData do not have the same number of columns!")
    }
    
    ## Check columns
    nams <- names(data)
    snams<- names(sdata)

    ## All columns in Data should be in SData
    for(n in nams) {
      if(any(n==snams)) {
      } else {
        return(paste(n," does not seem to be present in SData!"))
      }
    }

    ## All columns (except iter) should be in Data
    for(n in snams) {
      if(n == "iter") next
      if(any(n==nams)) {
      } else {
        return(paste(n," does not seem to be present in Data!"))
      }
    }

    ## We also need to check that the class definitions are the same
    for(n in nams) {
      if(class(data[,n]) != class(sdata[,n]))
        return(paste(n," does not have the same class definition in Data and SData!"))
    }

    ## Check that SData is an even multiple of Data
    dnro <- dim(data)[1]
    snro <- dim(sdata)[1]
    #  cat(snro)
    if(regexpr("\\.",as.character(snro/snro)) !=-1) {
      return("The lengths of Data and SData do not match!")
    }
    
  }
  
  return(TRUE)
}

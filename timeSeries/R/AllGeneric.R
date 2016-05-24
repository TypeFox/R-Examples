#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# GENERIC:                   DESCRIPTION
#  returns                    Computes returns
#  rowCumsums                 Computes row cumulated sums
#  series                     Extracts series data
#  series<-                   Assigns series data  
#  coredata                   Extracts series data
#  index                      deprecated
#  index <-                   deprecated
#  outlier                    Returns outliers
#  timeSeries                 Returns timeSeries
#  colCumsums                 Computes column cumulated sums
#  colCummaxs                 Computes column cumulated maxima
#  colCummins                 Computes column cumulated minima
#  colCumprods                Computes column cumulated products
#  colCumreturns              Computes column cumulated returns
################################################################################


setGeneric("returns", 
    function(x, ...)
    standardGeneric("returns"), package = "timeSeries")

    
setGeneric("rowCumsums", 
    function(x, na.rm = FALSE, ...)
    standardGeneric("rowCumsums"), package = "timeSeries")

    
setGeneric("series", 
    function(x) 
    standardGeneric("series"), package = "timeSeries")
    
    
setGeneric("series<-", 
    function(x, value)
    standardGeneric("series<-"), package = "timeSeries")


setGeneric("coredata", 
    function(x) 
    standardGeneric("coredata"), package = "timeSeries")
    
    
setGeneric("coredata<-", 
    function(x, value)
    standardGeneric("coredata<-"), package = "timeSeries")


    
## setGeneric("index", function(x, ...)
##      standardGeneric("index"), package = "timeSeries")


## setGeneric("index<-", function(x, value)
##      standardGeneric("index<-"), package = "timeSeries")


setGeneric("outlier",  
    function(x, sd = 5, complement = TRUE, ...)
    standardGeneric("outlier"))

           
setGeneric("timeSeries",
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    standardGeneric("timeSeries"))
    

setGeneric("colCumsums", 
    function(x, na.rm = FALSE, ...) 
    standardGeneric("colCumsums"))
    
    
setGeneric("colCummaxs", 
    function(x, na.rm = FALSE, ...) 
    standardGeneric("colCummaxs"))
    
    
setGeneric("colCummins", 
    function(x, na.rm = FALSE, ...) 
    standardGeneric("colCummins"))
    
    
setGeneric("colCumprods", 
    function(x, na.rm = FALSE, ...) 
    standardGeneric("colCumprods"))
    
    
setGeneric("colCumreturns",
    function(x, method = c("geometric", "simple"), na.rm = FALSE, ...)
    standardGeneric("colCumreturns"))


################################################################################


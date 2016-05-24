require(foreign)
#' @title Read SaTScan output files
#' 
#' @description
#' Reads a SaTScan output .dbf file.
#' 
#' @details
#' This is expected to be a purely internal function.  
#' It's called by \code{satscan()} with the location and file nmae provided 
#' to that function.  Since it's nothing more than \code{foreign::read.dbf()}, it's 
#' probably nor necessary to even have it as a function.
#' 
#' @param location A directory location, including the trailing "/"
#' @param file A file name, without the extension.
#' 
#' @return A data frame.
read.col = function(location,file) foreign::read.dbf(paste0(location,file,".col.dbf"))


#' @title Read SaTScan output files
#' 
#' @description
#' Reads a SaTScan output .dbf file.
#' 
#' @details
#' This is expected to be a purely internal function.  
#' It's called by \code{satscan()} with the location and file nmae provided 
#' to that function.  Since it's nothing more than \code{foreign::read.dbf()}, it's 
#' probably nor necessary to even have it as a function.
#' 
#' @param location A directory location, including the trailing "/"
#' @param file A file name, without the extension.
#' 
#' @return A data frame.
read.rr = function(location,file) foreign::read.dbf(paste0(location,file,".rr.dbf"))

#' @title Read SaTScan output files
#' 
#' @description
#' Reads a SaTScan output .dbf file.
#' 
#' @details
#' This is expected to be a purely internal function.  
#' It's called by \code{satscan()} with the location and file nmae provided 
#' to that function.  Since it's nothing more than \code{foreign::read.dbf()}, it's 
#' probably nor necessary to even have it as a function.
#' 
#' @param location A directory location, including the trailing "/"
#' @param file A file name, without the extension.
#' 
#' @return A data frame.
read.gis = function(location,file) foreign::read.dbf(paste0(location,file,".gis.dbf"))

#' @title Read SaTScan output files
#' 
#' @description
#' Reads a SaTScan output .dbf file.
#' 
#' @details
#' This is expected to be a purely internal function.  
#' It's called by \code{satscan()} with the location and file nmae provided 
#' to that function.  Since it's nothing more than \code{foreign::read.dbf()}, it's 
#' probably nor necessary to even have it as a function.
#' 
#' @param location A directory location, including the trailing "/"
#' @param file A file name, without the extension.
#' 
#' @return A data frame.
read.llr = function(location,file) foreign::read.dbf(paste0(location,file,".llr.dbf"))

#' @title Read SaTScan output files
#' 
#' @description
#' Reads a SaTScan output .dbf file.
#' 
#' @details
#' This is expected to be a purely internal function.  
#' It's called by \code{satscan()} with the location and file nmae provided 
#' to that function.  Since it's nothing more than \code{foreign::read.dbf()}, it's 
#' probably nor necessary to even have it as a function.
#' 
#' @param location A directory location, including the trailing "/"
#' @param file A file name, without the extension.
#' 
#' @return A data frame.
read.sci = function(location,file) foreign::read.dbf(paste0(location,file,".sci.dbf"))

#' @title Read SaTScan output files
#' 
#' @description
#' Reads a SaTScan output .dbf file.
#' 
#' @details
#' This is expected to be a purely internal function.  
#' It's called by \code{satscan()} with the location and file nmae provided 
#' to that function.  Since it's nothing more than \code{readLines()}, it's 
#' probably nor necessary to even have it as a function.
#' 
#' @param location A directory location, including the trailing "/"
#' @param file A file name, without the extension.
#' 
#' @return A data frame.
read.satscanmain = function(location,file) suppressWarnings(
  readLines(paste0(location,file,".txt")))

# need to get these to read dates correctly.  Currently reading as factors
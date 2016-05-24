#############################################################################
#
#   This file is part of the R package "RSNNS".
#
#   Author: Christoph Bergmeir
#   Supervisor: José M. Benítez
#   Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Library General Public
#   License as published by the Free Software Foundation; either
#   version 2 of the License, or (at your option) any later version.
# 
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Library General Public License for more details.
# 
#   You should have received a copy of the GNU Library General Public License
#   along with this library; see the file COPYING.LIB.  If not, write to
#   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#   Boston, MA 02110-1301, USA.
#
#############################################################################


#' This function contains a rudimentary parser for SNNS .res files. It
#' is completely implemented in R and doesn't make use of SNNS functionality.
#' 
#' @title Rudimentary parser for res files.
#' @param filename the name of the .res file 
#' @return a matrix containing the predicted values that were found in the .res file
#' @export
readResFile <- function(filename)  {
  
  allData <- scan(filename, what="list", multi.line=TRUE)
  
  table <- NULL
  mode <- "omitHeader"
  lineName <- NULL
  line <- NULL
  for(i in 1:length(allData)) {
    
    if(mode == "omitHeader")
    {
      firstSymbol <- substr(allData[i],1,1)
      if(!is.na(firstSymbol) && firstSymbol == "#")  {
        lineName <- allData[i]
        mode <- "parseData"
      }
    } else if(mode == "parseData")  {
      firstSymbol <- substr(allData[i],1,1)
      
      if(!is.na(firstSymbol) && firstSymbol == "#")  {
        table <- rbind(table, line)
        rownames(table)[nrow(table)] <- lineName
        
        lineName <- allData[i]
        line <- NULL
      } else  {
        line <- c(line, as.numeric(allData[i]))
      }    
    }
  }
  table <- rbind(table, line)
  rownames(table)[nrow(table)] <- lineName
  
  #table
  prediction <- table[,ncol(table)]
  return(prediction)
  
}


#' This function generates an \link{SnnsR-class} object, loads the given data there 
#' as a pattern set and then uses the functionality of SNNS to save the data as a .pat file. 
#' 
#' @title Save data to a pat file
#' @param inputs a matrix with input values 
#' @param targets a matrix with target values
#' @param filename the name of the .pat file 
#' @export
savePatFile <- function(inputs, targets, filename)  {

  snnsObject <- SnnsRObjectFactory()

  snnsObject$createNet(unitsPerLayer=c(ncol(as.matrix(inputs)),ncol(as.matrix(targets))), fullyConnectedFeedForward = TRUE)
  
  patSet <- snnsObject$createPatSet(inputs, targets)
  snnsObject$saveNewPatterns(filename, patSet$set_no)
  
  rm(snnsObject)
  
}


#' This function generates an \link{SnnsR-class} object, loads the given .pat file 
#' there as a pattern set and then extracts the patterns to a matrix, 
#' using \link{SnnsRObject$extractPatterns}. 
#' 
#' @title Load data from a pat file
#' @param filename the name of the .pat file
#' @return a matrix containing the data loaded from the .pat file. 
#' @export
readPatFile <- function(filename)  {
    
  snnsObject <- SnnsRObjectFactory()
  
  snnsObject$loadNewPatterns(filename)
  
  patterns <- snnsObject$extractPatterns()
  
  rm(snnsObject)
  
  return(patterns)
  
}


#' This function extracts all columns from a matrix whose column names begin with "in".
#' The example data of this package follows this naming convention. 
#' 
#' @title Get the columns that are inputs
#' @param patterns matrix or data.frame containing the patterns 
#' @export
inputColumns <- function(patterns)  {
  
  res <- which(substr(colnames(patterns),1,2) == "in")
  return(res)

}


#' This function extracts all columns from a matrix whose column names begin with "out".
#' The example data of this package follows this naming convention. 
#' 
#' @title Get the columns that are targets
#' @param patterns matrix or data.frame containing the patterns 
#' @export
outputColumns <- function(patterns)  {
  
  res <- which(substr(colnames(patterns),1,3) == "out")
  return(res)

}


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


#' SnnsR low-level function that extracts all patterns of the current pattern set and
#' returns them as a matrix. Columns are named with the prefix "in" or "out", respectively.
#'  
#' @title Extract the current pattern set to a matrix
#' @return a matrix containing the patterns of the currently loaded patern set.
#' @rdname SnnsRObject-extractPatterns
#' @name SnnsRObject$extractPatterns
#' @usage \S4method{extractPatterns}{SnnsR}()
#' @aliases extractPatterns,SnnsR-method SnnsR__extractPatterns
SnnsR__extractPatterns <- function(snnsObject)  {
  
  noPatterns <- snnsObject$getNoOfPatterns()
  
  inputs <- NULL
  outputs <- NULL
  
  for(i in 1:noPatterns)  {
    #INPUT: 1
    input <- snnsObject$getSubPatData(i-1, 0, 1)
    #OUTPUT: 2
    output <- snnsObject$getSubPatData(i-1, 0, 2)
    
    inputs <- rbind(inputs, input)
    outputs <- rbind(outputs, output)
  }
  
  colnames(inputs) <- paste("in", 1:ncol(inputs), sep="") 
  
  table <- inputs
  
  if(ncol(outputs) != 0)  {
    colnames(outputs) <- paste("out", 1:ncol(outputs), sep="")
    table <- cbind(table, outputs)
  } 
    
  
  #inputs
  #outputs
  
  rownames(table) <- paste("pattern", 1:nrow(table), sep="")
  
  return(table)
}
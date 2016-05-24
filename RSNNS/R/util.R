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


##Find out if a string has a certain suffix
endsWith <- function(myString, mySubString) {
  
  l1 <- nchar(myString)
  l2 <- nchar(mySubString)
  if(l1 < l2) return(FALSE)
  s <- substr(myString, l1-l2+1, l1)
  if(s == mySubString) return(TRUE)
  return(FALSE)
  
}

##Find out if a string has a certain prefix
beginsWith <- function(myString, mySubString) {
  
  l1 <- nchar(myString)
  l2 <- nchar(mySubString)
  if(l1 < l2) return(FALSE)
  s <- substr(myString, 1, l2)
  if(s == mySubString) return(TRUE)
  return(FALSE)
  
}

##Rotate a square matrix 90 degrees clockwise.
rot90 <- function(a) {
  n <- dim(a)[1]
  t(a[n:1, ])
}

#' Organize network activation as 2d map.
#'
#' The input to this function is a vector containing in each row an activation
#' pattern/output of a neural network. This function reorganizes the vector to
#' a matrix. Normally, only the number of rows \code{nrow} will be used.
#' 
#' @title Convert a vector to an activation map
#' @param v the vector containing the activation pattern
#' @param nrow number of rows the resulting matrices will have
#' @param ncol number of columns the resulting matrices will have
#' @return a matrix containing the 2d reorganized input
#' @export
#' @seealso \link{matrixToActMapList} \link{plotActMap}
vectorToActMap <- function(v, nrow=0, ncol=0) {
  
  if(nrow==0) nrow <- ncol(v) / ncol  
  return(matrix(v, nrow=nrow, byrow=TRUE))
  
}
  
#' Organize a matrix containing 1d vectors of network activations as 2d maps.
#'
#' The input to this function is a matrix containing in each row an activation
#' pattern/output of a neural network. This function uses \link{vectorToActMap} to 
#' reorganize the matrix to a list of matrices, whereby each row of the input matrix 
#' is converted to a matrix in the output list.
#' 
#' @title Convert matrix of activations to activation map list
#' @param m the matrix containing one activation pattern in every row
#' @param nrow number of rows the resulting matrices will have
#' @param ncol number of columns the resulting matrices will have
#' @return a list containing the activation map matrices
#' @export
#' @seealso \link{vectorToActMap} \link{plotActMap}
matrixToActMapList <- function(m, nrow=0, ncol=0) {
  
  actMapList <- apply(m, 1, function(x) { return(list(vectorToActMap(x,nrow,ncol)))})
  actMapList <- lapply(actMapList, function(x) {x[[1]]})
  actMapList
} 

#' Plot an activation map as a heatmap.
#'
#' @title Plot activation map
#' @param x the input data matrix 
#' @param ... parameters passed to \code{image}
#' @export
#' @seealso \link{vectorToActMap} \link{matrixToActMapList}
plotActMap <- function(x, ...) {
  image(rot90(x),...)
}

#' Get the function table of available SNNS functions.
#'
#' @title Get SnnsR function table
#' @return a data.frame with columns:
#' \item{name}{name of the function}
#' \item{type}{the type of the function (learning, init, update,...)}
#' \item{#inParams}{the number of input parameters of the function}
#' \item{#outParams}{the number of output parameters of the function}
#' @export
getSnnsRFunctionTable <- function() {

  snnsObject <- SnnsRObjectFactory()
  
  noFunc <- snnsObject$getNoOfFunctions()
  allFuncs <- data.frame()
  
  for(i in 1:noFunc) {
    fi <- snnsObject$getFuncInfo(i)
    fiInfo <- snnsObject$getFuncParamInfo(fi[[1]], fi[[2]])
    allFuncs <- rbind(allFuncs, cbind(fi$func_name, fi$func_type, fiInfo$no_of_input_params, fiInfo$no_of_input_params))
  }
  
  names(allFuncs) <- c("name", "type", "#inParams", "#outParams")
  
  rm(snnsObject)
  
  allFuncs
}

#' Get a define of the SNNS kernel from a defines-list.
#' All defines-lists present can be shown with \code{RSNNS:::SnnsDefines}. 
#' 
#' @title Get a define of the SNNS kernel
#' @param defList the defines-list from which to get the define from
#' @param defValue the value in the list
#' @return a string with the name of the define
#' @export           
#' @seealso \code{\link{resolveSnnsRDefine}}
#' @examples
#' getSnnsRDefine("topologicalUnitTypes",3)
#' getSnnsRDefine("errorCodes",-50)
getSnnsRDefine <- function(defList, defValue)  {
  defRow <- which(SnnsDefines[[defList]][,2] == toString(defValue))
  return(SnnsDefines[[defList]][defRow,1])
}


#' Resolve a define of the SNNS kernel using a defines-list.
#' All defines-lists present can be shown with \code{RSNNS:::SnnsDefines}.
#' 
#' @title Resolve a define of the SNNS kernel
#' @param defList the defines-list from which to resolve the define from
#' @param def the name of the define
#' @return the value of the define
#' @export
#' @seealso \code{\link{getSnnsRDefine}}           
#' @examples
#' resolveSnnsRDefine("topologicalUnitTypes","UNIT_HIDDEN")
resolveSnnsRDefine <- function(defList, def)  {
  defRow <- which(SnnsDefines[[defList]][,1] == toString(def))
  return(as.numeric(SnnsDefines[[defList]][defRow,2]))  
}

#SnnsDefines_showWarningFromSnnsError <- function(func, err) {
#  warning(paste("An error occured in ", func,": ", SnnsDefines_getDefine(SnnsDefines_errorCodes, err),sep=""))
#}

# Set the seed value used in all \code{SnnsR} objects.
# The seed value is used in the constructor of 
# every \code{SnnsCLib} object to set the seed of rand().
#' DEPRECATED, now just calls R's set.seed(), that should be used instead.
#'
#' @title DEPRECATED, Set the SnnsR seed value
#' @param seed the seed to use. If 0, a seed based on the system time is generated.
#' @export
setSnnsRSeedValue <- function(seed) {
  #.Call("setCurrentSeedVal", seed, package="RSNNS")  
  warning("Function setSnnsRSeedValue is deprecated. Now the R RNG is used, so use set.seed() instead")
  set.seed(seed)
}

getKrioTitle <- function(title_num) {
  .Call("getKrioTitle", title_num, package="RSNNS")  
}

is.nil <- function(ptr) {
  
  if (class(ptr)!="externalptr") stop("argument given is not a pointer")
  
  res <- .Call("isnil", ptr)
  
  res
}

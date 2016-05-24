#' @title Extract predictions from an object of class \code{mbl}
#' @description Extract predictions from an object of class \code{mbl}
#' @usage 
#' getPredictions(object)
#' @param object an object of class \code{mbl} as returned by \code{mbl}
#' @return a \code{data.frame} of predicted values according to either \code{k} or \code{k.dist}
#' @author Leonardo Ramirez-Lopez and Antoine Stevens
#' @seealso \code{\link{mbl}}
#' @examples
#' \dontrun{
#' require(prospectr)
#' 
#' data(NIRsoil)
#' 
#' Xu <- NIRsoil$spc[!as.logical(NIRsoil$train),]
#' Yu <- NIRsoil$CEC[!as.logical(NIRsoil$train)]
#' Yr <- NIRsoil$CEC[as.logical(NIRsoil$train)]
#' Xr <- NIRsoil$spc[as.logical(NIRsoil$train),]
#' 
#' Xu <- Xu[!is.na(Yu),]
#' Yu <- Yu[!is.na(Yu)]
#' 
#' Xr <- Xr[!is.na(Yr),]
#' Yr <- Yr[!is.na(Yr)] 
#'
#' ctrl <- mblControl(sm = "pls", 
#'                    pcSelection = list("opc", 40), 
#'                    valMethod = c("NNv"), 
#'                    scaled = TRUE, center = TRUE)
#'
#' ex1 <- mbl(Yr = Yr, Xr = Xr, Yu = NULL, Xu = Xu,
#'            mblCtrl = ctrl,
#'            distUsage = "predictors", 
#'            k = seq(30, 150, 15), 
#'            method = "wapls1",
#'            plsC = c(7, 20))
#'
#' getPredictions(ex1)
#' }
#' @export
######################################################################
# resemble
# Copyrigth (C) 2014 Leonardo Ramirez-Lopez and Antoine Stevens
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
######################################################################
getPredictions <- function(object){
  if(is.na(match("mbl", class(object))))
    stop("the object is not of class 'mbl'")
  
  ext.pred <- function(x,...){
    prediction <- x$pred
    return(prediction)
  }
  predictions <- as.data.frame(sapply(object$results, ext.pred))
  return(predictions)
}

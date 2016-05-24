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


#' This SnnsR low-level function masks the SNNS kernel function of the same name
#' to allow for both giving the initialization function
#' directly in the call or to use the one that is currently set.
#'  
#' @title Initialize the network
#' @param parameterInArray the parameters of the initialization function
#' @param initFunc the name of the initialization function
#' @rdname SnnsRObject-initializeNet
#' @name SnnsRObject$initializeNet
#' @usage \S4method{initializeNet}{SnnsR}(parameterInArray, initFunc)
#' @aliases initializeNet,SnnsR-method SnnsR__initializeNet
SnnsR__initializeNet <- function(snnsObject, parameterInArray, initFunc) {
  
  if(!missing(initFunc)) {
    .Call("SnnsCLib__setInitialisationFunc", snnsObject@variables$snnsCLibPointer, initFunc, package="RSNNS")
  }
  err <- .Call("SnnsCLib__initializeNet", snnsObject@variables$snnsCLibPointer, parameterInArray, package="RSNNS")
  err
}


#' This SnnsR low-level function masks the SNNS kernel function of the same name 
#' to allow both for giving the parameters directly or as a vector.
#' If the second parameter, \code{bias}, is missing, it is assumed 
#' that the first parameter should be interpreted as a vector containing
#' all parameters.
#' 
#' @title Set the unit defaults
#' @param act same as SNNS kernel function
#' @param bias idem
#' @param st idem
#' @param subnet_no idem
#' @param layer_no idem
#' @param act_func idem
#' @param out_func idem
#' @rdname SnnsRObject-setUnitDefaults
#' @name SnnsRObject$setUnitDefaults
#' @usage \S4method{setUnitDefaults}{SnnsR}(act, bias, st, subnet_no, layer_no, act_func, out_func)
#' @aliases setUnitDefaults,SnnsR-method SnnsR__setUnitDefaults
SnnsR__setUnitDefaults <- function(snnsObject, act, bias, st, subnet_no, layer_no, act_func, out_func) {
  
  #act <- c(1,0,1,0,1,"Act_Logistic","Out_Identity")
  
  if(missing(bias) && length(act)==7){

    bias <- as.numeric(act[2])
    st <- as.numeric(act[3])
    subnet_no <- as.numeric(act[4])
    layer_no <- as.numeric(act[5])
    act_func <- act[6]
    out_func <- act[7]
    act <- as.numeric(act[1])    
  }
    
  err <- .Call("SnnsCLib__setUnitDefaults", snnsObject@variables$snnsCLibPointer, act, bias, st, subnet_no, layer_no, 
      act_func, out_func, package="RSNNS")
  err
}
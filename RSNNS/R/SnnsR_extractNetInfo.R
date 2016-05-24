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


#' SnnsR low-level function to get all units in the net of a certain \code{ttype}. 
#' Possible \code{ttype} defined by SNNS are, among others:
#' "UNIT_OUTPUT", "UNIT_INPUT", and "UNIT_HIDDEN". For a full list, 
#' call \code{RSNNS:::SnnsDefines$topologicalUnitTypes}
#' As this is an SnnsR low-level function, you may want to have a look 
#' at \code{\link{SnnsR-class}} to find out how to properly use it.
#' 
#' @title Get all units in the net of a certain \code{ttype}.
#' @param ttype a string containing the \code{ttype}.
#' @return a vector with integer numbers identifying the units.
#' @rdname SnnsRObject-getAllUnitsTType
#' @name SnnsRObject$getAllUnitsTType
#' @usage \S4method{getAllUnitsTType}{SnnsR}(ttype)
#' @aliases getAllUnitsTType,SnnsR-method SnnsR__getAllUnitsTType
#' @seealso \link{SnnsRObject$getAllOutputUnits}, \link{SnnsRObject$getAllInputUnits}, \link{SnnsRObject$getAllHiddenUnits}
SnnsR__getAllUnitsTType <- function(snnsObject, ttype) {
  
  res <- NULL
  
  resolvedTType <- resolveSnnsRDefine("topologicalUnitTypes", ttype)
  
  nUnits <- snnsObject$getNoOfUnits()
  
  for(i in 1:nUnits)  {
    
    if(i==1)  unit <- snnsObject$getFirstUnit()
    else unit <- snnsObject$getNextUnit()
    
    type <- snnsObject$getUnitTType(unit)
    if(type == resolvedTType) res <- c(res, unit)
  }
  
  res
}


#' SnnsR low-level function to get all units from the net with the ttype "UNIT_OUTPUT".
#' This function calls \code{\link{SnnsRObject$getAllUnitsTType}} with the parameter "UNIT_OUTPUT".
#' 
#' @title Get all output units of the net.
#' @return a vector with integer numbers identifying the units.
#' @rdname SnnsRObject-getAllOutputUnits
#' @name SnnsRObject$getAllOutputUnits
#' @usage \S4method{getAllOutputUnits}{SnnsR}()
#' @aliases getAllOutputUnits,SnnsR-method SnnsR__getAllOutputUnits
#' @seealso \code{\link{SnnsRObject$getAllUnitsTType}}
SnnsR__getAllOutputUnits <- function(snnsObject) {
  c(snnsObject$getAllUnitsTType("UNIT_OUTPUT"), snnsObject$getAllUnitsTType("UNIT_SPECIAL_O"))  
}

#' SnnsR low-level function to get all units from the net with the ttype "UNIT_INPUT".
#' This function calls \code{\link{SnnsRObject$getAllUnitsTType}} with the parameter "UNIT_INPUT".
#' 
#' @title Get all input units of the net
#' @return a vector with integer numbers identifying the units.
#' @rdname SnnsRObject-getAllInputUnits
#' @name SnnsRObject$getAllInputUnits
#' @usage \S4method{getAllInputUnits}{SnnsR}()
#' @aliases getAllInputUnits,SnnsR-method SnnsR__getAllInputUnits
#' @seealso \link{SnnsRObject$getAllUnitsTType}
SnnsR__getAllInputUnits <- function(snnsObject) {
  c(snnsObject$getAllUnitsTType("UNIT_INPUT"), snnsObject$getAllUnitsTType("UNIT_SPECIAL_I"))  
}

#' SnnsR low-level function to get all units from the net with the ttype "UNIT_HIDDEN".
#' This function calls \code{\link{SnnsRObject$getAllUnitsTType}} with the parameter "UNIT_HIDDEN".
#' 
#' @title Get all hidden units of the net
#' @return a vector with integer numbers identifying the units.
#' @rdname SnnsRObject-getAllHiddenUnits
#' @name SnnsRObject$getAllHiddenUnits
#' @usage \S4method{getAllHiddenUnits}{SnnsR}()
#' @aliases getAllHiddenUnits,SnnsR-method SnnsR__getAllHiddenUnits
#' @seealso \link{SnnsRObject$getAllUnitsTType}
SnnsR__getAllHiddenUnits <- function(snnsObject) {
  c(snnsObject$getAllUnitsTType("UNIT_HIDDEN"), snnsObject$getAllUnitsTType("UNIT_SPECIAL_H")  )  
}


#' SnnsR low-level function to get the weight matrix between two sets of units.
#' 
#' @title Get the weight matrix between two sets of units
#' @param unitsSource a vector with numbers identifying the source units
#' @param unitsTarget a vector with numbers identifying the target units
#' @param setDimNames indicates, whether names of units are extracted and set as row/col names in the weight matrix
#' @return the weight matrix between the two sets of neurons 
#' @rdname SnnsRObject-getWeightMatrix
#' @name SnnsRObject$getWeightMatrix
#' @usage \S4method{getWeightMatrix}{SnnsR}(unitsSource, unitsTarget, setDimNames)
#' @aliases getWeightMatrix,SnnsR-method SnnsR__getWeightMatrix
#' @seealso \link{SnnsRObject$getAllUnitsTType}
SnnsR__getWeightMatrix <- function (snnsObject, unitsSource, unitsTarget, setDimNames=TRUE) {
  
  blank <- " "

  C <- matrix(nrow=length(unitsSource), ncol=length(unitsTarget))
  
  if(setDimNames) {
    
    unitsSourceNames <- NULL
    unitsTargetNames <- NULL
    
    for(i in 1:length(unitsSource)) {
      name <- snnsObject$getUnitName(unitsSource[i])
      if(is.null(name)) name <- blank
      
      unitsSourceNames <- c(unitsSourceNames, name)
    }
    
    for(j in 1:length(unitsTarget)) {
      name <- snnsObject$getUnitName(unitsTarget[j])
      if(is.null(name)) name <- blank
      
      unitsTargetNames <- c(unitsTargetNames, name)
      
    }

    rownames(C) <- unitsSourceNames    
    colnames(C) <- unitsTargetNames    
  }
  
  for(i in 1:length(unitsSource))
    for(j in 1:length(unitsTarget)) {
      res <- snnsObject$areConnectedWeight(unitsSource[i],unitsTarget[j])
      C[i,j] <- res$weight
    }
  C
}

#\link{getAllUnitsTType,SnnsR-method}

#' Set the activation function for all units of a certain ttype.
#' 
#' The function uses the function \code{\link{SnnsRObject$getAllUnitsTType}} to find all units of a certain
#' \code{ttype}, and sets the activation function of all these units to the given activation function.
#'  
#' @param ttype a string containing the \code{ttype}.
#' @param act_func the name of the activation function to set.
#' @rdname SnnsRObject-setTTypeUnitsActFunc
#' @name SnnsRObject$setTTypeUnitsActFunc
#' @usage \S4method{setTTypeUnitsActFunc}{SnnsR}(ttype, act_func)
#' @aliases setTTypeUnitsActFunc,SnnsR-method SnnsR__setTTypeUnitsActFunc
#' @seealso \code{\link{SnnsRObject$getAllUnitsTType}}
#' @examples
#' \dontrun{SnnsRObject$setTTypeUnitsActFunc("UNIT_HIDDEN", "Act_Logistic")}
SnnsR__setTTypeUnitsActFunc <- function(snnsObject, ttype, act_func) {
  
  units <- snnsObject$getAllUnitsTType(ttype)
  
  for(unit in units) {
    snnsObject$setUnitActFunc(unit, act_func)
  }

}

#' Get all units present in the net.
#' 
#' @return a vector with integer numbers identifying the units.
#' @rdname SnnsRObject-getAllUnits
#' @name SnnsRObject$getAllUnits
#' @usage \S4method{getAllUnits}{SnnsR}()
#' @aliases getAllUnits,SnnsR-method SnnsR__getAllUnits
SnnsR__getAllUnits <- function(snnsObject) {
  
  #res <- data.frame()
  res <- NULL
  
  nUnits <- snnsObject$getNoOfUnits()
  
  for(i in 1:nUnits)  {
    if(i==1)  unit <- snnsObject$getFirstUnit()
    else unit <- snnsObject$getNextUnit()
    
    #name <- snnsObject$getUnitName(unit)
    #res <- rbind(res, data.frame(name=name, unit=unit))
    res <- c(res, unit)
  }
  
  res
}

#' Get the complete weight matrix.
#' 
#' Get a weight matrix containing all weights of all neurons present in the net.
#' 
#' @param setDimNames indicates, whether names of units are extracted and set as row/col names in the weight matrix
#' @return the complete weight matrix
#' @rdname SnnsRObject-getCompleteWeightMatrix
#' @name SnnsRObject$getCompleteWeightMatrix
#' @usage \S4method{getCompleteWeightMatrix}{SnnsR}(setDimNames)
#' @aliases getCompleteWeightMatrix,SnnsR-method SnnsR__getCompleteWeightMatrix
SnnsR__getCompleteWeightMatrix <- function(snnsObject, setDimNames=TRUE) {
  
  allUnits <- snnsObject$getAllUnits()
  res <- snnsObject$getWeightMatrix(allUnits, allUnits, setDimNames=setDimNames)
  res
}


#' Find all units whose name begins with a given prefix.
#' 
#' @param prefix a prefix that the names of the units to find have.
#' @return a vector with integer numbers identifying the units.
#' @rdname SnnsRObject-getUnitsByName
#' @name SnnsRObject$getUnitsByName
#' @usage \S4method{getUnitsByName}{SnnsR}(prefix)
#' @aliases getUnitsByName,SnnsR-method SnnsR__getUnitsByName
SnnsR__getUnitsByName <- function(snnsObject, prefix) {
  
  res <- NULL
  
  nUnits <- snnsObject$getNoOfUnits()
  
  for(i in 1:nUnits)  {
    if(i==1)  unit <- snnsObject$getFirstUnit()
    else unit <- snnsObject$getNextUnit()
    
    name <- snnsObject$getUnitName(unit)
    if(beginsWith(name,prefix)) res <- c(res, unit)
  }
  
  res
}

#' Get an info header of the network.
#'  
#' @return a data frame containing some general characteristics of the network.
#' @rdname SnnsRObject-getInfoHeader
#' @name SnnsRObject$getInfoHeader
#' @usage \S4method{getInfoHeader}{SnnsR}()
#' @aliases getInfoHeader,SnnsR-method SnnsR__getInfoHeader 
SnnsR__getInfoHeader <- function(snnsObject) {
  
  NoOfUnits <- snnsObject$getNoOfUnits()
  netInfo <- snnsObject$getNetInfo()
  
  names <- NULL
  values <- NULL
  
  res <- NULL
  
  res <- rbind(res, data.frame(name=getKrioTitle(3), value=NoOfUnits, stringsAsFactors=FALSE))
  res <- rbind(res, data.frame(name=getKrioTitle(4), value=netInfo$no_of_links, stringsAsFactors=FALSE))
  res <- rbind(res, data.frame(name=getKrioTitle(5), value=netInfo$no_of_FTable_entries, stringsAsFactors=FALSE))
  res <- rbind(res, data.frame(name=getKrioTitle(6), value=netInfo$no_of_STable_entries, stringsAsFactors=FALSE))
  
  learnFunc <- snnsObject$getLearnFunc()
  updateFunc <- snnsObject$getUpdateFunc()
  
  res <- rbind(res, data.frame(name=getKrioTitle(7), value=learnFunc, stringsAsFactors=FALSE))
  res <- rbind(res, data.frame(name=getKrioTitle(16), value=updateFunc, stringsAsFactors=FALSE))
  
  if(learnFunc == "PruningFeedForward") {
    res <- rbind(res, data.frame(name=getKrioTitle(19), value=snnsObject$getPrunFunc(), stringsAsFactors=FALSE))
    res <- rbind(res, data.frame(name=getKrioTitle(20), value=snnsObject$getFFLearnFunc(), stringsAsFactors=FALSE))
  }
  
  res
}


#' Get the unit definitions of the network.
#'  
#' @return a data frame containing information about all units present in the network.
#' @rdname SnnsRObject-getUnitDefinitions
#' @name SnnsRObject$getUnitDefinitions
#' @usage \S4method{getUnitDefinitions}{SnnsR}()
#' @aliases getUnitDefinitions,SnnsR-method SnnsR__getUnitDefinitions 
SnnsR__getUnitDefinitions <- function(snnsObject) {
  
  blank <- " "
  snnsObject$getUnitDefaults()
  
  res <- NULL
  
  unit_no <- snnsObject$getFirstUnit()
  
  while(unit_no > 0)  {
    
    pos <- snnsObject$getUnitPosition( unit_no)
    
    u_name <- snnsObject$getUnitName( unit_no )
    
    if (is.null(u_name)) u_name <- blank
    
    u_type <- snnsObject$getUnitFTypeName( unit_no )
    
    if(is.null(u_type)) {
      
      no_Ftype <- TRUE
      
      u_type <- blank
      
      act_func <- snnsObject$getUnitActFuncName( unit_no )
      out_func <- snnsObject$getUnitOutFuncName( unit_no )
      
    } else {
      
      no_Ftype <- FALSE
      
      act_func <- NULL
      out_func <- NULL
      
    }
    
    if(is.null(act_func)) act_func <- blank
    if(is.null(out_func)) out_func <- blank
    
    unit_act <- snnsObject$getUnitActivation( unit_no ) 
    unit_bias <- snnsObject$getUnitBias( unit_no )
    unit_ttype <- snnsObject$getUnitTType( unit_no )
    
    unit_ttype <- getSnnsRDefine("topologicalUnitTypes", unit_ttype)
    
    
    sites <- NULL
    
    if ( no_Ftype )  {
      
      site_no <- snnsObject$setFirstSite()
      
      while(site_no > 0) {
        
        sites <- c(sites, snnsObject$getSiteName())
        site_no <- snnsObject$setNextSite()
      }
    }
    
    if(is.null(sites)) {
      siteNames <- blank
    } else {
      siteNames <- paste(sites, sep="", collapse=", ")  
    }
    
    
    res <- rbind(res, data.frame(unitNo=unit_no, unitName=u_name,
            unitAct=unit_act,
            unitBias=unit_bias,
            type=unit_ttype,
            posX=pos$x, posY=pos$y, posZ=pos$z, 
            actFunc=act_func, outFunc=out_func, sites=siteNames, 
            stringsAsFactors=FALSE))
    
    unit_no <- snnsObject$getNextUnit()  
  }
  
  res  
}


#' Get the sites definitions of the network.
#'  
#' @return a data frame containing information about all sites present in the network.
#' @rdname SnnsRObject-getSiteDefinitions
#' @name SnnsRObject$getSiteDefinitions
#' @usage \S4method{getSiteDefinitions}{SnnsR}()
#' @aliases getSiteDefinitions,SnnsR-method SnnsR__getSiteDefinitions 
SnnsR__getSiteDefinitions <- function(snnsObject) {
  
  res <- NULL
  
  tableEntry <- snnsObject$getFirstSiteTableEntry()
  
  while(tableEntry$ret) {
    
    res <- rbind(res, data.frame(siteName=tableEntry$site_name, 
            siteFunc=tableEntry$site_func, stringsAsFactors=FALSE))
    
    tableEntry <- snnsObject$getNextSiteTableEntry()
  }
  
  res  
}

#' Get the FType definitions of the network.
#'  
#' @return a data frame containing information about FType units present in the network.
#' @rdname SnnsRObject-getTypeDefinitions
#' @name SnnsRObject$getTypeDefinitions
#' @usage \S4method{getTypeDefinitions}{SnnsR}()
#' @aliases getTypeDefinitions,SnnsR-method SnnsR__getTypeDefinitions 
SnnsR__getTypeDefinitions <- function(snnsObject) {
  
  blank <- " "
  
  res <- NULL
  
  isEntry <- snnsObject$setFirstFTypeEntry()
  
  while(isEntry) {
    
    sites <- NULL
    
    isSite <- snnsObject$setFirstFTypeSite()
    
    while(isSite) {
      
      sites <- c(sites, snnsObject$getFTypeSiteName())
      isSite <- snnsObject$setNextFTypeSite()
    }
    
    if(is.null(sites)) {
      siteNames <- blank
    } else {
      siteNames <- paste(sites, sep="", collapse=", ")  
    }
    
    res <- rbind(res, data.frame(typeName=snnsObject$getFTypeName(), 
            typeActFuncName=snnsObject$getFTypeActFuncName(), 
            typeOutFuncName=snnsObject$getFTypeOutFuncName(),
            sites=siteNames,
            stringsAsFactors=FALSE))
    
    isEntry <- snnsObject$setNextFTypeEntry()
  }
  
  res  
}

SnnsR__getConnectionDefs <- function(snnsObject) {
  
  res <- NULL
  
  # this should be a reimplementation of krio_writeConnectionDefs() in kr_io.cpp
  # however, it is currently not implemented. 
  
  res
}

SnnsR__getSubnetDefs <- function(snnsObject) {
  
  res <- NULL
  
  # this should be a reimplementation of krio_writeSubnetDefs() in kr_io.cpp
  # however, it is currently not implemented. 
  
  res
}

SnnsR__getLayerDefs <- function(snnsObject) {
  
  res <- NULL
  
  # this should be a reimplementation of krio_writeLayerDefs() in kr_io.cpp
  # however, it is currently not implemented. 
  
  res
}

SnnsR__getTimeDelayDefs <- function(snnsObject) {
  
  res <- NULL
  
  # this should be a reimplementation of krio_writeTimeDelayDefs() in kr_io.cpp
  # however, it is currently not implemented. 
  
  res
}

#' Get characteristics of the network.
#' 
#' The returned list has three members: 
#' \itemize{
#' \item infoHeader general information about the network
#' \item unitDefinitions information about the units
#' \item fullWeightMatrix weight matrix of the connections
#' }
#'  
#' @return a list of data frames containing information extracted from the network.
#' @rdname SnnsRObject-extractNetInfo
#' @name SnnsRObject$extractNetInfo
#' @usage \S4method{extractNetInfo}{SnnsR}()
#' @aliases extractNetInfo,SnnsR-method SnnsR__extractNetInfo
SnnsR__extractNetInfo <- function(snnsObject) {
  
  res <- list()
  
  res[["infoHeader"]] <- snnsObject$getInfoHeader()
  res[["typeDefinitions"]] <- snnsObject$getTypeDefinitions()
  res[["unitDefinitions"]] <- snnsObject$getUnitDefinitions()
  res[["siteDefinitions"]] <- snnsObject$getSiteDefinitions()

  res[["fullWeightMatrix"]] <- snnsObject$getWeightMatrix(snnsObject$getAllUnits(), snnsObject$getAllUnits())
  
  res[["connectionDefs"]] <- snnsObject$getConnectionDefs()
  res[["subnetDefs"]] <- snnsObject$getSubnetDefs()
  res[["layerDefs"]] <- snnsObject$getLayerDefs()
  res[["timeDelayDefs"]] <- snnsObject$getTimeDelayDefs()
  
  res
}


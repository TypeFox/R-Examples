# Copyright 2006 Laboratoire de Biologie et de Biometrie Appliquée 
# (UMR 5558);CNRS; Univ. Lyon 1, 43 bd 11 nov, 69622, 
# Villeurbanne Cedex, France.
#
# This file is part of MareyMapGUI.
#
# MareyMapGUI is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# MareyMapGUI is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MareyMapGUI; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


setGeneric("updateLocalRateQuery", function(qry, fil, pos) standardGeneric("updateLocalRateQuery"))
setMethod("updateLocalRateQuery", "ANY", function(qry, fil, pos) {
  curMap <- getVar("curMap")
	curSet <- getVar("curSet")
	coll <- getVar("coll")
  map <- regEval("coll[[curSet, curMap]]")
	if((!is.null(map)) && length(interpolations(map))) {
    tkconfigure(qry, state = "normal")
    tkconfigure(fil, state = "normal")
    tkconfigure(pos, state = "normal")
  } else {
    tkconfigure(qry, state = "disabled")
    tkconfigure(fil, state = "disabled")
    tkconfigure(pos, state = "disabled")
  }
})


onSelectPositionFile <- function(resFrm) {
  if(length(regEval("interpolations(coll[[curSet, curMap]])")) == 0) {
    tkmessageBox(title = "Error", message = "   Use an interpolation method before  ")
  } else {
  	fildial <- paste("tk_getOpenFile",
      			"-filetypes { 	{\"text files\" {.txt}} \
  					{\"All files\" *} }")
    ret <- .Tcl(fildial)
  	if(tclvalue(ret) == "")
  		tkmessageBox(title = "Error", message = "   Problem with the source file  ")
  	else
      onQueryLocalRate(ret, resFrm)
	}
}


decodePosition <- function(pos) {
  sapply(pos, function(p) {
    p <- gsub("KB", "*1000", p)		
    p <- gsub("Kb", "*1000", p)
    p <- gsub("kb", "*1000", p)
    p <- gsub("MB", "*1000000", p)
    p <- gsub("Mb", "*1000000", p)
    p <- gsub("mb", "*1000000", p)
    val <- -1
    try(val <- as.numeric(eval(parse(text = p))))
	})
}


setGeneric("onQueryLocalRate", function(pos, resFrm) standardGeneric("onQueryLocalRate"))
setMethod("onQueryLocalRate", "ANY", function(pos, resFrm) {
  .Tcl(paste("eval destroy [winfo children ",	as.character(resFrm), "]"))
  readFromFile <- FALSE
  pos <- tclvalue(pos)
  if(length(regEval("interpolations(coll[[curSet, curMap]])")) == 0) {
    tkmessageBox(title = "Error", message = "Use an interpolation method before")
    
  } else {
  	if(file.exists(pos)) {
  		pos <- read.table(pos, header = TRUE)
      if(!("map" %in% names(pos)))
        tkmessageBox(title = "Error", message = "   The text file must contain a 'map' column   ")
      if(!("phys" %in% names(pos)))
        tkmessageBox(title = "Error", message = "   The text file must contain a 'phys' column   ")
      if(length(levels(pos$map)) > 1)
        tkmessageBox(title = "Caution", message = paste("Only the current map ('", getVar("curMap"),"' in '", getVar("curSet"),"') is selected in the file.", sep = ""))
      pos <- pos[pos$map == getVar("curMap"), ]
      readFromFile <- TRUE
    } else {
      # pos it is not a filename. try to use it as a ':' separated list of positions
      pos <- strsplit(pos, ":")[[1]]
      pos <- data.frame("map" = rep(getVar("curMap"), length(pos)), "phys" = pos)
    }
    
    if(!("set" %in% names(pos)))
      pos <- cbind("set" = rep(getVar("curSet"), nrow(pos)), pos)
    
    pos$phys <- decodePosition(pos$phys)
    sets <- split(pos, pos$set, drop = TRUE)
    setres <- lapply(sets, function(set) {
      maps <- split(set, set$map, drop = TRUE)
      mapres <- lapply(maps, function(mappos) {
        map <- regEval(paste("coll[[\"", as.vector(mappos$set)[[1]], "\",\"", as.vector(mappos$map)[[1]], "\"]]", sep = ""))
        if(is.null(map)) 
          return()
        if(length(interpolations(map)) == 0)
          return()
        qry <- query(map, mappos$phys)
        if(!is.matrix(qry))
          qry <- t(qry)
        cbind(mappos, qry)
      })
    
      mapmerged <- NULL
      for(i in mapres) {
        if(!is.null(i)) {
          if(is.null(mapmerged))
            mapmerged <- i
          else
            mapmerged <- merge(mapmerged, i, all = T)
        }
      }
    	mapmerged
  	})
  
  	setmerged <- NULL
    for(i in setres) {
      if(!is.null(i)) {
      	if(is.null(setmerged))
        	setmerged <- i
        else
          setmerged <- merge(setmerged, i, all = T)
      }
    }
    
    if(is.null(setmerged))
      return()
    if(!readFromFile && nrow(setmerged) == 1) {
      ## display result in the main window
      frmpos <- tkframe(resFrm)
  		tkpack(tklabel(frmpos, text = paste(setmerged$phys[[1]], "bp")))
  		tkpack(tklabel(frmpos))
      map <- regEval("coll[[curSet, curMap]]")
      for(itr in interpolations(map)) {
      	tkpack(tklabel(frmpos, text = paste(name(itr), ": ", setmerged[[name(itr)]][[1]], " cM/Mb", sep = "")))
      }
  		tkpack(tklabel(frmpos))
  		tkpack(frmpos, side = "top")
      setVar("curEst", setmerged$phys)
    } else {
      if(!is.null(getVar("curEst")))
        setVar("curEst",NULL)
      
      #open a separate window to display the results
      resdial <- tktoplevel()
      tktitle(resdial) <- "Results"
  		resyscr <- tkscrollbar(resdial, command = function(...) tkyview(restxt, ...))
  		resxscr <- tkscrollbar(resdial, orient = "horizontal", command = function(...) tkxview(restxt, ...))
  		restxt <- tktext(resdial, yscrollcommand = function(...) tkset(resyscr, ...), xscrollcommand = function(...) tkset(resxscr, ...))
      butFrm <- tkframe(resdial)  ## frame where the results are written
  		clobut <- tkbutton(butFrm, text = "Close", command = function() tkdestroy(resdial))
      savbut <- tkbutton(butFrm, text = "Save as text file", command = function() {
  			fildial <- paste("tk_getSaveFile", "-filetypes { 	{\"text files\" {.txt}} \ {\"All files\" *} }") ## 'tk_getSaveFile' open a french written window
  			ret <- .Tcl(fildial)
        if(tclvalue(ret) == "")
        	tkdestroy(resdial)
        write.table(setmerged, file = tclvalue(ret), row.names = FALSE)
  			tkdestroy(resdial)
      })    
      tkgrid(tklabel(butFrm), savbut, tklabel(butFrm), clobut, tklabel(butFrm))
      tkgrid(tklabel(resdial), restxt, resyscr, tklabel(resdial))
      tkgrid(tklabel(resdial), resxscr, tklabel(resdial), tklabel(resdial))
      tkgrid(tklabel(resdial))
      tkgrid(tklabel(resdial), butFrm, tklabel(resdial))
      tkgrid.configure(resyscr, sticky = "ns")
      tkgrid.configure(resxscr, sticky = "ew")
  		for(elt in names(setmerged)) {
        tkinsert(restxt, "end", as.character(elt))
        tkinsert(restxt, "end", "\t")
      }
  		tkinsert(restxt, "end", "\n")
      for(i in 1:nrow(setmerged)) {
        for(elt in setmerged[i, ]) {
          tkinsert(restxt, "end", as.character(elt))
          tkinsert(restxt, "end", "\t")
        }
        tkinsert(restxt, "end", "\n")
      }
    }
  }
})


setGeneric("LocalRateQuery", function(parent) standardGeneric("LocalRateQuery"))
setMethod("LocalRateQuery", "ANY", function(parent) {
	regVar("curEst", NULL)
  regVar("qryPos", NULL)
  
	rateFrm <- tkframe(parent, relief = "sunken", borderwidth = 2)
	rateTopFrm <- tkframe(rateFrm, relief = "groove", borderwidth = 2)
	tkpack(rateTopFrm, side = "top", fill = "x")
	tkpack(tklabel(rateTopFrm, text = "LOCAL \nRECOMBINATION RATE"), side = "top")
	rateCtrFrm <- tkframe(rateFrm, relief = "flat") 
	tkpack(rateCtrFrm, side = "top", fill = "x")
	rateBtmFrm <- tkframe(rateFrm, relief = "flat") 
	tkpack(rateBtmFrm, side = "top", fill = "x")
	rateBtmFrm2 <- tkframe(rateFrm, relief = "flat")
	tkpack(rateBtmFrm2, side = "top", fill = "x")
	posVar <- tclVar(0)
	posEntry <- tkentry(rateBtmFrm, width = 12, textvariable = posVar, state = "disabled")
  
	qryBut <- tkbutton(rateBtmFrm, text = "Query", state = "disabled", command = function() onQueryLocalRate(posVar, rateCtrFrm))
	tkgrid(posEntry, tklabel(rateBtmFrm, text = ""), qryBut)		
	filBut <- tkbutton(rateBtmFrm2, text = "Read positions from file", state = "disabled", command = function() onSelectPositionFile(rateCtrFrm))
	tkgrid(filBut)
	attachCallback("clearlocalrate", function() {
		.Tcl(paste("eval destroy [winfo children ",	as.character(rateCtrFrm), "]"))
		#updating curReg without using setvar("curEst") to avoid plotting twice
		regEval("curEst=NULL")
	}, "curMap")
  attachCallback("external query", function() {
    onQueryLocalRate(tclVar(getVar("qryPos")), rateCtrFrm)
  }, "qryPos")

	attachCallback("updateLocalRateQuery", function() updateLocalRateQuery(qryBut, filBut, posEntry), "curMap")
	rateFrm
})

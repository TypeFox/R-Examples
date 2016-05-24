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



setGeneric("updateMkrInfo", function(name, phys, gen, val, valvar, del, qry) standardGeneric("updateMkrInfo"))
setMethod("updateMkrInfo", "ANY", function(name, phys, gen, val, valvar, del, qry) {
  curMkr <- getVar("curMkr")
  if(curMkr == 0) {
	  tkconfigure(name, text = "")
	  tkconfigure(phys, text = "")
	  tkconfigure(gen, text = "")
	  tkconfigure(val, state = "disabled")
	  tkconfigure(del, state = "disabled")
    tkconfigure(qry, state = "disabled")
  } else {
	  coll <- getVar("coll")
	  curSet <- getVar("curSet")
	  curMap <- getVar("curMap")
	  mkr <- coll[[curSet, curMap]][[curMkr]]
	  tkconfigure(name, text = mkr[["name"]])
	  tkconfigure(phys, text = mkr[["phys"]])
	  tkconfigure(gen, text = mkr[["gen"]])
	  if(mkr[["valid"]] == T)
		  tclvalue(valvar) <- 1
	  else
		  tclvalue(valvar) <- 0
	  tkconfigure(val, state = "normal")
	  tkconfigure(del, state = "normal")
    tkconfigure(qry, state = "normal")
  }
})


setGeneric("onMkrVal", function(x, y) standardGeneric("onMkrVal"))
setMethod("onMkrVal", c("character", "character"), function(x, y) {
	coll <- getVar("coll")
	curSet <- getVar("curSet")
	curMap <- getVar("curMap")
	curMkr <- getVar("curMkr")
	if(curMkr == "0") {
    tkmessageBox(title = "Error", message = "   Please select one marker on the map   ")
  } else {
  	mkr <- coll[[curSet, curMap]][[curMkr]]
  	if(mkr[["valid"]])
  		mkr[["valid"]] <- F
  	else
  		mkr[["valid"]] <- T
  	coll[[curSet, curMap]][[curMkr]] <- mkr
  	setVar("coll", coll)
  	setVar("curMkr", 0)
	}
})


setGeneric("onMkrDel", function(dummy) standardGeneric("onMkrDel"))

setMethod("onMkrDel", "character", function(dummy){
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Confirm")
	frm1 <- tkframe(dialog)
	msglbl <- tklabel(frm1, text = "  Informations about this marker will be   \n\ permanently lost, continue ?")
	tkpack(msglbl, side = "left")
	
	frm3 <- tkframe(dialog)
	butno <- tkbutton(frm3, text = "No", command = function() {tkdestroy(dialog)})
	butyes <- tkbutton(frm3, text = "Yes", command = function() {
		getVar("curMkr")
		curSet <- getVar("curSet")
		curMap <- getVar("curMap")
		regEval("coll <- coll + removeMarker(coll[[curSet, curMap]], curMkr)")
		setVar("curMkr", 0)
		tkdestroy(dialog)
	})
	tkpack(butno, side = "left")
	tkpack(butyes, side = "right")
	
	tkpack(tklabel(dialog), frm1, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2, fill = "x")
})


onMkrQuery <- function(){
  mkr <- regEval("coll[[curSet, curMap]][[curMkr]]")
  setVar("qryPos", mkr[["phys"]])
}


setGeneric("MarkerInfo", function(parent) standardGeneric("MarkerInfo"))

setMethod("MarkerInfo", "ANY", function(parent) {
	frm <- tkframe(parent, relief = "sunken", borderwidth = 2)
	markerLbl <- tklabel(frm, text = "MARKER", relief = "groove", borderwidth = 2, width = 30)
	tkpack(markerLbl, side = "top", fill = "x")

	frm1 <- tkframe(frm)
	mkrName <- tklabel(frm1, text = "")
	tkpack(tklabel(frm1, text = "Name: "), mkrName, side = "left")
	tkpack(frm1, side = "top", fill = "x")
	
	frm2 <- tkframe(frm)
	mkrPhys <- tklabel(frm2, text = "")
	tkpack(tklabel(frm2, text = "Physical Position: "), mkrPhys, side = "left")
	tkpack(frm2, side = "top", fill = "x")
	
	frm3 <- tkframe(frm)
	mkrGen <- tklabel(frm3, text = "")
	tkpack(tklabel(frm3, text = "Genetical Position: "), mkrGen, side = "left")
	tkpack(frm3, side = "top", fill = "x")
	
	frm4 <- tkframe(frm)
	mkrValVar <- tclVar("0")
	mkrVal <- tkcheckbutton(frm4, state = "disabled", variable = mkrValVar)
	tkbind(mkrVal, "<ButtonRelease-1>", onMkrVal)
	tkpack(tklabel(frm4, text = "Valid: "), mkrVal, side = "left")
	tkpack(frm4, side = "top", fill = "x")
	
	mkrDel <- tkbutton(frm, text = "Delete", state = "disabled", command = onMkrDel)
	tkpack(mkrDel, side = "top")

  mkrQry <- tkbutton(frm, text = "Query recombination rate", state = "disabled", command = onMkrQuery)
  tkpack(mkrQry, side = "top")

	attachCallback("updateMkrInfo", function() updateMkrInfo(mkrName, mkrPhys, mkrGen, mkrVal, mkrValVar, mkrDel, mkrQry), "curMkr")
	frm
})


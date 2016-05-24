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


setGeneric("startMareyMapGUI", function(dummy) standardGeneric("startMareyMapGUI"))

setMethod("startMareyMapGUI", "missing", function(dummy) {
    
  regVar("curSet", "")
  regVar("curMap", "")
  regVar("curMkr", 0)
  mareymaps <- MapCollection()
  
	## top level
	mainWin <- tktoplevel(height = 700, width = 1000)
	tkwm.title(mainWin, "Marey Map")	
	
  regVar("coll", mareymaps)
  tkconfigure(mainWin, menu = MenuBar(mainWin))
	tkgrid(tklabel(mainWin, text = ""))
  tkgrid(namlbl <- tklabel(mainWin, text = "Load a data set (File - Open)"), columnspan = 7)
  
	## left frame : 'leftFrm'
	leftFrm <- tkframe(mainWin, relief = "flat")
	tkpack(MapList(leftFrm), side = "top", fill = "x")
	tkpack(tklabel(leftFrm, text = ""), side = "top")
	tkpack(MarkerInfo(leftFrm), side = "top")
	
  ## right frame : 'rightFrm'
	rightFrm <- tkframe(mainWin, relief = "flat")
	tkpack(InterpolationList(rightFrm), side = "top")
	tkpack(tklabel(rightFrm, text = ""), side = "top")
	tkpack(LocalRateQuery(rightFrm), side = "top", fill = "x")
	
	## arange widgets in a grid
	tkgrid(tklabel(mainWin, text = ""))
	tkgrid(tklabel(mainWin, text = ""), leftFrm, tklabel(mainWin, text = ""),	PlotArea(mainWin), tklabel(mainWin, text = ""), rightFrm, tklabel(mainWin, text = ""))
	tkgrid(tklabel(mainWin, text = ""))
  
	tkgrid.configure(leftFrm, sticky = "ns")  ## the frame will be stretched to fill the entire height (n=north, s=south)
	tkgrid.configure(rightFrm, sticky = "ns") ## the frame will be stretched to fill the entire height (n=north, s=south)
  
  ## 'attachCallback' is a method created in 'Registry.R'
	attachCallback("updateSpcLbl", function() {tkconfigure(namlbl, text = getVar("curSet"))}, "curSet")
	
	attachCallback("updateMapLbl", function() {
    map <- regEval("coll[[curSet, curMap]]") 
    if(length(map) == 0) {
      if(getVar("curSet") != "")
        tkconfigure(namlbl, text = getVar("curSet"))
      else
        tkconfigure(namlbl, text = "Please select a set from the menu ('Data')")
    } else {
	    tkconfigure(namlbl, text = paste(setName(map), " - ", mapName(map), ", ", length(physicalPositions(map)), " markers", sep = ""))
    }
  }, "curMkr")


})

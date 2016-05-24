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



setGeneric("buildMapList", function(maplst) standardGeneric("buildMapList"))
setMethod("buildMapList", "ANY", function(maplst) {
	tkdelete(maplst, 0, tksize(maplst))
	curSet <- getVar("curSet")
	curMap <- getVar("curMap")
	coll <- getVar("coll")
	if(nchar(curSet) != 0) {
		lapply(regEval("mapNames(coll[[curSet]])"), function(x) {
			tkinsert(maplst, "end", x)
			if(curMap == x)
				tkselection.set(maplst, "end")
		})
	}
})


setGeneric("MapList", function(parent) standardGeneric("MapList"))
setMethod("MapList", "ANY", function(parent) {
  
	frm <- tkframe(parent, relief = "sunken", borderwidth = 2)
	mapLbl <- tklabel(frm, text = "MAPS", relief = "groove", borderwidth = 2)
	mapscr <- tkscrollbar(frm, repeatinterval = 5, command = function(...) tkyview(maplst, ...))
	maplst <- tklistbox(frm, selectmode = "single", yscrollcommand = function(...) tkset(mapscr,...), borderwidth = 0, width = 30)
	tkpack(mapLbl, side = "top", fill = "x")
	tkpack(maplst, side = "left", fill = "both")
	tkpack(mapscr, side = "right", fill = "y")
	tkbind(maplst, "<ButtonRelease-1>", function() {
    curSet <- getVar("curSet")
    if(curSet != "") {
      sel <- as.character(tclvalue(tkget(maplst, tkcurselection(maplst))))
	    if(sel != getVar("curMap")) {
		    setVar("curMap", sel)
		    setVar("curMkr", 0)
	    }
    } else
      tkmessageBox(title = "Error", message = "   Please select one set on the collection   ")
	})
	
	attachCallback("buildMapList", function() buildMapList(maplst), "curSet")
	attachCallback("buildMapList", function() buildMapList(maplst), "curMap") ## etait commente
	attachCallback("buildMapList", function() buildMapList(maplst), "coll")
	frm
})


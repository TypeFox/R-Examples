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


#--- onSave
setGeneric("onSave", function(dummy) standardGeneric("onSave"))
setMethod("onSave", "character", function(dummy) {
  coll <- getVar("coll")
  if(length(coll) == 0) { ## Never occurs
    tkmessageBox(title = "Error", message =  "  No collection opened to save.   
Load a collection previously.")
    
  } else {
		dialog <- tktoplevel()
  	tkgrab.set(dialog)
  	tkwm.title(dialog, "Save")
    
  	curSet <- getVar("curSet")
  	curMap <- getVar("curMap")
  
    ## frame1
  	frm1 <- tkframe(dialog)
  	filepath <- tclVar("")
  	entfil <- tkentry(frm1, width = "20", textvariable = filepath)
  	butbro <- tkbutton(frm1, text = "Browse", command = function() {
  		fildial <- paste("tk_getSaveFile","-filetypes { \
  						{\"R binary files\" {.rda}} \
  						{\"text files\" {.txt}} }")
  		tclvalue(filepath) <- .Tcl(fildial)
  	})
  	tkpack(entfil, side = "left")
  	tkpack(tklabel(frm1), side = "left")
  	tkpack(butbro, side = "left")
  
    ## frame2
  	frm2 <- tkframe(dialog)
  	rbValue <- tclVar("coll")
    
    rbcoll <- tkradiobutton(frm2)
    tkconfigure(rbcoll, variable = rbValue, value = "coll", text = "All data")
    tkpack(rbcoll, side = "left")
    tkpack(tklabel(frm2), side = "left")
    
  	if(nchar(curSet) > 0) {
    	rbValue <- tclVar(curSet)
      rbspec <- tkradiobutton(frm2)
      tkconfigure(rbcoll, variable = rbValue, value = "coll", text = "All data")
  		tkconfigure(rbspec, variable = rbValue, value = curSet, text = curSet)
      tkpack(rbspec, side = "left")
      
  		if(nchar(curMap) > 0) {
        rbValue <- tclVar(curMap)
  			rbmap <- tkradiobutton(frm2)
        tkconfigure(rbcoll, variable = rbValue, value = "coll", text = "All data")
  		  tkconfigure(rbspec, variable = rbValue, value = curSet, text = curSet)
  			tkconfigure(rbmap, variable = rbValue, value = curMap, text = curMap)
  			tkpack(rbmap, side = "left")
  		}
  	}
  
    ## frame3
  	frm3 <- tkframe(dialog)
  	butcan <- tkbutton(frm3, text = "Cancel", command = function() {tkdestroy(dialog)})
  	butsav <- tkbutton(frm3, text = "Save", command = function() {
      path <- tclvalue(filepath)
      if(path == "") {
        tkmessageBox(title = "Error", message = "   Please fill a file path   ")
      } else {
        ext <- substr(path, nchar(path) - 2, nchar(path))
        curSet <- getVar("curSet")
        curMap <- getVar("curMap")
        val <- tclvalue(rbValue)
        
        if(ext == "txt") {
        	if(val == curSet) {
        		set <- regEval("coll[[curSet]]")
        		textFile(set, path)
        	} else if(val == curMap) {
        		map <- regEval("coll[[curSet, curMap]]")
        		textFile(map, path)
    			} else if(val == "coll") {
        		coll <- getVar("coll")
        		textFile(coll, path)
        	}
          
        } else if(ext == "rda") {
    			if(val == curSet) {
    				set <- regEval("coll[[curSet]]")
    				save("set", file = path)
    			} else if(val == curMap) {
    				map <- regEval("coll[[curSet,curMap]]")
    				save("map", file = path)
    			} else if(val == "coll") {
    				coll <- getVar("coll")
    				save("coll", file = path)
    			}
          
    		} else
    			tkmessageBox(title = "Error", message = "   File extension must be 'rda' or 'txt'  ")
        
        tkdestroy(dialog)
      }
  	})
  	tkpack(butcan, side = "left")
  	tkpack(butsav, side = "right")
  	tkpack(tklabel(dialog), frm1, tklabel(dialog), pady = 0.2)
  	tkpack(tklabel(dialog), frm2, tklabel(dialog), pady = 0.2)
  	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2, fill = "x")
  }
})


#--- onOpen
setGeneric("onOpen", function(dummy) standardGeneric("onOpen"))
setMethod("onOpen", "character", function(dummy) {
	fildial <- paste("tk_getOpenFile", "-initialdir MareyMap/data",
    			"-filetypes { 	{\"R binary files\" {.rda .Rdata .Rda}} \
					{\"text files\" {.txt}}}")
  filename <- tclvalue(.Tcl(fildial))
  if(filename == "") {
    tkmessageBox(title = "Error", message = "No file was selected!")
    stop("Process killed because there is no data file dowloaded. Please restart.", call. = FALSE)
  } else {
    ext <- file_ext(filename)
    if(ext == "txt")
		  readTextFile(filename)
    else
      readRDataFile(filename)
  }
})


#--- onExport
setGeneric("onExport", function(dummy) standardGeneric("onExport"))
setMethod("onExport", "character", function(dummy) {
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Export")
  
  ## frame1
	frm1 <- tkframe(dialog)
	filepath <- tclVar("")
	entfil <- tkentry(frm1, width = "20", textvariable = filepath)
	butbro <- tkbutton(frm1, text = "Browse", command = function() {
		fildial <- paste("tk_getSaveFile", "-filetypes { \
			{\".eps\" {.eps}} \
			{\".pdf\" {.pdf}} \
			{\".png\" {.png}} \
			{\".jpeg\" {.jpg}} }")
		tclvalue(filepath) <- .Tcl(fildial)
	})

	tkpack(entfil, side = "left")
	tkpack(tklabel(frm1), side = "left")
	tkpack(butbro, side = "left")

  ## frame2
	frm2 <- tkframe(dialog)
	rbmap <- tkradiobutton(frm2)
	rbrate <- tkradiobutton(frm2)
	rbboth <- tkradiobutton(frm2)

  if(length(regEval("interpolations(coll[[curSet, curMap]])")) == 0) {
    rbval <- tclVar("map")
	  tkconfigure(rbmap, variable = rbval, value = "map", text = "Marey Map")
	  tkconfigure(rbrate, variable = rbval, value = "rates", text = "Recombination Rates", state = "disabled")
	  tkconfigure(rbboth, variable = rbval, value = "both", text = "Both", state = "disabled")
  } else {
    rbval <- tclVar("both")
	  tkconfigure(rbmap, variable = rbval, value = "map", text = "Marey Map")
	  tkconfigure(rbrate, variable = rbval, value = "rates", text = "Recombination Rates")
	  tkconfigure(rbboth, variable = rbval, value = "both", text = "Both")
  }
    
	tkpack(rbmap, side = "left")
	tkpack(tklabel(frm2), side = "left")
	tkpack(rbrate, side = "left")
	tkpack(tklabel(frm2), side = "left")
	tkpack(rbboth, side = "left")

  ## frame3
	frm3 <- tkframe(dialog)
	butcan <- tkbutton(frm3, text = "Cancel", command = function() {tkdestroy(dialog)})
	butsav <- tkbutton(frm3, text = "Export", command = function() {
		path <- tclvalue(filepath)
		if(path == "")
			tkmessageBox(title = "Error", message = "   Please fill a file path   ")
		else {
      ext <- substr(path, nchar(path) - 2, nchar(path))
  		if(tclvalue(rbval) == "map") {
  			if(ext == "eps")
  				postscript(file = path, horizontal = FALSE, onefile = FALSE, width = 6, height = 4.8, paper = "special", title = paste(getVar("curSet"), " - ", getVar("curMap"), sep = ""))
  			else if(ext == "pdf")
  				pdf(file = path, width = 6, height = 4.8, title = paste(getVar("curSet"), " - ", getVar("curMap"), sep = ""))
  			else if(ext == "png")
  				png(filename = path, width = 480, height = 384)
  			else if(ext == "jpg")
  				jpeg(filename = path, width = 480, height = 384)
  			else
  				tkmessageBox(title = "Error", message = paste("   The file extension ", ext, " is not allowed   ", sep = ""))
  			drawPlt1()
  			dev.off()
  		}
  		if(tclvalue(rbval) == "rates") {
  			if(ext == "eps")
  				postscript(file = path, horizontal = FALSE, onefile = FALSE, width = 6, height = 2, paper = "special", title = paste(getVar("curSet"), " - ", getVar("curMap"), sep = ""))
  			else if(ext == "pdf")
  				pdf(file = path, width = 6, height = 2, title = paste(getVar("curSet"), " - ", getVar("curMap"), sep = ""))
  			else if(ext == "png")
  				png(filename = path, width = 480, height = 144)
  			else if(ext=="jpg")
  				jpeg(filename = path, width = 480, height = 144)
  			else
  				tkmessageBox(title = "Error", message = paste("   The file extension ", ext, " is not allowed   ", sep = ""))
  			drawPlt2()
  			dev.off()
  		}
  		if(tclvalue(rbval) == "both") {
  			if(ext == "eps")
  				postscript(file = path, horizontal = FALSE, onefile = FALSE, width = 6, height = 7.2, paper = "special", title = paste(getVar("curSet"), " - ", getVar("curMap"), sep = ""))
  			else if(ext == "pdf")
  				pdf(file = path, width = 6, height = 7.2, title = paste(getVar("curSet"), " - ", getVar("curMap"), sep = ""))
  			else if(ext == "png")
  				png(filename = path, width = 480, height = 576)
  			else if(ext == "jpg")
  				jpeg(filename = path, width = 480, height = 576)
  			else
  				tkmessageBox(title = "Error", message = paste("   The file extension ", ext, " is not allowed   ", sep = ""))
  			lo <- layout(matrix(c(1, 2), 2), heights = c(2, 1))
  			drawPlt1()
  			drawPlt2()
  			dev.off()
  		}
  		tkdestroy(dialog)
    }
	})
	tkpack(butcan, side = "left")
	tkpack(butsav, side = "right")
	tkpack(tklabel(dialog), frm1, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm2, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2, fill = "x")
})


#--- readTextFile
setGeneric("readTextFile", function(file) standardGeneric("readTextFile"))
setMethod("readTextFile", "character", function(file) {
	ll <- readLines(file)
	ll <- ll[grep("@@@", ll)]
	ll <- lapply(ll, function(x) {
    x <- unlist(strsplit(x, "@@@"))
  })

	tab <- read.table(as.character(file), header = T)
  if(!"set" %in% names(tab))
    tkmessageBox(title = "Error", message = "   A 'set' column is missing   ")
  if(!"map" %in% names(tab))
    tkmessageBox(title = "Error", message = "   A 'map' column is missing   ")
  if(!"mkr" %in% names(tab))
    tkmessageBox(title = "Error", message = "   A 'mkr' column is missing   ")
  if(!"phys" %in% names(tab))
    tkmessageBox(title = "Error", message = "   A 'phys' column is missing   ")
  if(!"gen" %in% names(tab))
    tkmessageBox(title = "Error", message = "   A 'gen' column is missing   ")
  if(!"valid" %in% names(tab))
    tab$valid <- rep("TRUE", nrow(tab))
  
	sets <- split(tab, tab$set, drop = TRUE)
	coll <- getVar("coll")
	lapply(sets, function(set) {
		maps <- split(set, set$map, drop = TRUE)
		lapply(maps, function(map) {
			newmap <- as(map, "MareyMap")
			lapply(ll, function(x) {
				if(x[[2]] == setName(newmap) && x[[3]] == mapName(newmap)) {
					eval(parse(text = paste("itr<-", x[[5]])))
					newmap <<- newmap + x[[5]] 
				}
			})
			coll <<- coll + newmap
		})
	})
	setVar("coll", coll)
})


onRenameMap <- function() {
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Rename map")
  
  ## frame1
	frm1 <- tkframe(dialog)
	filepath <- tclVar(getVar("curMap"))
	entfil <- tkentry(frm1, width = "20", textvariable = filepath)
	tkpack(tklabel(frm1, text = "name: "), entfil, side = "left")

  ## frame3
	frm3 <- tkframe(dialog)
	butcan <- tkbutton(frm3, text = "Cancel", command = function() {tkdestroy(dialog)})
	butsav <- tkbutton(frm3, text = "Rename", command = function() {
		nn <- tclvalue(filepath)
		coll <- getVar("coll")
		curMap <- getVar("curMap")
		curSet <- getVar("curSet")
		set <- coll[[curSet]]
		map <- set[[curMap]]
		set <- set - curMap
		mapName(map)<- nn
		set <- set + map
		coll <- coll + set
		setVar("coll", coll)
		setVar("curMap", nn)
		setVar("curMkr", 0)
		tkdestroy(dialog)
	})

	tkpack(butcan, tklabel(frm3), side = "left")
	tkpack(butsav, side = "right")
	tkpack(tklabel(dialog), frm1, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2)
}


onRemoveMap <- function() {
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Confirm")
  
  ## frame1
	frm1 <- tkframe(dialog)
	msglbl <- tklabel(frm1, text = paste("   Remove the map '", getVar("curMap"), "' from the collection ?   ", sep = ""))
	tkpack(msglbl, side = "left")

  ## frame3
	frm3 <- tkframe(dialog)
	butno <- tkbutton(frm3, text = "No", command = function() {tkdestroy(dialog)})
	butyes <- tkbutton(frm3, text = "Yes", command = function() {
		coll <- getVar("coll")
		curMap <- getVar("curMap")
		curSet <- getVar("curSet")
		set <- coll[[curSet]]
		set <- set - curMap
		coll <- coll + set
		setVar("coll", coll)
		setVar("curMap", "")
		setVar("curMkr", 0)
		tkdestroy(dialog)
	})

	tkpack(butno, side = "left")
	tkpack(butyes,side = "right")
	tkpack(tklabel(dialog),	frm1, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2, fill = "x")
}


onRenameSet <- function() {
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Rename set")
  
  ## frame1
	frm1 <- tkframe(dialog)
	filepath <- tclVar(getVar("curMap"))
	entfil <- tkentry(frm1, width = "20", textvariable = filepath)
	tkpack(tklabel(frm1, text = "name: "), entfil, side = "left")

  ## frame3
	frm3 <- tkframe(dialog)
	butcan <- tkbutton(frm3, text = "Cancel", command = function() {tkdestroy(dialog)})
	butsav <- tkbutton(frm3, text = "Rename", command = function() {
		nn <- tclvalue(filepath)
		coll <- getVar("coll")
		curSet <- getVar("curSet")
		set <- coll[[curSet]]
    coll  <- coll - curSet
		setName(set)<- nn
		coll <- coll + set
		setVar("coll", coll)
    setVar("curSet", nn)
    setVar("curMap", getVar("curMap"))
		tkdestroy(dialog)
	})

	tkpack(butcan, tklabel(frm3), side = "left")
	tkpack(butsav, side = "right")
	tkpack(tklabel(dialog), frm1, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2)
}


onRemoveSet <- function() {
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Confirm")
  
  ## frame1
	frm1 <- tkframe(dialog)
	msglbl <- tklabel(frm1, text = paste("   Remove the set '", getVar("curSet"), "' from the collection ?   ", sep = ""))
	tkpack(msglbl, side = "left")

  ## frame3
	frm3 <- tkframe(dialog)
	butno <- tkbutton(frm3, text = "No", command = function() {tkdestroy(dialog)})
	butyes <- tkbutton(frm3, text = "Yes", command = function() {
		coll <- getVar("coll")
		curMap <- getVar("curMap")
		curSet <- getVar("curSet")
		coll <- coll - curSet
    setVar("curSet", "")
		setVar("curMap", "")
		setVar("curMkr", 0)
		setVar("coll", coll)
    tkdestroy(dialog)
	})

	tkpack(butno, side = "left")
	tkpack(butyes,side = "right")
	tkpack(tklabel(dialog), frm1, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2, fill = "x")
}


onClearCollection <- function() {
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Confirm")
  
  ## frame1
	frm1 <- tkframe(dialog)
	msglbl <- tklabel(frm1, text = "   Clear the collection ?   ")
	tkpack(msglbl, side = "left")

  ## frame3
	frm3 <- tkframe(dialog)
	butno <- tkbutton(frm3, text = "No", command = function() {tkdestroy(dialog)})
	butyes <- tkbutton(frm3, text = "Yes", command = function() {
    setVar("curSet", "")
		setVar("curMap", "")
		setVar("curMkr", 0)
		setVar("coll", MapCollection())
    tkdestroy(dialog)
	})

	tkpack(butno, side = "left")
	tkpack(butyes, side = "right")
	tkpack(tklabel(dialog), frm1, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2, fill = "x")
}


#--- readRDataFile
setGeneric("readRDataFile", function(file) standardGeneric("readRDataFile"))
setMethod("readRDataFile", "character", function(file) {
	objname <- load(file)
	on.exit(return)
	obj <- get(objname)
  
	if(length(objname) == 1) {
		if(is(obj, "MapCollection")) {
      if(length(getVar("coll")) == 0) { ## no collection already opened
        setVar("coll", obj)
        tkmessageBox(title = "Close", message = paste("Collection '", objname, "' has been added.", sep = ""))
      } else {
        dialog <- tktoplevel()
			  tkwm.title(dialog, "")
			  frmbut <- tkframe(dialog)
			  merbut <- tkbutton(frmbut, text = "Merge", command = function() {
          coll <- getVar("coll") + obj
          setVar("coll", coll)
          tkdestroy(dialog)
        })
			  repbut <- tkbutton(frmbut, text = "Replace", command = function() {
          setVar("coll", obj)
          tkdestroy(dialog)
        })
			  tkpack(tklabel(frmbut, text = ""), merbut, side = "left")
			  tkpack(repbut, tklabel(frmbut, text = ""), side = "right")
			  tkpack(tklabel(dialog, text = ""))
			  tkpack(tklabel(dialog, text = "The file contains a collection of maps.\
   Do you want this colection to replace the current   \
one or merge the two collections together ?"))
			  tkpack(tklabel(dialog, text = ""))
			  tkpack(frmbut, fill = "x")
			  tkpack(tklabel(dialog, text = ""))
			  return()
      }

		} else if(is(obj, "MapSet")) {
			coll <- getVar("coll") + obj
			setVar("coll", coll)
			setVar("curSet", "")
			setVar("curMap", "")
			setVar("curMkr", 0)
      tkmessageBox(title = "Close", message = paste("Species \"", setName(obj),	"\" has been added to the collection.", sep = ""))

		} else if(is(obj, "MareyMap")) {
			coll <- getVar("coll") + obj
			setVar("coll", coll)
			setVar("curSet", "")
			setVar("curMap", "")
			setVar("curMkr", 0)
      tkmessageBox(title = "Close", message = paste("Map \"", setName(obj), ", ", mapName(obj), "\" has been added to the collection.", sep = ""))
      
		} else {
      tkmessageBox(title = "Error", message = "Invalid object type")
    }
	}
})


setGeneric("onLoadCollRda", function(dummy) standardGeneric("onLoadCollRda"))
setMethod("onLoadCollRda", "character", function(dummy) {
	fildial <- paste("tk_getOpenFile", "-filetypes { {\"R binary files\" {.rda .Rdata .Rda}} {{All files} *} }")
	ret <- .Tcl(fildial)
	if(nchar(tclvalue(ret)) == 0)
		return
	objname <- load(tclvalue(ret))
	obj <- get(objname)
	if(!is(obj, "MapCollection")) {
    dialog <- tktoplevel()
		tkwm.title(dialog, "Error")
		clobut <- tkbutton(dialog, text = "Close", command = function() tkdestroy(dialog))
		tkgrid(tklabel(dialog, text = ""))
		tkgrid(tklabel(dialog, text = "The selected file does not contain an object of type MapCollection"), padx = 10)
		tkgrid(tklabel(dialog, text = ""))
		tkgrid(clobut)
		tkgrid(tklabel(dialog, text = ""))
	} else {
		setVar("coll", obj)
		setVar("curSet", "")
		setVar("curMap", "")
		setVar("curMkr", 0)
		dialog <- tktoplevel()
		tkwm.title(dialog, "")
		clobut <- tkbutton(dialog, text = "Close", command = function() tkdestroy(dialog))
		tkgrid(tklabel(dialog, text = ""))
		tkgrid(tklabel(dialog, text = "Your working map collection has been changed."), padx = 10)
		tkgrid(tklabel(dialog, text = ""))
		tkgrid(clobut)
		tkgrid(tklabel(dialog, text = ""))
	}
})


setGeneric("onQuit", function(win) standardGeneric("onQuit"))
setMethod("onQuit", "ANY", function(win) { 
  dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Quit")
      
  frm1 <- tkframe(dialog)
  msglbl <- tklabel(frm1, text = "   Close the software ?   ")
	tkpack(msglbl, side = "left")
    
  frm3 <- tkframe(dialog)
  butno <- tkbutton(frm3, text = "No", command = function() {tkdestroy(dialog)})
  butyes <- tkbutton(frm3, text = "Yes", command = function() {
  	tkdestroy(dialog)
    tkdestroy(win)
  })
    
	tkpack(butno, side = "left")
	tkpack(butyes, side = "right")
	tkpack(tklabel(dialog), frm1, tklabel(dialog), pady = 0.2)
	tkpack(tklabel(dialog), frm3, tklabel(dialog), pady = 0.2, fill = "x")
})


onAbout <- function() {
	conn <- file(system.file("about.txt", package = "MareyMap"))
	abtdial <- tktoplevel()
	tkwm.title(abtdial, "About MareyMap")
	restxt <- tktext(abtdial)
	clobut <- tkbutton(abtdial, text = "Close", command = function() tkdestroy(abtdial))
	tkgrid(tklabel(abtdial))
	tkgrid(tklabel(abtdial), restxt, tklabel(abtdial))
	tkgrid(tklabel(abtdial))
	tkgrid(tklabel(abtdial), clobut, tklabel(abtdial))
	for(line in readLines(conn)) {
		tkinsert(restxt, "end", line)
		tkinsert(restxt, "end", "\n")
	}
	tkconfigure(restxt, state = "disabled")
}


onLicense <- function() {
	conn <- file(system.file("license.txt", package = "MareyMap"))
	abtdial <- tktoplevel()
	tkwm.title(abtdial, "License")
  restxt <- tktext(abtdial)
	clobut <- tkbutton(abtdial, text = "Close", command = function() tkdestroy(abtdial))
	tkgrid(tklabel(abtdial))
	tkgrid(tklabel(abtdial), restxt, tklabel(abtdial))
	tkgrid(tklabel(abtdial))
	tkgrid(tklabel(abtdial), clobut, tklabel(abtdial))
	for(line in readLines(conn)) {
		tkinsert(restxt, "end", line)
		tkinsert(restxt, "end", "\n")
	}
	tkconfigure(restxt, state = "disabled")
}


setGeneric("buildSetMenu", function(setMenu) standardGeneric("buildSetMenu"))
setMethod("buildSetMenu", "ANY", function(setMenu) {
	tkdelete(setMenu, 0, "end")
	coll <- getVar("coll")
	curSet <- getVar("curSet")
	names <- regEval("setNames(coll)")
	lapply(names, function(x) {
		tkadd(setMenu, "command", label = x, command = function() {
			if(x != getVar("curSet")) {
				setVar("curSet", x)
				setVar("curMap", "")
				setVar("curMkr", 0)
			}
		})
	})		
})


setGeneric("MenuBar", function(parent) standardGeneric("MenuBar"))
setMethod("MenuBar", "ANY", function(parent) {

  menubar <- tkmenu(parent)
  
# Data menu
	fileMenu <- tkmenu(menubar, tearoff = FALSE)
	tkadd(menubar, "cascade", label = "File", menu = fileMenu)
	tkadd(fileMenu, "command", label = "Open", command = onOpen)
	tkadd(fileMenu, "command", label = "Save", command = onSave, state = "disabled")
	tkadd(fileMenu, "command", label = "Export", command = onExport, state = "disabled")
	sep <- tkadd(fileMenu, "separator")
	tkadd(fileMenu, "command", label = "Quit", command = function() onQuit(parent))
  
# Edit menu
	editMenu <- tkmenu(menubar, tearoff = FALSE )
	tkadd(menubar, "cascade", label = "Edit", menu = editMenu)
	tkadd(editMenu, "command", label = "Rename set", command = onRenameSet, state = "disabled")
	tkadd(editMenu, "command", label = "Remove set", command = onRemoveSet, state = "disabled")
	sep <- tkadd(editMenu, "separator")
	tkadd(editMenu, "command", label = "Rename map", command = onRenameMap, state = "disabled")
	tkadd(editMenu, "command", label = "Remove map", command = onRemoveMap, state = "disabled")
	sep <- tkadd(editMenu, "separator")
	tkadd(editMenu, "command", label = "Clear map collection", command = onClearCollection, state = "disabled")
		
# Sets menu
	setMenu <- tkmenu(menubar, tearoff = FALSE)
	tkadd(menubar, "cascade", label = "Data", menu = setMenu)
	buildSetMenu(setMenu)

# Help menu
	hlpMenu <- tkmenu(menubar, tearoff = FALSE)
	tkadd(menubar, "cascade", label = "Help", menu = hlpMenu)
	tkadd(hlpMenu, "command", label = "Using MareyMap pdf tutorial", command = function() {
		f <- system.file("Using_MareyMap.pdf", package = "MareyMap")
		viewer <- getOption("pdfviewer")
		if(!is.null(viewer)) {
			if(file.exists(viewer))
				system(paste(viewer, f, "&"))
			else
				tkmessageBox(title = "Error", message = paste("your\
					current pdf viewer: ", viewer, "does not exist.
					You can change it by setting the variable
					R_PDFVIEWER in you ~/.Renviron file to a correct location.
					Ex: R_PDFVIEWER=/usr/bin/evince"))
		}
	})
	if(is.null(getOption("pdfviewer")))
		.Tcl(paste(fileMenu, " entryconfigure ", 0," -state disabled", sep = ""))
	tkadd(hlpMenu, "command", label = "About", command = onAbout)
	tkadd(hlpMenu, "command", label = "License", command = onLicense)

# Maps menu
	attachCallback("updateMapsAndSetMenu", function() {
		buildSetMenu(setMenu)
	}, "curSet")
  
	attachCallback("updateMapsAndSetMenu", function() {
		buildSetMenu(setMenu)
	}, "coll")
  
	attachCallback("activecoll", function() {
		if(length(getVar("coll")) == 0) {
			.Tcl(paste(editMenu, " entryconfigure ", 6," -state disabled", sep = ""))
      .Tcl(paste(fileMenu, " entryconfigure ", 1," -state disabled", sep = ""))
		} else {
			.Tcl(paste(editMenu, " entryconfigure ", 6," -state normal", sep = ""))
      .Tcl(paste(fileMenu, " entryconfigure ", 1," -state normal", sep = ""))
    }
	}, "coll")

	attachCallback("activeMap", function() {
		if(getVar("curMap") != "") {
			.Tcl(paste(fileMenu, " entryconfigure ", 2," -state normal", sep = ""))
			.Tcl(paste(editMenu, " entryconfigure ", 3," -state normal", sep = ""))
			.Tcl(paste(editMenu, " entryconfigure ", 4," -state normal", sep = ""))
		} else {
			.Tcl(paste(fileMenu, " entryconfigure ", 2," -state disabled", sep = ""))
			.Tcl(paste(editMenu, " entryconfigure ", 3," -state disabled", sep = ""))
			.Tcl(paste(editMenu, " entryconfigure ", 4," -state disabled", sep = ""))
		}
	}, "curMap")

	attachCallback("activeSet", function() {
		if(getVar("curSet") != "") {
			.Tcl(paste(editMenu, " entryconfigure ", 0," -state normal", sep = ""))
			.Tcl(paste(editMenu, " entryconfigure ", 1," -state normal", sep = ""))
		} else {
			.Tcl(paste(editMenu, " entryconfigure ", 0," -state disabled", sep = ""))
			.Tcl(paste(editMenu, " entryconfigure ", 1," -state disabled", sep = ""))
		}
	}, "curSet")

	menubar
})


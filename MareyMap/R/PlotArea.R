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



setGeneric("toRCoords", function(tkCoords) standardGeneric("toRCoords"))
setMethod("toRCoords", "vector", function(tkCoords) {
	plotParam <- getVar("plotParam")
	plotArea <- getVar("plotArea")
	
	usr <- plotParam$usr
	plt <- plotParam$plt
	
	tk_width <- as.numeric(tclvalue(tkwinfo("reqwidth", plotArea)))
	tk_height <- as.numeric(tclvalue(tkwinfo("reqheight", plotArea)))
	
	tk_xMouse <- as.numeric(tkCoords[[1]])
	tk_yMouse <- tk_height - as.numeric(tkCoords[[2]])
	
	tk_xOrig <- plt[[1]] * tk_width
	tk_xMax <- plt[[2]] * tk_width
	tk_yOrig <- plt[[3]] * tk_height
	tk_yMax <- plt[[4]] * tk_height
	
	tkRxratio <- (tk_xMax - tk_xOrig) /	(usr[[2]] - usr[[1]])
	tkRyratio <- (tk_yMax - tk_yOrig) / (usr[[4]] - usr[[3]])

	tk_xMouse <- tk_xMouse - tk_xOrig
	tk_yMouse <- tk_yMouse - tk_yOrig
	
	R_xMouse <- tk_xMouse / tkRxratio + usr[[1]]
	R_yMouse <- tk_yMouse / tkRyratio + usr[[3]]
	
	c(R_xMouse, R_yMouse)
})


setGeneric("onPlotClicked", function(x, y) standardGeneric("onPlotClicked"))
setMethod("onPlotClicked", c("character", "character"), function(x, y) {
	plotParam <- getVar("plotParam")
	curMkr <- getVar("curMkr")
	coords <- toRCoords(c(x, y))
	usr <- plotParam$usr
	map <- regEval("coll[[curSet, curMap]]")
  
  if(length(map) == 0) {
    tkmessageBox(title = "Error", message = "   Please select one map on a set   ")
  } else {
    sqdist <- apply(as(map, "data.frame"), 1, function(mkr) {
	        xDot <- as.numeric(mkr[["phys"]]) / 1000000
	        yDot <- as.numeric(mkr[["gen"]])
	        ((coords[[1]] - xDot) / (usr[[2]] - usr[[1]])) ^ 2 + ((coords[[2]] - yDot) / (usr[[4]] - usr[[3]])) ^ 2
      })
    if(which.min(sqdist) == curMkr)
	      setVar("curMkr", 0)
    else
	      setVar("curMkr", which.min(sqdist))  
  }
})


setGeneric("onPlotScrollUp", function(x, y) standardGeneric("onPlotScrollUp"))
setMethod("onPlotScrollUp", c("character", "character"), function(x, y) {
	plotParam <- getVar("plotParam")
	zoomRatio <- getVar("zoomRatio")
	mapBounds <- getVar("mapBounds")
	usr <- plotParam$usr
	plt <- plotParam$plt
	
	coords <- toRCoords(c(x, y))
  map <- regEval("coll[[curSet, curMap]]")
  
  if(length(map) == 0) {
    tkmessageBox(title = "Error", message = "   Please select one map on a set   ")
  } else {
  	xMouse <- coords[[1]]
  	yMouse <- coords[[2]]
  	
  	xMousePosRatio <- (xMouse - usr[[1]]) / (usr[[2]] - usr[[1]])
  	
  	xDec <- (usr[[2]] - usr[[1]]) * zoomRatio
  	usr[[1]] <- usr[[1]] + xDec * xMousePosRatio
  	usr[[2]] <- usr[[2]] - xDec * (1 - xMousePosRatio)				
  	
  	yMousePosRatio <- (yMouse - usr[[3]]) / (usr[[4]] - usr[[3]])
  	
  	yDec <- (usr[[4]] - usr[[3]]) * zoomRatio
  	usr[[3]] <- usr[[3]] + yDec * yMousePosRatio
  	usr[[4]] <- usr[[4]] - yDec * (1 - yMousePosRatio)	
  		
  	if(usr[[1]] < mapBounds[[1]])
  		usr[[1]] <- mapBounds[[1]]
  	if(usr[[2]] > mapBounds[[2]])
  		usr[[2]] <- mapBounds[[2]]
  	if(usr[[3]] < mapBounds[[3]])
  		usr[[3]] <- mapBounds[[3]]
  	if(usr[[4]] > mapBounds[[4]])
  		usr[[4]] <- mapBounds[[4]]
  	
  	setVar("pltLim", usr)
  }
})


setGeneric("onPlotScrollDown", function(x, y) standardGeneric("onPlotScrollDown"))
setMethod("onPlotScrollDown", c("character", "character"), function(x, y) {
	plotParam <- getVar("plotParam")
	zoomRatio <- getVar("zoomRatio")
	mapBounds <- getVar("mapBounds")
	usr <- plotParam$usr
	plt <- plotParam$plt

	coords <- toRCoords(c(x, y))
  map <- regEval("coll[[curSet, curMap]]")
  
  if(length(map) == 0) {
    tkmessageBox(title = "Error", message = "   Please select one map on a set   ")
  } else {
    xMouse <- coords[[1]]
  	yMouse <- coords[[2]]
  	
  	xMousePosRatio <- (xMouse - usr[[1]]) / (usr[[2]] - usr[[1]])
  	
  	xDec <- -(usr[[2]] - usr[[1]]) * zoomRatio
  	usr[[1]] <- usr[[1]] + xDec * xMousePosRatio
  	usr[[2]] <- usr[[2]] - xDec * (1 - xMousePosRatio)				
  	
  	yMousePosRatio <- (yMouse - usr[[3]]) / (usr[[4]] - usr[[3]])
  
  	yDec <- -(usr[[4]] - usr[[3]]) * zoomRatio
  	usr[[3]] <- usr[[3]] + yDec * yMousePosRatio
  	usr[[4]] <- usr[[4]] - yDec * (1 - yMousePosRatio)
  
  	if(usr[[1]] < mapBounds[[1]])
  		usr[[1]] <- mapBounds[[1]]
  	if(usr[[2]] > mapBounds[[2]])
  		usr[[2]] <- mapBounds[[2]]
  	if(usr[[3]] < mapBounds[[3]])
  		usr[[3]] <- mapBounds[[3]]
  	if(usr[[4]] > mapBounds[[4]])
  		usr[[4]] <- mapBounds[[4]]
  		
  	setVar("pltLim", usr)
  }
})


##--- draw in the top area plot
setGeneric("drawPlt1", function(dummy) standardGeneric("drawPlt1"))
setMethod("drawPlt1", "missing", function(dummy) {
	curMkr <- getVar("curMkr")
	pltLim <- getVar("pltLim")
	map <- regEval("coll[[curSet,curMap]]")
	
	if(is.null(map)) {
		par(bg = "white", col = "white", col.lab = "white", col.axis = "white", tck = 0)
		plot(0, 0)
	} else {
    par(bg = "white")
		cl <- match.call()
		if(!is.na(pltLim[1])) {
			cl$xlim <- pltLim[1:2]
			cl$ylim <- pltLim[3:4]
		} else {
			cl$ylim <- c(min(geneticDistances(map), na.rm = T), max(geneticDistances(map), na.rm = T))
			cl$xlim <- c(min(physicalPositions(map) / 1000000, na.rm = T), max(physicalPositions(map) / 1000000, na.rm = T))
		}
		cl$object <- map
		cl[[1]] <- as.name("plotModels")
		eval(cl)
		par(new = T)
		cl[[1]] <- as.name("plotMarkers")
		eval(cl)
		if(curMkr > 0) {
			mkr <- map[[curMkr]]
			points(as.numeric(mkr[[2]]) / 1000000, as.numeric(mkr[[3]]), col = "#FF0000", pch = 19, type = "o")
		}
		lapply(getVar("curEst"), function(est) {
			abline(v = est / 1000000, col = "red")
		})
		setVar("plotParam", par())
	}
})


##--- draw in the bottom area plot
setGeneric("drawPlt2", function(dummy) standardGeneric("drawPlt2"))
setMethod("drawPlt2", "missing", function(dummy) {
	curMkr <- getVar("curMkr")
	pltLim <- getVar("pltLim")
	map <- regEval("coll[[curSet,curMap]]")
	
	if(is.null(map)) {
		par(bg = "white", col = "white", col.lab = "white", col.axis = "white", tck = 0)
		plot(0, 0)
	} else {
    par(bg = "white") 
		cl <- match.call()
		if(!is.na(pltLim[1])) {
			cl$xlim <- pltLim[1:2]
		} else {
			cl$xlim <- c(min(physicalPositions(map) / 1000000, na.rm = T), max(physicalPositions(map) / 1000000, na.rm = T))
		}
		cl$object <- map
		cl[[1]] <- as.name("plotRates")
		eval(cl)
    if(length(interpolations(map)) > 0)
			abline(h = 0, col = "red", lty = 2)
		lapply(getVar("curEst"), function(est) {
			abline(v = est / 1000000, col = "red")
		})
	}
})


setGeneric("plotMap", function(dummy) standardGeneric("plotMap"))
setMethod("plotMap", "missing", function(dummy) {
	curMkr <- getVar("curMkr")
	pltLim <- getVar("pltLim")
	map <- regEval("coll[[curSet, curMap]]")
	
	if(is.null(map)) {
		par(bg = "white", col = "white", col.lab = "white", col.axis = "white", tck = 0)
		plot(0, 0)
	} else {
		par(mar = c(4, 5, 1, 1) + 0.1)
		cl <- match.call()
		if(!is.na(pltLim[1])) {
			cl$xlim <- pltLim[1:2]
			cl$ylim <- pltLim[3:4]
		} else {
			cl$ylim <- c(0, max(geneticDistances(map), na.rm = T))
			cl$xlim <- c(0, max(physicalPositions(map) / 1000000, na.rm = T))
		}
		
		cl$xaxt <- "n"
		cl$yaxt <- "n"
		cl$xlab <- ""
		cl$ylab <- ""
		cl[[1]] <- as.name("plotModels")
		lapply(interpolations(map), function(itr) {
			if(visible(itr)) {
				cl$object <- itr
				eval(cl)
				par(new = T)
			}
		})
		
		cl$xaxt <- NULL
		cl$yaxt <- NULL
		cl$xlab <- NULL
		cl$ylab <- NULL
		cl$object <- NULL
		cl$x <- map
		cl[[1]] <- as.name("plot")
		eval(cl)
	}

	if(curMkr > 0) {
		mkr <- map[[curMkr]]
		points(as.numeric(mkr[[2]]) / 1000000, as.numeric(mkr[[3]]), col = "#FF0000", pch = 19, type = "o")
	}
	setVar("plotParam", par())
})


setGeneric("pltRates", function(dummy) standardGeneric("pltRates"))
setMethod("pltRates", "missing", function(dummy) {
	par(mar = c(4, 5, 0, 1) + 0.1)
	curMkr <- getVar("curMkr")
	pltLim <- getVar("pltLim")

	map <- regEval("coll[[curSet,curMap]]")
	
	if(is.null(map) || length(interpolations(map)) < 1) {
		par(bg = "white", col = "white", col.lab = "white", col.axis = "white", tck = 0)
		plot(0, 0)
	} else {
		cl <- match.call()
		if(!is.na(pltLim[1]))
			cl$xlim <- pltLim[1:2]
		else
			cl$xlim <- c(0, max(physicalPositions(map) / 1000000, na.rm = T))			
		
    cl$ylim <- c(0, max(unlist(lapply(map@interpolations, function(itr) {max(rates(itr), na.rm = T)}))) + 0.5)
    
		## setting up axis
		cl[[1]] <- as.name("plot")
		cl$object <- NULL
		cl$xlab <- "Physical positions (Mb)"
		cl$ylab <- "Recomb. rates\n (cM/Mb)"
		cl$x <- 0
		cl$y <- 0
		cl$col <- "white"
		eval(cl)
		cl$xaxt <- "n"
		cl$yaxt <- "n"
		cl$ylab <- ""
		cl$xlab <- ""
		cl[[1]] <- as.name("plotRates")
		lapply(interpolations(map), function(itr) { 
			if(visible(itr)) {
				par(new = T)
				cl$object <- itr
				eval(cl)
			}
		})
	}
})


setGeneric("PlotArea", function(parent) standardGeneric("PlotArea"))
setMethod("PlotArea", "ANY", function(parent) {
	if(class(parent) != "tkwin")
		return()
	regVar("plotParam")
	regVar("xLim", NA)
	regVar("yLim", NA)
	regVar("pltLim", NA)
	regVar("zoomRatio", 0.1)
	regVar("mapBounds", c(-10, 10, -10, 10))
	regVar("savBounds", T)
	
	frm <- tkframe(parent, relief = "sunken", borderwidth = 2)
	plt1 <- tkrplot(frm, fun = drawPlt1, hscale = 1, vscale = 0.8)
	plt2 <- tkrplot(frm, fun = drawPlt2, hscale = 1, vscale = 0.4)
	
	regVar("plotArea", plt1)
	tkgrid(plt1)
	tkgrid(plt2)
	tkbind(plt1, "<ButtonRelease-5>", onPlotScrollDown)
	tkbind(plt1, "<ButtonRelease-4>", onPlotScrollUp)
	tkbind(plt1, "<ButtonRelease-1>", onPlotClicked)
	
	attachCallback("tkrreplot", function() {
		tkrreplot(plt1)
		tkrreplot(plt2)		
		if(getVar("savBounds")) {
			setVar("mapBounds", getVar("plotParam")$usr)
			setVar("savBounds", F)
		}
	}, "curMkr")

	attachCallback("tkrreplot", function() {
		if(!is.na(getVar("pltLim")[1])) {
			tkrreplot(plt1)
			tkrreplot(plt2)
		}
	},"pltLim")

	attachCallback("updateLimsAndRecMapBounds", function() {
		setVar("pltLim", NA)
		setVar("savBounds", T)
	}, "curMap")

	attachCallback("drawest", function() {
		tkrreplot(plt2)
		tkrreplot(plt1)
	}, "curEst")
	frm
})


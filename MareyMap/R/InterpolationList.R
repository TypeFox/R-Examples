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


setGeneric("onChangeColor",	function(name, tag) standardGeneric("onChangeColor"))
setMethod("onChangeColor", "character", function(name, tag) {
	coll <- getVar("coll")
	curSet <- getVar("curSet")
	curMap <- getVar("curMap")
	
	itr <- interpolation(coll[[curSet, curMap]], name)
	col <- tclvalue(tcl("tk_chooseColor", initialcolor = color(itr), title = "color"))
	if(nchar(col) != 0) {
		color(itr) <- col
		coll[[curSet, curMap]] <- coll[[curSet, curMap]] + itr
	}
	tkconfigure(tag, bg = color(itr))
	setVar("coll", coll)
	## force map redrawing 
	setVar("curMkr", getVar("curMkr"))
})


setGeneric("onToggleVisibility", function(name) standardGeneric("onToggleVisibility"))
setMethod("onToggleVisibility", "character", function(name) {
	coll <- getVar("coll")
	curSet <- getVar("curSet")
	curMap <- getVar("curMap")
	if(visible(interpolation(coll[[curSet, curMap]], name)))
		visible(interpolation(coll[[curSet, curMap]], name)) <- F
	else
		visible(interpolation(coll[[curSet, curMap]], name)) <- T
	setVar("coll", coll)
	## force map redrawing 
	setVar("curMkr", getVar("curMkr"))
})


setGeneric("onTogglePersistency",	function(name) standardGeneric("onTogglePersistency"))
setMethod("onTogglePersistency", "character", function(name) {
	coll <- getVar("coll")
	curSet <- getVar("curSet")
	curMap <- getVar("curMap")
	if(persistent(interpolation(coll[[curSet, curMap]], name)))
		persistent(interpolation(coll[[curSet, curMap]], name)) <- F
	else
		persistent(interpolation(coll[[curSet, curMap]], name)) <- T
	setVar("coll", coll)
  setVar("curSet", curSet)
  setVar("curMap", curMap)
})


setGeneric("buildInterList", function(frm) standardGeneric("buildInterList"))
setMethod("buildInterList", "ANY", function(frm) {
  frm <- getVar("interCenterFrm")
	.Tcl(paste("eval destroy [winfo children ", as.character(frm), "]"))
	map <- regEval("coll[[curSet, curMap]]")
	if(!is.null(map)) {
		lapply(interpolations(map), function(x) {
			frm <- tkframe(frm)
			coltag <- tkcanvas(frm, width = "10", height = "10", bg = color(x))
			tkbind(coltag, "<ButtonRelease-1>", function() onChangeColor(name(x), coltag))
			if(visible(x))
				visVar <- tclVar("1")
			else
				visVar <- tclVar("0")
			visTick <- tkcheckbutton(frm,	variable = visVar) #, font=tkfont.create(size=10)
			tkbind(visTick, "<ButtonRelease-1>", function() onToggleVisibility(name(x)))
			if(persistent(x))
				savVar <- tclVar("1")
			else
				savVar <- tclVar("0")
			savTick <- tkcheckbutton(frm, variable = savVar) #, font=tkfont.create(size=10))
			tkbind(savTick, "<ButtonRelease-1>", function() onTogglePersistency(name(x)))
			namlbl <- tklabel(frm, text = name(x))
			tkpack(frm, side = "top", anchor = "w")
			tkpack(tklabel(frm, text = " "), coltag, tklabel(frm, text = " "), visTick, savTick, tklabel(frm, text = " "), namlbl, side = "left")
		})
	}
})


setGeneric("interParamGui",	function(inter) standardGeneric("interParamGui"))
setMethod("interParamGui", "Interpolation", function(inter) {
    
  #--- interParamGui callbacks
  onParamHlp <- function(x) {
	  hlpwin <- tktoplevel()
	  tkgrab.set(hlpwin)
	  tkwm.protocol(hlpwin, "WM_DELETE_WINDOW", function() {
		  tkgrab.set(dialog)
		  tkdestroy(hlpwin)
	  })
	  tkwm.title(hlpwin, paramName(x))
	  tkgrid(tklabel(hlpwin, text = ""))
	  desc <- paramDesc(x)
	  tkgrid(tklabel(hlpwin, text = desc), sticky = "ew")
	  but <- tkbutton(hlpwin, text = "Discard",	command = function() {
		 tkgrab.set(dialog)
		 tkdestroy(hlpwin)
		})
	  tkgrid(but, padx = 10, pady = 10)
  }
	#---  end of interParamGui callbacks
    
	coll <- getVar("coll")
	curSet <- getVar("curSet")
	curMap <- getVar("curMap")
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	oldname <- name(inter)
	tkwm.title(dialog, "Parameters")
	tkgrid(tklabel(dialog, text = " "))
	bottom <- tkframe(dialog)
	finbut <- tkbutton(bottom, text = "Finish", command = function() {
			tkdestroy(dialog)
			if(getVar("calcall") == T || getVar("calcset") == T) {
        progdial <- tktoplevel()
        tkwm.title(progdial, "Please wait...")
        tkpack(tklabel(progdial, text = ""))
        proglbl <- tklabel(progdial, text = "", font = "courier")
        tkpack(tklabel(progdial), proglbl, tklabel(progdial))
        tkpack(tklabel(progdial, text = ""))
        counter <- 0

				if(getVar("calcall") == T) {
					total <- length(coll)
					sets <- coll@sets
				}
				if(getVar("calcset") == T) {
					total <- length(coll[[curSet]])
					sets <- list(coll[[curSet]])
				}
        
        tkgrab.set(progdial)
				lapply(sets, function(set) {lapply(set@maps, function(map) {
						interpolations(map)[[oldname]] <- NULL
						map <- map + inter
						coll <<- coll + map
            counter <<- counter + 1
            prog <- ""

            if(counter %% 4 == 0)
              prog <- "|"
            else if (counter %% 4 == 1)
              prog <- "/"
            else if (counter %% 4 == 2)
              prog <- "-"
            else
							prog <- "\\"
            prog <- paste(prog, "[", sep = "")
            blocks <- round(counter / total * 50)

            for(i in 1:blocks) {
              prog <- paste(prog, "=", sep = "")
            }
            for(i in 1:(50 - blocks)) {
              prog <- paste(prog, " ", sep = "")
            }
            
            prog <- paste(prog, "]", sep = "")
            tkconfigure(proglbl, text = prog)
            tcl("update")
				})})
  			tkdestroy(progdial)
        
			} else {
				interpolations(coll[[curSet, curMap]])[[oldname]] <- NULL
				coll[[curSet, curMap]] <- coll[[curSet, curMap]] + inter
			}
			setVar("coll", coll)
			## force UI update
			setVar("calcall", F)
			setVar("calcset", F)
			setVar("curMap", getVar("curMap"))
			setVar("curMkr", getVar("curMkr"))
			buildInterList()
	})
  
	bakbut <- tkbutton(bottom, text = "Cancel", command = function() {
		setVar("calcset", F)
		setVar("calcall", F)
		tkdestroy(dialog)
	})

	tkpack(finbut, side = "right")
	tkpack(bakbut, side = "left")
	env <- environment()
	lapply(userParam(inter), function(x) {
		frm <- tkframe(dialog)
		lbl <- tklabel(frm, text = paste(paramName(x), " : "))
		tcl("image", "create", "photo", "stock_help-agent",	file = system.file("stock_help-agent.gif", package = "MareyMap"))
		hlp <- tkbutton(frm, image = "stock_help-agent", command = function() onParamHlp(x))
		entry <- 0
		type <- paramType(x)
    
		if(type == "logical") {
			updatefun <- function() {
				txtval <- "F"
				if(tclvalue(val) == "1")
					txtval <- "T"
				expr <- paste(paramFun(x), "(inter) <- ", txtval, sep = "")
				eval(parse(text = expr), envir = env)
			}
			startval <- eval(parse(text = paste(paramFun(x), "(inter)", sep = "")))
			val <- tclVar(as.logical(startval))
			entry <- tkcheckbutton(frm, variable = val)
			tkbind(entry, "<ButtonRelease-1>", updatefun)
      
		}	else if(type == "character") {
			updatefun <- function() {
				if(is.na(text <- as.character(tclvalue(txtval))))
					text <- paramDefault(x)
				else if((!is.null(paramMin(x))) && nchar(text) < paramMin(x)) {
					text <- paramDefault(x)
				} else if((!is.null(paramMax(x)))	&& nchar(text) > paramMax(x)) {
					text <- strtrim(text, paramMax(x))
				}
				tclvalue(txtval) <- as.character(text)
				expr <- paste(paramFun(x), "(inter)<-\"", text, "\"", sep = "")
				eval(parse(text = expr), envir = env)
			}
			startval <- eval(parse(text =	paste(paramFun(x), "(inter)", sep = "")))
			txtval <- tclVar(as.character(startval))
			w <- paramMax(x)
			if(is.null(w) || (w > 30))
				w <- 30
			entry <- tkentry(frm, textvariable = txtval, width = w)
			tkbind(entry, "<FocusOut>", updatefun)
			.Tcl(paste("bind ", .Tk.ID(finbut), " <ButtonRelease-1> ", "{+ ", .Tcl.callback(updatefun), " }", sep = ""))
      
		}	else if(type == "numeric") {
			updatefun <- function() {
				if(is.na(numval <- as.numeric(tclvalue(txtval))))
					numval <- paramDefault(x)
				else if((!is.null(paramMin(x))) && numval < paramMin(x)) {
					numval <- paramMin(x)
				} else if((!is.null(paramMax(x))) && numval > paramMax(x)) {
					numval <- paramMax(x)
				}
				tclvalue(txtval) <- as.character(numval)
				expr <- paste(paramFun(x), "(inter)<-", numval, sep = "")
				eval(parse(text = expr), envir = env)
			}
			startval <- eval(parse(text = paste(paramFun(x), "(inter)", sep = "")))
			txtval <- tclVar(as.character(startval))
			entry <- tkentry(frm, width = "10", textvariable = txtval)
			tkbind(entry, "<FocusOut>", updatefun)
			.Tcl(paste("bind ", .Tk.ID(finbut), " <ButtonRelease-1> ", "{+ ", .Tcl.callback(updatefun), " }", sep = ""))
		
    }	else if(type == "integer") {
			updatefun <- function() {
				if(is.na(intval <- as.integer(tclvalue(txtval))))
					intval <- paramDefault(x)
				else if((!is.null(paramMin(x)))	&& intval < paramMin(x)) {
					intval <- paramMin(x)
				} else if((!is.null(paramMax(x))) && intval > paramMax(x)) {
					intval <- paramMax(x)
				}
				tclvalue(txtval) <- as.character(intval)
				intval <- as.integer(intval)
				expr <- paste(paramFun(x), "(inter)<-", intval, sep = "")
				eval(parse(text = expr), envir = env)
			}		
			startval <- eval(parse(text = paste(paramFun(x), "(inter)", sep = "")))
			txtval <- tclVar(as.character(startval))
			entry <- tkentry(frm, width = "10", textvariable = txtval)
			tkbind(entry, "<FocusOut>", updatefun)
			.Tcl(paste("bind ", .Tk.ID(finbut), " <ButtonRelease-1> ", "{+ ", .Tcl.callback(updatefun), " }", sep = ""))
		
    }	else if(type == "color") {			
			updatefun <- function() {
				expr <- paste(paramFun(x), "(inter)", sep = "")
				col <- eval(parse(text = expr), envir = env)
				col <- tclvalue(tcl("tk_chooseColor",	initialcolor = col,	title = paramName(x)))
				if(nchar(col) != 0) {
					tkconfigure(entry, bg = col)
					expr <- paste(paramFun(x), "(inter)<-\"",	col, "\"", sep = "")
					eval(parse(text = expr), envir = env)
				}
			}
			startval <- eval(parse(text =	paste(paramFun(x), "(inter)", sep = "")))
			entry <- tklabel(frm, text = " ",	bg = startval, width = 3)
			tkbind(entry, "<ButtonRelease-1>", updatefun)
		
    } else if(type == "list") {
      updatefun <- function() {
	      val <- paramValues(x)[[as.integer(tkcurselection(entry)) + 1]]
  				if(is.numeric(val))
          	expr <- paste(paramFun(x), "(inter)<-", val, sep = "")
          else
            expr <- paste(paramFun(x), "(inter)<-\"", val, "\"", sep = "")
				eval(parse(text = expr), envir = env)
			}
      entry <- tklistbox(frm, selectmode = "single")
      lapply(paramValues(x), function(pp) {
      	tkinsert(entry, "end", pp)
        if(pp == paramDefault(x))
          tkselection.set(entry, "end")
      })
      tkbind(entry, "<ButtonRelease-1>", function() {
        if(tclvalue(tksize(entry)) > 0)
        	updatefun()
      })
  	}
  
		tkpack(lbl, entry, side = "left")
		tkpack(hlp, tklabel(frm, text = "\t"), side = "right")
		tkgrid(frm, padx = 10, pady = 5, sticky = "ew")
	})
  
	tkgrid(bottom, padx = 10, pady = 10, sticky = "ew")
})


##--------------- To add a new interpolation ---------------------------------------------##
setGeneric("onAddInter", function(dummy) standardGeneric("onAddInter"))
setMethod("onAddInter", "character", function(dummy) {

	## onAddInter callbacks
	## called when the 'select' button is clicked
	onSelectInterpolationMethod <- function() {
		itrlst <- getInterList()
		itrname <- names(itrlst)[[as.integer(tkcurselection(lst)) + 1]]
		itrclass <- itrlst[[as.integer(tkcurselection(lst)) + 1]]
		inter <- new(itrclass)
		i <- 1
		basename <- name(inter)
		while(name(inter) %in% lapply(regEval("interpolations(coll[[curSet, curMap]])"), function(x) name(x))) {
	    i <- i + 1
      name(inter) <- paste(basename, "-", i, sep = "")
		}
		color(inter) <- substr(rainbow(100, start = 0.1, end = 0.9)[runif(1, 2, 100)[1]], 1, 7)
		tkdestroy(dialog)
		interParamGui(inter)
	}
  
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Methods")
	lbl <- tklabel(dialog, text = "   Choose an interpolation method:   ")
	lst <- tklistbox(dialog, selectmode = "single")
	lapply(names(getInterList()),	function(x) {tkinsert(lst, "end", x)})
	tkbind(lst, "<ButtonRelease-1>", function() {
		if(tclvalue(tksize(lst)) > 0)
			tkconfigure(selbut, state = "normal")
	})
	valcalcset <- tclVar("0")
	valcalcall <- tclVar("0")
	setVar("calcset", F)
	setVar("calcall", F)
	
	calcsetbut <- tkcheckbutton(dialog, text = "Apply to every map in the set", command = function() {
		if(tclvalue(valcalcset) != 0) {
		  setVar("calcset", T)
			setVar("calcall", F)
			tclvalue(valcalcall) <- 0
    } else
      setVar("calcset", F)
  })
	tkconfigure(calcsetbut, variable = valcalcset)
	
#  calcallbut <- tkcheckbutton(dialog, text = "Apply to every map in the set", command = function() {
#    if(tclvalue(valcalcall) != 0) {
#      setVar("calcall", T)
#      setVar("calcset", F)
#      tclvalue(valcalcset) <- 0
#    } else
#      setVar("calcall", F)
#  })
#  tkconfigure(calcallbut, variable = valcalcall)

	bottom <- tkframe(dialog)
	canbut <- tkbutton(bottom, text = "Cancel", command = function() {tkdestroy(dialog)})
	selbut <- tkbutton(bottom, text = "Select >>", state = "disabled", command = onSelectInterpolationMethod)
	
	tkpack(canbut, side = "left")
	tkpack(selbut, side = "right")	
	tkgrid(lbl, padx = 10, pady = 10)
	tkgrid(lst, padx = 10)
	tkgrid(calcsetbut, padx = 10, pady = 5)
#  tkgrid(calcallbut, padx = 10, pady = 5)
	tkgrid(bottom, padx = 10, pady = 10)
	tkgrid.configure(lbl, sticky = "ew")
	tkgrid.configure(lst, sticky = "nsew")
	tkgrid.configure(bottom, sticky = "ew")
})


##--------------- To remove a selected interpolation ---------------------------------------------##
setGeneric("onRemInter", function(dummy) standardGeneric("onRemInter"))
setMethod("onRemInter", "character", function(dummy) {
	coll <- getVar("coll")
	curSet <- getVar("curSet")
	curMap <- getVar("curMap")
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Remove an interpolation")
	lbl <- tklabel(dialog, text = "Choose the interpolation to be removed: ")
	lst <- tklistbox(dialog, selectmode = "single")
	tkbind(lst, "<ButtonRelease-1>", function() {
		if(tclvalue(tksize(lst)) > 0)
			tkconfigure(rembut, state = "normal")
	})
	lapply(interpolations(coll[[curSet, curMap]]), function(x) {
		tkinsert(lst, "end", name(x))
	})

  valcalcset <- tclVar("0")
	valcalcall <- tclVar("0")
	setVar("calcset", F)
	setVar("calcall", F)
	
	calcsetbut <- tkcheckbutton(dialog, text = "Apply to every map in the set", command = function() {
    if(tclvalue(valcalcset) != 0) {
    	setVar("calcset", T)
			setVar("calcall", F)
			tclvalue(valcalcall) <- 0
     } else
      setVar("calcset", F)
  })
	tkconfigure(calcsetbut, variable = valcalcset)
	
#	calcallbut <- tkcheckbutton(dialog, text = "Apply to every map in the set", command = function() {
#    if(tclvalue(valcalcall) != 0) {
#    	setVar("calcall", T)
#		  setVar("calcset", F)
#			tclvalue(valcalcset) <- 0
#     } else
#      setVar("calcall", F)
#  })
#	tkconfigure(calcallbut, variable = valcalcall)
  
	bottom <- tkframe(dialog)
	canbut <- tkbutton(bottom, text = "Cancel", command = function() {
		setVar("calcall", F)
		setVar("calset", F)
		tkdestroy(dialog)
	})
	rembut <- tkbutton(bottom, text = "Remove", state = "disabled",	command = function() {
  	coll <- getVar("coll")
  	idx <- as.integer(tkcurselection(lst)) + 1
  	if(getVar("calcall") == T) {
  		itrname <- names(interpolations(coll[[curSet, curMap]]))[[idx]]
  		for(set in setNames(coll)) {
  			for(map in mapNames(coll[[set]])) {
  				interpolations(coll[[set, map]])[[itrname]] <- NULL
  		}}
  		setVar("coll", coll)
  	} else if(getVar("calcset") == T) {
  		itrname <- names(interpolations(coll[[curSet, curMap]]))[[idx]]
  		for(map in mapNames(coll[[curSet]])) {
  			interpolations(coll[[curSet, map]])[[itrname]] <- NULL
  		}
  		setVar("coll", coll)
  	} else {
  		interpolations(coll[[curSet, curMap]])[[idx]] <- NULL
  		setVar("coll", coll)
  	}
  	setVar("calcset", F)
  	setVar("calcall", F)
  	#force UI update
  	setVar("curMap", getVar("curMap"))
  	setVar("curMkr", getVar("curMkr"))
  	tkdestroy(dialog)
	})

	tkpack(canbut, side = "left")
	tkpack(rembut, side = "right")
	tkpack(tklabel(dialog, text = ""), side = "top")
	tkpack(lbl, padx = 10, side = "top")
	tkpack(lst, pady = 10, side = "top")
	tkpack(calcsetbut, padx = 10, pady = 5)
#	tkpack(calcallbut, padx = 10, pady = 5)
	tkpack(bottom, padx = 10, side = "top", fill = "x")
	tkpack(tklabel(dialog, text = ""), side = "top")
})


##--------------- To manage interpolation configuration ---------------------------------------------##
setGeneric("onCfgInter",function(dummy) standardGeneric("onCfgInter"))
setMethod("onCfgInter", "character", function(dummy) {
	dialog <- tktoplevel()
	tkgrab.set(dialog)
	tkwm.title(dialog, "Configure an interpolation")
	lbl <- tklabel(dialog, text = "Choose the interpolation to be configured: ")
	lst <- tklistbox(dialog, selectmode = "single")
	tkbind(lst, "<ButtonRelease-1>", function() {
		if(tclvalue(tksize(lst)) > 0)
			tkconfigure(selbut, state = "normal")
	})
	lapply(regEval("interpolations(coll[[curSet, curMap]])"), function(x) {tkinsert(lst, "end", name(x))})
  valcalcset <- tclVar("0")
	valcalcall <- tclVar("0")
	setVar("calcset", F)
	setVar("calcall", F)
	
	calcsetbut <- tkcheckbutton(dialog, text = "Apply to every map in the set", command = function() {
    if(tclvalue(valcalcset) != 0) {
    	setVar("calcset", T)
      setVar("calcall", F)
			tclvalue(valcalcall) <- 0
     } else
       setVar("calcset", F)
  })
	tkconfigure(calcsetbut, variable = valcalcset)
	
#	calcallbut <- tkcheckbutton(dialog, text = "Apply to every map in the set", command = function() {
#    if(tclvalue(valcalcall) != 0) {
#    	setVar("calcall", T)
#			setVar("calcset", F)
#			tclvalue(valcalcset) <- 0
#    } else
#      setVar("calcall", F)
#  })
#	tkconfigure(calcallbut, variable = valcalcall)

	bottom <- tkframe(dialog)
	canbut <- tkbutton(bottom, text = "Cancel",	command = function() {tkdestroy(dialog)})
	selbut <- tkbutton(bottom, text = "Select", state = "disabled", command = function() {
		idx <- as.integer(tkcurselection(lst)) + 1
		interParamGui(regEval("interpolations(coll[[curSet, curMap]])")[[idx]])			
		tkdestroy(dialog)
	})
	tkpack(canbut, side = "left")
	tkpack(selbut, side = "right")
	tkpack(tklabel(dialog, text = ""), side = "top")
	tkpack(lbl, padx = 10, side = "top")
	tkpack(lst, pady = 10, side = "top")
  tkpack(calcsetbut, padx = 10, side = "top")
#	tkpack(calcallbut, padx = 10, side = "top")
	tkpack(bottom, padx = 10, side = "top", fill = "x")
	tkpack(tklabel(dialog, text = ""), side = "top")
})


setGeneric("InterpolationList",	function(parent) standardGeneric("InterpolationList"))
setMethod("InterpolationList", "ANY", function(parent) {
	regVar("calcall", F)
	regVar("calcset", F)
	interFrm <- tkframe(parent, relief = "sunken", borderwidth = 2)
  
  ##------------- Top Frame in Interpolation
	interTopFrm <- tkframe(interFrm, relief = "groove", borderwidth = 2)
	tkpack(interTopFrm, side = "top", fill = "x")
  
	## colors icon
	tcl("image", "create", "photo", "stock_color", file = system.file("stock_color.gif", package = "MareyMap"))
	interColorsIcn <- tklabel(interTopFrm, image = "stock_color", state = "disabled")
  
	## visible icon
	tcl("image", "create", "photo", "stock_show-all", file = system.file("stock_show-all.gif", package = "MareyMap"))		
	interVisibleIcn <- tklabel(interTopFrm, image = "stock_show-all", state = "disabled")
		
  ## save icon
	tcl("image", "create", "photo", "stock_save", file = system.file("stock_save.gif", package = "MareyMap"))
	interSaveIcn <- tklabel(interTopFrm, image = "stock_save", state = "disabled")
	
  ## name label
	interNameLbl <- tklabel(interTopFrm, text = "name", state = "disabled")
	
	interlbl <- tklabel(interTopFrm, text = "INTERPOLATIONS", width = "25")
	tkpack(interlbl, side = "top")
	tkpack(tklabel(interTopFrm, text = ""), interColorsIcn,
				tklabel(interTopFrm, text = " "),	interVisibleIcn,
				tklabel(interTopFrm, text = " "),	interSaveIcn,
				tklabel(interTopFrm, text = "  "), interNameLbl, side = "left")

  ##------------- Center Frame in Interpolation    
	interCenterFrm <- tkframe(interFrm)
	tkpack(interCenterFrm, side = "top", fill = "both", anchor = "nw")
	regVar("interCenterFrm", interCenterFrm)
  
  ##------------- Bottom Frame in Interpolation
	interBottomFrm <- tkframe(interFrm)
	tkpack(interBottomFrm, side = "bottom", fill = "x")
  
	## add interpolation button
	tcl("image", "create", "photo", "gtk-add", file = system.file("gtk-add.gif", package = "MareyMap"))
	addItrBut <- tkbutton(interBottomFrm, image = "gtk-add", command = onAddInter, state = "disabled")
  
	## remove interpolation button
	tcl("image", "create", "photo", "stock_delete", file = system.file("stock_delete.gif", package = "MareyMap"))
	remItrBut <- tkbutton(interBottomFrm, image = "stock_delete", command = onRemInter, state = "disabled")
  
	## config interpolation button
	tcl("image", "create", "photo", "gnome-settings", file = system.file("gnome-settings.gif", package = "MareyMap"))
	cfgItrBut <- tkbutton(interBottomFrm, image = "gnome-settings", command = onCfgInter, state = "disabled")
	
	tkgrid(addItrBut, remItrBut, cfgItrBut, sticky = "ew")
	
	attachCallback("updateInterList", function() {
		curMap <- getVar("curMap")
		curSet <- getVar("curSet")
		coll <- getVar("coll")
		if(curMap != "")
			tkconfigure(addItrBut, state = "normal")
		else
			tkconfigure(addItrBut, state = "disabled")
		map <- regEval("coll[[curSet, curMap]]")
		if((!is.null(map)) && length(interpolations(map))) {
			tkconfigure(remItrBut, state = "normal")
			tkconfigure(cfgItrBut, state = "normal")
      tkconfigure(interColorsIcn, state = "normal")
      tkconfigure(interVisibleIcn, state = "normal")
      tkconfigure(interSaveIcn, state = "normal")
      tkconfigure(interNameLbl, state = "normal")
		} else {
			tkconfigure(remItrBut, state = "disabled")
			tkconfigure(cfgItrBut, state = "disabled")
      tkconfigure(interColorsIcn, state = "disabled")
      tkconfigure(interVisibleIcn, state = "disabled")
      tkconfigure(interSaveIcn, state = "disabled")
      tkconfigure(interNameLbl, state = "disabled")
		}
		buildInterList(interCenterFrm)
	}
	,"curMap")

	interFrm
})


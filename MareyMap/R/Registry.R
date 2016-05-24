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


reg <- new.env()
assign(".hooks", list(), envir = reg)
assign("inter.list", list(), envir = reg)


## assign a value to a variable
setGeneric("regVar", function(varname, val = NULL) standardGeneric("regVar"))
setMethod("regVar", "character", function(varname, val = NULL) {
    assign(varname, val, envir = reg)
  })


## remove a variable
setGeneric("unregVar", function(varname) standardGeneric("unregVar"))
setMethod("unregVar", "character", function(varname) {
    if(exists(varname, envir = reg))
      rm(varname, envir = reg)
  })


## get the value of a variable
setGeneric("getVar", function(varname) standardGeneric("getVar"))
setMethod("getVar", "character", function(varname) {
	  if(exists(varname, envir = reg))
		  return(get(varname, envir = reg))
  })


## modify a variable
setGeneric("setVar", function(varname, val) standardGeneric("setVar"))
setMethod("setVar", "character", function(varname, val) {
	  if(exists(varname, envir = reg)) {
		  assign(varname, val, envir = reg)
		  hooks <- get(".hooks", envir = reg)
		  if(!is.null(funhooks <- hooks[[varname]]))
			  lapply(funhooks, function(fun) if(!is.null(fun)) fun())
	  }
  })


setGeneric("regEval", function(str) standardGeneric("regEval"))
setMethod("regEval", "character", function(str) {
	  ret <- eval(parse(text = str), envir = reg)
	  ret
  })


setGeneric("attachCallback", function(funname, fun, varname) standardGeneric("attachCallback"))
setMethod("attachCallback", c("character", "function", "character"),
	function(funname, fun, varname) {
		hooks <- get(".hooks", envir = reg)
		if(exists(varname, envir = reg)) {
			if(is.null(hooks[[varname]]))
				hooks[[varname]] <- list()
			hooks[[varname]][[funname]] <- fun
			assign(".hooks", hooks, envir = reg)
		}
	})


setGeneric("detachCallback", function(funname, varname) standardGeneric("detachCallback"))
setMethod("detachCallback", c("character", "character"),
	function(funname, varname) {
    hooks <- get(".hooks", envir = reg)
		if(exists(varname, envir = reg)) {
			if(!is.null(hooks[[varname]]))
				hooks[[varname]][[funname]] <<- NULL
		}
	}
)


setMethod("registerInterpolationMethod",
	c("character", "character"), function(name, classname) { 
		itrlst <- get("inter.list", envir = reg)
		itrlst[[name]]<- classname
		assign("inter.list", itrlst, envir = reg)
	}
)

setMethod("getInterList", "missing", function(dummy) get("inter.list", envir = reg))



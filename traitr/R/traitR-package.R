##  Copyright (C) 2010 John Verzani
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; either version 2 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  A copy of the GNU General Public License is available at
##  http://www.r-project.org/Licenses/

##' An interface for GUI creation using gWidgets
##'
##' This package provides an alternate interface for creating graphical user interfaces. The design was
##' inspired by the Traits UI module for python developed by enthought.com.
##' The implementation uses the MVC design pattern in the background, although the user need not be
##' aware of this.
##'
##' For basic use, the user creates a bunch of items (the model), specifies how these will be
##' layed out in a simple manner (the view), specifies actions to happen (the controller) and then creates a dialog.
##' See \code{\link{aDialog}} for examples.
##'
##' Creating basic dialogs requires no actual GUI programming
##' knowledge. One specifies the items by type of variable
##' (\code{numericItem} or \code{stringItem}, say) and the "action"
##' through a method call.
##'
##' The package uses the \pkg{proto} package so at some level, the R
##' user must use that OO syntax. In particular, methods calls are
##' done with the notation \code{obj$method_name} and method
##' definitions have an initial argument \code{.} for passing in a
##' reference to the \code{proto} object. (The name \code{.} is a
##' convention, but can be changed to \code{self} or \code{this}, if
##' that naming convention is preferred. Methods for \code{proto}
##' objects are documented (to some degree anyways). Their help page
##' is shown through the method \code{show_help}, as in
##' \code{obj$show_help()}.
##' 
##' @name traitR-package
##' @docType package
roxygen()

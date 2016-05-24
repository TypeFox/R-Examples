### tk2dialogs.R - Support for additional Tk dialogs
### Copyright (c), Philippe Grosjean (phgrosjean@sciviews.org)
### Licensed under LGPL 3 or above
###
### Changes:
### - 2007-02-07: fisrt version (for tcltk2_1.0-2)
###
### To do:
### - A good dialog box for choosing R colors (by number, name or RGB value)

tk2chooseFont <- function (...)
{
	if (!is.tk()) stop("Package Tk is required but not loaded")
	tclRequire("choosefont")
	tcl("choosefont::choosefont", ...)
}

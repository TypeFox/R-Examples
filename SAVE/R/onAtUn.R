##########################################################################
## start-up and clean-up functions
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2013-present Jesus Palomo, Gonzalo Garcia-Donato, 
##							  and Rui Paulo
##    
##########################################################################

.onAttach <- function(...) {
	
	date <- date()
	x <- regexpr("[0-9]{4}", date)
	this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
	
	# echo output to screen
	packageStartupMessage("#########")
	packageStartupMessage("##\n## SAVE Package for the statistical analysis of complex computer models (SAVE)")
	packageStartupMessage("## Copyright (C) 2013-", this.year,
			" Jesus Palomo, Gonzalo Garcia-Donato and Rui Paulo", sep="")
	packageStartupMessage("#########")
}

.onUnload <- function(libpath) {
	library.dynam.unload("SAVE", libpath)
}
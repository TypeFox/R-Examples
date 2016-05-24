#!/usr/bin/Rscript  --vanilla 
#
# SVMBridge 
#		(C) 2015, by Aydin Demircioglu
#
#		zzz.R
# 
# SVMBridge is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SVMBridge is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# Please do not use this software to destroy or spy on people, environment or things.
# All negative use is prohibited.
#


#' @useDynLib yakmoR
#' @importFrom Rcpp sourceCpp


# create own environment
yakmoREnv = new.env(parent = emptyenv())


.onLoad <- function(libname, pkgname) {
	op <- options()
	op.devtools <- list(
		devtools.path = "~/R-dev",
		devtools.install.args = "",
		devtools.name = "Aydin Demircioglu",
		devtools.desc.author = '"Aydin Demircioglu <aydin.demircioglu@ini.rub.de> [aut, cre]"',
		devtools.desc.license = "LGPL-3 + file LICENSE",
		devtools.desc.suggests = NULL,
		devtools.desc = list()
	)
	toset <- !(names(op.devtools) %in% names(op))
	if(any(toset)) options(op.devtools[toset])

	#  find_rtools()

	invisible()
}


.onAttach <- function (libname, pkgname) {
	packageStartupMessage("yakmoR loaded.")
}


.onUnload <- function (libpath) {
	library.dynam.unload("yakmoR", libpath)
}

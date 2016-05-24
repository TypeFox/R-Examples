#!/usr/bin/Rscript  --vanilla 
#
# LiblineaR.ACF
#		(C) 2015, by Aydin Demircioglu, Tobias Glasmachers, Urun Dogan,
#							Thibault Helleputte, Pierre Gramme
#
#		zzz.R
# 
# LiblineaR.ACF is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# LiblineaR.ACF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# Please do not use this software to destroy or spy on people, environment or things.
# All negative use is prohibited.
#


#' @useDynLib LiblineaR.ACF


# create own environment
LiblineaRACFEnv = new.env(parent = emptyenv())


.onLoad <- function(libname, pkgname) {
	op <- options()
	op.devtools <- list(
		devtools.path = "~/R-dev",
		devtools.install.args = "",
		devtools.name = "Aydin Demircioglu",
		devtools.desc.author = '"Aydin Demircioglu <aydin.demircioglu@ini.rub.de> [aut, cre]"',
		devtools.desc.license = "GPL-2",
		devtools.desc.suggests = NULL,
		devtools.desc = list()
	)
	toset <- !(names(op.devtools) %in% names(op))
	if(any(toset)) options(op.devtools[toset])

	invisible()
}


.onAttach <- function (libname, pkgname) {
	packageStartupMessage("LiblineaR.ACF loaded.")
}


.onUnload <- function (libpath) {
	library.dynam.unload("LiblineaR.ACF", libpath)
}

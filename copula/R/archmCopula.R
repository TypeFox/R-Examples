## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


archmCopula <- function(family, param = NA_real_, dim = 2L, ...) {
    family <- tolower(family)
    dim <- as.integer(dim)
    if(family == "amh" && dim != 2L)
	stop("'amh' is not yet available for dim > 2")
    switch(family,
	   "clayton" = claytonCopula(param, dim = dim, ...),
	   "frank"   = frankCopula  (param, dim = dim, ...),
	   "amh"     = amhCopula    (param, dim = dim, ...),
	   "gumbel"  = gumbelCopula (param, dim = dim, ...),
	   "joe"     = joeCopula    (param, dim = dim, ...),
	   ## otherwise:
	   { fams <- sub("Copula$", '', names(getClass("archmCopula")@subclasses))
	     stop("Valid family names are ", paste(dQuote(fams), collapse=", "))
	 })
}


tauArchmCopula <- function(copula) {
  1 + 4 * integrate(function(x) iPsi(copula, x) / diPsi(copula, x),
                    0, 1)$value
}

setMethod("tau", signature("archmCopula"), tauArchmCopula)

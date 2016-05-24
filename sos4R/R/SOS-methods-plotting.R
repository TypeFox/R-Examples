################################################################################
# Copyright (C) 2010 by 52 North                                               #
# Initiative for Geospatial Open Source Software GmbH                          #
#                                                                              #
# Contact: Andreas Wytzisk                                                     #
# 52 North Initiative for Geospatial Open Source Software GmbH                 #
# Martin-Luther-King-Weg 24                                                    #
# 48155 Muenster, Germany                                                      #
# info@52north.org                                                             #
#                                                                              #
# This program is free software; you can redistribute and/or modify it under   #
# the terms of the GNU General Public License version 2 as published by the    #
# Free Software Foundation.                                                    #
#                                                                              #
# This program is distributed WITHOUT ANY WARRANTY; even without the implied   #
# WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU #
# General Public License for more details.                                     #
#                                                                              #
# You should have received a copy of the GNU General Public License along with #
# this program (see gpl-2.0.txt). If not, write to the Free Software           #
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or #
# visit the Free Software Foundation web page, http://www.fsf.org.             #
#                                                                              #
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)                          #
# Created: 2011-02-09                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

#
#
#
plot.SosObservationOffering <- function(x, y, ..., add = FALSE) {
	.off.spatial <- NULL
	tryCatch({.off.spatial <- as(x, "Spatial") },
			error = function(e) { 
				warning(paste("Cannot not coerce offering", toString(sosId(x)),
								"to Spatial for plotting -- Error: ", e))
				
				return()
			})
	
	if(is.null(.off.spatial))
		warning("Cannot plot NULL offering!")
	else plot(x = .off.spatial, add = add, ...)
}
setMethod("plot", signature(x = "SosObservationOffering", y = "missing"),
		plot.SosObservationOffering)

#
#
#
plot.SOS <- function(x, y, ..., border.color.pal = sosDefaultColorPalette) {
	.offs <- sosOfferings(x)
	
	.args <- list(...)
	if(!is.null(.args[["add"]]))
		.addGiven <- TRUE
	else .addGiven <- FALSE
	
	for (i in seq(1, length(.offs))) {
		.off <- .offs[[i]]
		.add <- i != 1 # do not 'add' the first time
			
		if(!any(is.na(border.color.pal))) {
			.col <- border.color.pal[[(i %% length(border.color.pal)) + 1]]
			
			if(.addGiven) plot(x = .off, border = .col, ...)
			else plot(x = .off, border = .col, add = .add, ...)
		}
		else {
			if(.addGiven) plot(x = .off, ...)
			else plot(x = .off, add = .add, ...)
		}
	}
}
setMethod("plot", signature(x = "SOS", y = "missing"), plot.SOS)

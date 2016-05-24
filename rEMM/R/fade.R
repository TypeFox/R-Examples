#######################################################################
# rEMM - Extensible Markov Model (EMM) for Data Stream Clustering in R
# Copyrigth (C) 2011 Michael Hahsler
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


## build uses its own implementation of fade! See build.R

.fade <- function(x, t=1, lambda = NULL) {
	if(is.null(lambda)) lambda_factor <- x@tnn_d$lambda_factor
	else lambda_factor <- 2^(-lambda)

	## fade counts (tNN)
	x@tnn_d$counts <- x@tnn_d$counts * lambda_factor^t 

	## fade transition counts (TRACDS)
	x@tracds_d$mm <- smc_fade(x@tracds_d$mm, lambda_factor^t)
	
	invisible(x)
}


setMethod("fade", signature(x = "EMM", t= "numeric", lambda = "missing"),
	function(x, t, lambda) .fade(x, t)
)

setMethod("fade", signature(x = "EMM", t= "missing", lambda = "missing"),
	function(x, t, lambda) .fade(x, t=1)
)

setMethod("fade", signature(x = "EMM", t= "numeric", lambda = "numeric"),
	function(x, t, lambda) .fade(x, t, lambda) 
)

setMethod("fade", signature(x = "EMM", t= "missing", lambda = "numeric"),
	function(x, t, lambda) .fade(x, , lambda) 
)

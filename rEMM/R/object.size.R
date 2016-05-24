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


### This is needed to calculate the sice of the elements in the environments
### correctly

setMethod("object.size", signature(x = "TRACDS"),
	function(x) {
	    s <- sum(sapply(ls(x@tracds_d),  function(y)
			    object.size(get(y, envir = x@tracds_d))))

	    class(s) <- "object_size"
	    s
	})


setMethod("object.size", signature(x = "tNN"),
	function(x) {
	    s <- sum(sapply(ls(x@tnn_d),  function(y)
			    object.size(get(y, envir = x@tnn_d))))

	    class(s) <- "object_size"
	    s
	})

setMethod("object.size", signature(x = "EMM"),
	function(x) {
	    s <- sum(sapply(ls(x@tracds_d),  function(y)
			    object.size(get(y, envir = x@tracds_d))))
	    
	    s <- s + sum(sapply(ls(x@tnn_d),  function(y)
			    object.size(get(y, envir = x@tnn_d))))

	    s <- s+object.size(unclass(x))

	    class(s) <- "object_size"
	    s
	})




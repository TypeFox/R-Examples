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


setMethod("c", signature(x = "EMM"),
	function(x, ..., copy=TRUE, recursive = FALSE) {
	    args <- list(...)

	    if(copy) x <- copy(x)

	    if (recursive)
		args <- unlist(args)
	    for (y in args) {
		if (!is(y, "EMM"))
		    stop("can combine EMM only")


		## combine tracds
		nx <- smc_size(x@tracds_d$mm)
		ny <- smc_size(y@tracds_d$mm)
		n <- nx + ny
		
		new_labels <- as.character(1:n)
		#new_labels <- c(clusters(x), paste("A", clusters(y)))

		x@tracds_d$mm <- new("SimpleMC",
			unused = rep(as.integer(NA), n),
			top = 0L,
			counts = structure(rbind( 
				cbind(
					smc_countMatrix(x@tracds_d$mm), 
					matrix(0, ncol=ny, nrow=nx)),
				cbind(
					matrix(0, ncol=nx, nrow=ny),
					smc_countMatrix(y@tracds_d$mm)
				    )),
				dimnames=list(new_labels,new_labels)
				),

			initial_counts = structure(
				c(smc_initialCounts(x@tracds_d$mm),
				smc_initialCounts(y@tracds_d$mm)),
				names=new_labels)
			)
		
		x@tracds_d$current_state <- NA

		## combine tNN
		x@tnn_d$centers <- structure(
			rbind(x@tnn_d$centers, y@tnn_d$centers),
			dimnames = list(new_labels, colnames(x@tnn_d$centers))) 
		
		x@tnn_d$counts <- structure(
			c(x@tnn_d$counts, y@tnn_d$counts),
			names = new_labels
			)
		
		x@tnn_d$var_thresholds <- structure(
			c(x@tnn_d$var_thresholds, 
				y@tnn_d$var_thresholds),
			names = new_labels
			)

		x@tnn_d$last <- NA
	    }
	
	
	    invisible(x)
	})


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


## constructor
EMM <- function(threshold=0.2, measure="euclidean", distFun = NULL, 
	centroids=identical(tolower(measure), "euclidean"), 
	lambda=0, data=NULL) {

	emm <- new("EMM", measure=measure, distFun = distFun, 
		threshold=threshold, centroids=centroids, lambda=lambda)

	if(!is.null(data)) build(emm, data)
	emm

    }

### deep copy
setMethod("copy", signature(x = "EMM"),
	function(x) {

	    r <- new("EMM",
		    threshold = x@threshold,
		    measure = x@measure,
		    distFun = x@distFun,
		    centroids = x@centroids,
		    lambda = x@lambda,
		    lambda_factor = x@lambda_factor)

	    ## copy environments
	    r@tnn_d <- as.environment(as.list(x@tnn_d))
	    r@tracds_d <- as.environment(as.list(x@tracds_d))

	    r
	})


## show
setMethod("show", signature(object = "EMM"),
        function(object) {
            cat("EMM with", size(object), "states/clusters.\n", 
                    "Measure:", object@measure, "\n",
                    "Threshold:", object@threshold, "\n",
                    "Centroid:", object@centroids, "\n",
                    "Lambda:", object@lambda, "\n"
                    )
            invisible(NULL)
        })

## size delegates to nclusters (in tNN)
setMethod("size", signature(x = "EMM"),
	        function(x) nclusters(x)
		)



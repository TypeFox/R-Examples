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

setClass("StreamClustering")

setClass("tNN",
	contains = ("StreamClustering"),
	representation(
		measure		= "character",
		distFun		= "ANY",
		centroids	= "logical",
		threshold	= "numeric",
		lambda		= "numeric",
		lambda_factor	= "numeric",
	        ### this is all in an environment now	
		#centers		= "matrix",	## row names are cluster names
		#counts		= "numeric",
		#var_thresholds	= "numeric",
		#last		= "character"
		tnn_d		= "environment"
	)

	## FIXME: Implement check
	#validity= function(object) {}
)

setMethod("initialize", "tNN", function(.Object, 
		threshold = 0.2, 
		measure = "euclidean", 
		distFun = NULL,
		centroids = identical(tolower(measure), "euclidean"), 
		lambda=0, ...){
	    
	    .Object@threshold <- threshold
	    .Object@measure <- measure
	    .Object@centroids <- centroids
	    .Object@lambda <- lambda
	    .Object@lambda_factor <- 2^(-lambda)
		
	    .Object@tnn_d <- new.env()
	    assign("centers", matrix(nrow=0, ncol=0), envir = .Object@tnn_d)
	    assign("counts", numeric(), envir = .Object@tnn_d)
	    assign("var_thresholds", numeric(), envir = .Object@tnn_d)
	    assign("last", as.character(NA), envir = .Object@tnn_d)
	    
	    ### get dist function from proxy registry
	    ### Note: this makes calling dist in proxy faster (no more lockup) 
	    if(!is.null(distFun)) .Object@distFun <- distFun
	    else .Object@distFun <- pr_DB[[measure]]
	    
	    #validObject(.Object)
	    #.Object <- callNextMethod(.Object, ...)
	    
	    .Object
	})



.smc_size <- 10L
setClass("SimpleMC",
	representation(
		unused      = "integer", ## list of unused cols/rows
		top         = "integer", ## top of unused
		counts      = "matrix",
		initial_counts = "numeric" 
		),
	
	prototype(
		unused	    = .smc_size:1,
		top	    = .smc_size,
		counts	    = matrix(0, ncol=.smc_size, nrow=.smc_size),
		initial_counts = structure(rep(0, .smc_size), 
			names=rep(NA, .smc_size))  ## also holds cluster names
		)

	## FIXME: Implement check
	#validity= function(object) {
	#}
	)

setClass("TRACDS",
	representation(
		lambda		= "numeric",
		lambda_factor	= "numeric",
		tracds_d    	= "environment"
		#mm		= "SimpleMC", 
		#current_state	= "character"
		),

	prototype(
		lambda		= 0,
		lambda_factor	= 1,
		tracds_d    	= emptyenv()
		#mm		= new("SimpleMC"),
		#current_state	= as.character(NA)
		),

	## FIXME: Implement check
	#validity= function(object) {
	#}
	)

setMethod("initialize", "TRACDS", function(.Object, lambda=0, ...){
	    
	    .Object@lambda <- lambda
	    .Object@lambda_factor <- 2^(-lambda)

	    .Object@tracds_d <-  new.env()
	    assign("mm", new("SimpleMC"), envir = .Object@tracds_d)
	    assign("current_state", as.character(NA), envir = .Object@tracds_d)
	    
	    #cat("TRACDS initializes.\n")
	    #validObject(.Object)
	    
	    #.Object <- callNextMethod(.Object, ...)
	    
	    .Object
	})


### EMM is a subclass of tNN and TRACDS
setClass("EMM", contains = c("TRACDS", "tNN"))

## S4 BUG!!! R seems to only call initialize for the first class in contains
setMethod("initialize", "EMM", function(.Object, ...){
	    #threshold, measure, lambda,
	    #	centroids, ...){
	    
	    #.Object <- callNextMethod(.Object, ...)
	    
	    #.Object@threshold <- threshold
	    #.Object@measure <- measure
	    #.Object@centroids <- centroids
	    #.Object@lambda <- lambda
	    #.Object@lambda_factor <- 2^(-lambda)

	    .Object <- getMethod("initialize", "TRACDS")(.Object, ...)	 
	    .Object <- getMethod("initialize", "tNN")(.Object, ...)	 

	    validObject(.Object)
	    
	    .Object
	})

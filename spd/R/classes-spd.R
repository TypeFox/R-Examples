#################################################################################
##
##   R package spd by Alexios Ghalanos Copyright (C) 2008-2013
##   This file is part of the R package spd.
##
##   The R package spd is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package spd is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################
# Developer Notes: SPD is the superclass which can hold different tail distributions
setClass("SPD","VIRTUAL")

setClass("GPDTAILS",
		representation(
				call = "call",
				method = "character",
				kernel = "character",
				data = "matrix",
				threshold = "list",
				ptails = "list",
				fit = "list",
				title = "character",
				description = "character"),
		contains="SPD")

setClass("GEVTAILS",
		representation(
				call = "call",
				method = "character",
				kernel = "character",
				data = "matrix",
				threshold = "list",
				ptails = "list",
				fit = "list",
				title = "character",
				description = "character"),
		contains="SPD")
#---------------------------------------------------------------------
setMethod("show",
		signature(object="SPD"),
		function(object){
			cat("\nTitle:\n", object@title)
			# Function Call:
			cat(paste("\nTail Fit: ",substr(is(object)[1],1,3),sep=""))
			cat(paste("\nKernel Fit: ",object@kernel,sep=""))
			# Estimation Type:
			cat("\nTail Estimation Method: ",object@method,"\n")
			# Estimated Parameters:
			cat("\nUpper Tail:\n")
			cat("-----------------------------------------\n")
			cat("Estimated Parameters:\n")
			print(round(object@fit$upperFit@fit$par.est,5))
			cat(paste("Threshold: ",round(object@threshold$upper,5),sep=""),"\n")
			cat("-----------------------------------------\n")
			cat("\nLower Tail:\n")
			cat("-----------------------------------------\n")
			cat("Estimated Parameters:\n")
			print(round(object@fit$lowerFit@fit$par.est,5))
			cat(paste("Threshold: ",round(object@threshold$lower,5),sep=""),"\n")
			cat("-----------------------------------------\n")
			# Desription:
			cat("\nDescription\n ", object@description, "\n\n")

			# Return Value:
			invisible(object)

		})
#---------------------------------------------------------------------
setClass("GPDFIT",
		representation(
				call = "call",
				method = "character",
				parameter = "list",
				data = "list",
				fit = "list",
				residuals = "numeric",
				title = "character",
				description = "character"
		)
)

setMethod("show",
		signature(object="GPDFIT"),
		function(object){
			cat("\nTitle:\n ", object@title, "\n")
			# Function Call:
			cat("\nCall:\n ")
			cat(paste(deparse(object@call), sep = "\n",
							collapse = "\n"), "\n", sep = "")
			# Estimation Type:
			cat("\nEstimation Method:\n ", object@method, "\n")
			# Estimated Parameters:
			cat("\nEstimated Parameters:\n")
			print(object@fit$par.ests)
			# Desription:
			cat("\nDescription\n ", object@description, "\n\n")
			# Return Value:
			invisible(object)
		})
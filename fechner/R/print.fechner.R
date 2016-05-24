###############################################
##print method for objects of class "fechner"##
###############################################
print.fechner <-
function(x, ...){
	# x: an object of class "fechner", obtained from a call to function fechner
	# ...: further arguments to be passed to or from other methods; they are ignored in this function

	if (attr(x, which = "computation", exact = TRUE) == "short"){  # corresponds to "fechner" object resulting from short computation
		# overall Fechnerian distances and geodesic loops printed are of the first kind (or, in the first observation area); in particular,
		# geodesic loops must be read accordingly, from left to right for the first kind and from right to left for the second kind
		# (or, in the second observation area)
		cat("\n")
		cat("overall Fechnerian distances:\n")
		print(x$overall.Fechnerian.distances)
		cat("\n")
		cat("geodesic loops:\n")
		print(x$geodesic.loops)
		cat("\n")
		invisible(x)
	} else
		if (attr(x, which = "computation", exact = TRUE) == "long"){  # corresponds to "fechner" object resulting from long computation
			# overall Fechnerian distances and geodesic loops printed are of the first kind (or, in the first observation area); in particular,
			# geodesic loops must be read accordingly, from left to right for the first kind and from right to left for the second kind
			# (or, in the second observation area)
			cat("\n")
			cat("overall Fechnerian distances:\n")
			print(x$overall.Fechnerian.distances.1)
			cat("\n")
			cat("geodesic loops:\n")
			print(x$geodesic.loops.1)
			cat("\n")
			invisible(x)
		} else
			stop("object attribute computation must have value \"short\" or \"long\"")
}

###############################################
##plot method for objects of class "fechner"###
###############################################
plot.fechner <-
function(x, level = 2, ...){
	# x: an object of class "fechner", obtained from a call to function fechner
	# level: the 'level' of comparison of overall Fechnerian distance, G, and S-index, S; it refers to the minimum number
	#        of links in geodesic loops; that is, choosing level n means that comparison involves only those G and
	#        and S values that have geodesic loops containing not less than n links
	# ...: further arguments to be passed to or from other methods; they are ignored in this function

	if (mode(level) != "numeric")
		stop("level must be number")
	if (!is.finite(level))
		stop("level must be finite, i.e., not be NA, NaN, Inf, or -Inf")
	if (as.integer(level) != level)
		stop("level must be integer")
	if (level < 2)
		stop("level must be greater than or equal to 2")

	if (attr(x, which = "computation", exact = TRUE) == "short"){  # corresponds to "fechner" object resulting from short computation
		G <- x$overall.Fechnerian.distances  # obtained from computations of the first kind
		S <- x$S.index
		number.links <- x$graph.lengths.of.geodesic.loops  # obtained from computations of the first kind
	} else
		if (attr(x, which = "computation", exact = TRUE) == "long"){  # corresponds to "fechner" object resulting from long computation
			G <- x$overall.Fechnerian.distances.1  # obtained from computations of the first kind
			S <- x$S.index
			number.links <- x$graph.lengths.of.geodesic.loops.1  # obtained from computations of the first kind
		} else
			stop("object attribute computation must have value \"short\" or \"long\"")

	pairs <- (number.links[upper.tri(number.links)] >= level)
	if (!any(pairs))
		stop(paste("plot is not possible: there are no (off-diagonal) pairs of stimuli with geodesic loops containing at least ",
					 level, " links", sep = ""))
	G.level <- G[upper.tri(G)][pairs]
	S.level <- S[upper.tri(S)][pairs]
	plot(S.level, G.level, xlim = c(0, ceiling(max(S.level, G.level))), ylim = c(0, ceiling(max(S.level, G.level))),
		 main = paste("Scatterplot \"(overall) Fechnerian distance G versus S-Index\"\n(for comparison level ", level,
		 			  ", with diagonal line y = x)", sep = ""), xlab = "S-index", ylab = "Fechnerian distance G", pch = 19)
	rug(jitter(S.level, amount = 0.01), ticksize = 0.0225, col = gray(0.075))
	rug(jitter(G.level, amount = 0.01), side = 2, ticksize = 0.0225, col = gray(0.075))
	abline(0, 1, col = "lightblue")

	invisible(NULL)
}

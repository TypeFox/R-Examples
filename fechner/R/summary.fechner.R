#################################################
##summary method for objects of class "fechner"##
#################################################
summary.fechner <-
function(object, level = 2, ...){
	# object: an object of class "fechner", obtained from a call to function fechner
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

	if (attr(object, which = "computation", exact = TRUE) == "short"){  # corresponds to "fechner" object resulting from short computation
		G <- object$overall.Fechnerian.distances  # obtained from computations of the first kind
		S <- object$S.index
		number.links <- object$graph.lengths.of.geodesic.loops  # obtained from computations of the first kind
	} else
		if (attr(object, which = "computation", exact = TRUE) == "long"){  # corresponds to "fechner" object resulting from long computation
			G <- object$overall.Fechnerian.distances.1  # obtained from computations of the first kind
			S <- object$S.index
			number.links <- object$graph.lengths.of.geodesic.loops.1  # obtained from computations of the first kind
		} else
			stop("object attribute computation must have value \"short\" or \"long\"")

	pairs <- (number.links[upper.tri(number.links)] >= level)
	if (!any(pairs))
		stop(paste("summary is not possible: there are no (off-diagonal) pairs of stimuli with geodesic loops containing at least ",
				   level, " links", sep = ""))
	G.level <- G[upper.tri(G)][pairs]
	S.level <- S[upper.tri(S)][pairs]

	correlation <- cor(S.level, G.level)  # Pearson's correlation coefficient
	if (is.na(correlation))
		correlation <- "Pearson's correlation coefficient is not defined"
	C.index <- ((2 * sum((S.level - G.level)^2)) / (sum((S.level)^2) + sum((G.level)^2)))  # C-index

	# data frame of the pairs of stimuli and the corresponding S and G values used for comparison
	stimuli.pairs <- paste(rownames(S)[row(S)[upper.tri(S)][pairs]], ".", colnames(S)[col(S)[upper.tri(S)][pairs]], sep = "")
	comparison.pairs <- data.frame(stimuli.pairs = stimuli.pairs, S.index = S.level, Fechnerian.distance.G = G.level, stringsAsFactors = FALSE)

	results <- list(pairs.used.for.comparison = comparison.pairs,
					Pearson.correlation = correlation,
					C.index = C.index,
					comparison.level = level)

	# introduces the class "summary.fechner", to provide a print method to print summary information about objects of class "fechner"
	class(results) <- "summary.fechner"

	return(results)
}

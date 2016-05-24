##########################################################################
##provides all Fechnerian computations; the main function of the package##
##########################################################################
fechner <-
function(X, format = c("probability.different", "percent.same", "general"), compute.all = FALSE, check.computation = FALSE){
	# X: the data
	# format: the data format that is used; "probability.different" and "percent.same" are for datasets in the probability-different
	#         and percent-same formats, and in the latter case, the data are automatically transformed prior to the analysis using
	#         transformation (100 - X) / 100; "general" is to be used for datasets that are properly in the general data format;
	#         the default is "probability.different"
	#         NOTE: for "percent.same", the input data must satisfy regular maximality, for "probability.different" and "general",
	#               regular minimality (otherwise respective messages are produced); in particular, data
	#               in the general format may possibly need to be transformed manually prior to calling function fechner
	# compute.all: the default value FALSE corresponds to short computation which yields the main Fechnerian computations;
	#              TRUE corresponds to long computation which additionally yields intermediate results and also allows for a check
	#              of computations if check.computation is set TRUE
	# check.computation: check for whether the overall Fechnerian distance of the first kind (or, in the first observation area),
	#                    overall.G.1, is equal to the overall Fechnerian distance of the second kind (or, in the second observation area),
	#                    overall.G.2, which holds by theory; the difference overall.G.1 - overall.G.2 is computed; the default is FALSE
	#                    NOTE: requires compute.all to be set TRUE

	if ((!is.logical(compute.all) || is.na(compute.all)) || (!is.logical(check.computation) || is.na(check.computation)))
		stop("specify FALSE or TRUE for logical formal arguments")
	if ((!compute.all) && check.computation)
		stop("it is not possible to check computations without computing all")

	format <- match.arg(format)
	if (format == "probability.different"){
		Canonical <- check.regular(X, type = "probability.different")
		M <- Canonical$canonical.representation
	}
	if (format == "percent.same"){
		Canonical <- check.regular(X, type = "percent.same")
		M <- ((100 - Canonical$canonical.representation) / 100)
	}
	if (format == "general"){
		Canonical <- check.regular(X, type = "reg.minimal")
		M <- Canonical$canonical.representation
	}
	if (is.null(Canonical))
		invisible(NULL)
	else{
		PSE <- Canonical$canonical.transformation
		n <- dim(M)[1]

		# matrices of the psychometric increments of the first and second kind (or, in the first and second observation areas),
		# phi.1 and phi.2, respectively
		phi.1 <- phi.2 <- matrix(nrow = n, ncol = n)
		for (i in 1:n){
			for (j in 1:n){
				phi.1[i, j] <- (M[i, j] - M[i, i])
				phi.2[i, j] <- (M[j, i] - M[i, i])
			}
		}
		dimnames(phi.1) <- dimnames(phi.2) <- dimnames(M)

		# generalized "Shepardian" dissimilarity index (or S-index); it is defined as the psychometric length of the loop
		# between the row stimulus and the column stimulus containing only these two stimuli
		S.index <- matrix(nrow = n, ncol = n)
		for (i in 1:n){
			for (j in 1:n){
				S.index[i, j] <- (M[i, j] + M[j, i] - M[i, i] - M[j, j])
			}
		}
		dimnames(S.index) <- dimnames(M)

		if (!compute.all){
			information.1 <- shortest.paths.information(phi.1)

			# matrix of the oriented Fechnerian distances of the first kind (or, in the first observation area), oriented.G.1
			oriented.G.1 <- information.1$weight.distances
			# overall Fechnerian distances of the first kind (or, in the first observation area), overall.G.1
			overall.G.1 <- (oriented.G.1 + t(oriented.G.1))
			# matrix of the predecessors of the column stimuli in shortest paths from the row stimuli (as source vertices) to the column stimuli
			# (as target vertices); of the first kind (or, in the first observation area)
			predecessors.1 <- information.1$predecessors
			# geodesic loops of the first kind (or, in the first observation area), geodesic.loops.1; must be read from left to right
			# for the first kind (or, first observation area), and from right to left for the second kind (or, second observation area)
			geodesic.loops.1 <- geodesic.chains.loops(predecessors.1)$geodesic.loops
			# matrix of the graph-theoretic (edge/link based) lengths of the shortest paths from source vertices (stimuli) to target vertices (stimuli);
			# of the first kind (or, in the first observation area)
			graph.theoretic.lengths.1 <- information.1$edge.distances

			results <- list(points.of.subjective.equality = PSE,
							canonical.representation = M,
							overall.Fechnerian.distances = overall.G.1,
							geodesic.loops = geodesic.loops.1,
							graph.lengths.of.geodesic.loops = (graph.theoretic.lengths.1 + t(graph.theoretic.lengths.1)),
							S.index = S.index)
		} else{
			information.1 <- shortest.paths.information(phi.1)
			information.2 <- shortest.paths.information(phi.2)

			# matrices of the oriented Fechnerian distances of the first and second kind (or, in the first and second observation area),
			# oriented.G.1 and oriented.G.2
			oriented.G.1 <- information.1$weight.distances
			oriented.G.2 <- information.2$weight.distances
			# overall Fechnerian distances of the first and second kind (or, in the first and second observation area),
			# overall.G.1 and overall.G.2
			overall.G.1 <- (oriented.G.1 + t(oriented.G.1))
			overall.G.2 <- (oriented.G.2 + t(oriented.G.2))
			# matrices of the predecessors of the column stimuli in shortest paths from the row stimuli (as source vertices) to the column stimuli
			# (as target vertices); of the first and second kind (or, in the first and second observation area)
			predecessors.1 <- information.1$predecessors
			predecessors.2 <- information.2$predecessors

			geodesics.1 <- geodesic.chains.loops(predecessors.1)
			geodesics.2 <- geodesic.chains.loops(predecessors.2)

			# geodesic chains of the first and second kind (or, in the first and second observation area), geodesic.chains.1 and geodesic.chains.2
			geodesic.chains.1 <- geodesics.1$geodesic.chains
			geodesic.chains.2 <- geodesics.2$geodesic.chains
			# geodesic loops of the first and second kind (or, in the first and second observation area), geodesic.loops.1 and geodesic.loops.2
			geodesic.loops.1 <- geodesics.1$geodesic.loops
			geodesic.loops.2 <- geodesics.2$geodesic.loops
			# matrices of the graph-theoretic (edge/link based) lengths of the shortest paths from source vertices (stimuli) to target vertices (stimuli);
			# of the first and second kind (or, in the first and second observation area)
			graph.theoretic.lengths.1 <- information.1$edge.distances
			graph.theoretic.lengths.2 <- information.2$edge.distances

			# checking computations; difference, and in terms of computational accuracy (i.e., up to machine precision)
			if (check.computation){
				difference <- (overall.G.1 - overall.G.2)  # differences
				are.nearly.equal <- isTRUE(all.equal(difference, matrix(0, nrow = n, ncol = n), check.attributes = FALSE))  # up to machine precision
				check <- list(difference = difference, are.nearly.equal = are.nearly.equal)
			} else
				check <- "computation check was not requested"

			results <- list(points.of.subjective.equality = PSE,
							canonical.representation = M,
							psychometric.increments.1 = phi.1,
							psychometric.increments.2 = phi.2,
							oriented.Fechnerian.distances.1 = oriented.G.1,
							overall.Fechnerian.distances.1 = overall.G.1,
							oriented.Fechnerian.distances.2 = oriented.G.2,
							overall.Fechnerian.distances.2 = overall.G.2,
							check = check,
							geodesic.chains.1 = geodesic.chains.1,
							geodesic.loops.1 = geodesic.loops.1,
							graph.lengths.of.geodesic.chains.1 = graph.theoretic.lengths.1,
							graph.lengths.of.geodesic.loops.1 = (graph.theoretic.lengths.1 + t(graph.theoretic.lengths.1)),
							geodesic.chains.2 = geodesic.chains.2,
							geodesic.loops.2 = geodesic.loops.2,
							graph.lengths.of.geodesic.chains.2 = graph.theoretic.lengths.2,
							graph.lengths.of.geodesic.loops.2 = (graph.theoretic.lengths.2 + t(graph.theoretic.lengths.2)),
							S.index = S.index)
		}

		# sets a named attribute indicating whether short or long computation was used
		attr(results, "computation") <- ifelse(compute.all, "long", "short")
		# introduces the class "fechner", to provide plot, print, and summary methods
		# NOTE: objects of class "fechner" are assumed to have the attribute computation, with value "short" or "long")
		class(results) <- "fechner"

		return(results)
	}
}

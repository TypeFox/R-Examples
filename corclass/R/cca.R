###############################################################################
#################### CORRELATIONAL CLASS ANALYSIS  (CCA) #####################
# This is an implementation of the CCA algorithm [1], which improves the RCA 
#  algorithm from Goldberg (2011).  Like RCA, it divides a population into groups 
#  or "schematic classes" based on similarity of correlation patterns.  It differs 
#  from RCA in using absolute row correlations instead of relationalities.  
#  This greatly increases the accuracy of the algorithm[1].  It also speeds it 
#  up a great deal, and makes it more analytically tractable.
# 
# Description of the algorithm and why it improves RCA's accuracy: 
# [1]  Boutyline, Andrei. Under review.
#		"Correlational Class Analysis Draws a More Accurate Map 
#			(Comment on Goldberg)." 
#       Available from www.ocf.berkeley.edu/~andrei.
#
# This is the first public release, from 11.25.2013.
#
#    Andrei Boutyline
#    University of California, Berkeley
#    andrei.boutyline@gmail.com
###############################################################################


# package corclass

## MAIN FUNCTION:
cca <- function (dtf, filter.significance = TRUE, 
		filter.value = 0.01,
		zero.action = c("drop", "ownclass"), 
		verbose = TRUE) {
	# CORRELATIONAL CLASS ANALYSIS (CRCA) 
	# PURPOSE:  Divide observations into schematic classes based on row correlations.
	# PARAMETERS:
	# 	* dtf : the data frame containing the variables of interest
	#   * filter.significance : Significance filtering sets "insignificant" ties to 0 to decrease noise 
	#     and increase stability.  Simulation results show that this greatly increases accuracy in many settings.
	#     Set filter.significance = FALSE to disable this.
	#   * zero.action : what to do with 0-variance rows before partitioning the graph. 
	#       zero.action = "drop"     : drop rows with 0 variance from the analyses (default).
	#       zero.action = "ownclass" : the correlations between 0-var rows and all other rows
	#     	is set to 0, and the correlations between all pairs of 0-var rows are set to 1. 
	#       This effectively creates a "zero class".
	#   * verbose : whether to print details of what CCA is doing to the screen
	#  RETURN VALUES:
	#   * Class memberships are returned in $membership (see help(leading.eigenvector.community) for details).
	#   * The whole correlation matrix is in $cormat (minus the 0-variance rows, if they were dropped).  
	#       It contains a "dtf" attribute which holds the whole dataframe (minus dropped rows)
	#       It also contains a "zeros" attribute which holds the indexes of the dropped rows.
	#   * For convenience, the dataframe is separated into the modules found by the algorithm.  
	#         These individual modules (classes) are in $modules.  
	#         The rows corresponding for each module i are in $modules[[i]]$dtf.
	#         The column correlations for each module i are in $modules[[i]]$cormat.
	#         These modules can be plotted via the S3 method below by simply calling plot(module).

	
	if (verbose) echo <- cat 
	else echo <- c
		
	if (packageVersion("igraph") < 0.7) {
		warning(paste("Your igraph version is ", as.character(packageVersion("igraph")), 
						". CCA produces more accurate results with igraph >= 0.7.\n", sep = ""))
	}
	
	cormat <- .make.cormat(dtf, zero.action)
	if (filter.significance == TRUE) {
			echo("Filtering out correlations for which Pr(|r| != 0) > ", filter.value, "\n")
			cormat <- .filter.insignif (cormat, ncol(dtf), pcutoff = filter.value)
	} else {
			echo("Not filtering significances.\n")
	}
	
	graph <- .cormat.to.igraph (cormat, absolute.value = TRUE)
	comm <- igraph::leading.eigenvector.community(graph)
			
	modules <- .separate(attr(cormat, "dtf"), comm$membership) 
	
	val <- list (membership = comm$membership, modules = modules, cormat = cormat)
	class (val) <- "cca"

	echo (paste(capture.output(print(val)), collapse = "\n"), "\n")
	
	return (invisible(val))
}

print.cca <- function (x, ...) {
	cat("CCA found", length(unique(x$membership)), "schematic classes. Sizes:", table(x$membership), "\n")
	degen <- sapply(x$modules, function (m1) m1$degenerate)
	if (any(degen)) {
		cat("NOTE: result contains ", sum(degen), " degenerate class(es): ", paste("#", which(degen), sep = "", collapse = " "), ".\n", sep="")
	}
}

plot.cca <- function (x, module.index, cutoff = 0.05, LAYOUT = igraph::layout.kamada.kawai, 
		drop.neg.ties.for.layout = TRUE, bw = FALSE, main = NULL, file = NULL, ...) {
	## PLOTTING FUNCTION for objects returned by cca().
	##             It is invoked when plot () is called on an object of the "cca" class.                             
	## PARAMETERS: 
	##    x  : the cca object returned by cca().  
	##    module.index : index of the module to plot.
	##    cutoff : the minimum absolute value of correlations to plot
	##    LAYOUT : either an igraph layout function, or a static layout in the same format
	##             as created by such a function.
	##    drop.neg.ties.for.layout : whether to drop negative ties for the purpose of layout.
	##       This is necessary because some layout algorithms crash if negative ties are present.
	##       Note that the dropped ties are still included in the actual plot.
	##    bw     : whether to print in color for screen viewing (FALSE), or in b&w with
	##             dashed lines for negative ties for a journal manuscript (TRUE)
	##    main  : caption at the top of the graph.  The default is the module number
	##    file  : if a filename is provided, the graph is saved as a pdf with that filename.
	##            Note that this requires the Cairo library.
	##  RETURN VALUE:
	##    The static layout that was used for this plot is returned.  This allows the same exact layout to be reproduced in the future.
	##    E.g., layout1 <- plot(cca.object, 1);  plot(cca.object, 1, LAYOUT = layout1).
	if (missing(module.index)) {
		warning("No module index provided. Plotting first module.")
		module.index <- 1
	}
	
	if(x$modules[[module.index]]$degenerate == TRUE) {
		stop(paste("Module #", module.index, " is degenerate (one or more column correlations are undefined).", sep = ""))
	}
	
	cormat <- x$modules[[module.index]]$cormat
	cormat[abs(cormat) < 0.05] <- 0
	
	G <- .cormat.to.igraph (cormat, FALSE)
	
	if (drop.neg.ties.for.layout == TRUE)
		G.forlayout <- 	igraph::delete.edges(G, igraph::E(G)[igraph::E(G)$weight < 0])	
	else
		G.forlayout <- G
	
	if (is.function(LAYOUT))
		G$layout <- LAYOUT(G.forlayout, weights = igraph::E(G.forlayout)$weight)
	else
		G$layout <- LAYOUT
	
	if (is.null(main)) {
		main <- paste("Module #", module.index, sep = "")
	}
	
	if (bw) {
		edge.lty <- c(2,1)[1 + (igraph::E(G)$weight > 0)]
		edge.color <- "black"
		vertex.color <- "white"
		vertex.frame.color <- "black"
		label.color <- "black"
	}
	else {
		edge.lty <- 1
		edge.color <- c(rgb(0.86,0.08,0.24,0.5), rgb(0.44,0.44,0.78,0.5))[1 + (igraph::E(G)$weight > 0)]
		vertex.color <- rgb(1,1,1,0.7)
		vertex.frame.color <- rgb(0.44,0.44,0.78, 1)
		label.color <- "black"
	}
	edge.width <- (igraph::E(G)$weight) * 10
	
	plot(G, edge.width = edge.width, edge.color = edge.color, edge.lty = edge.lty, 
			vertex.color = vertex.color, vertex.frame.color = vertex.frame.color, 
			vertex.size = 20, vertex.label.color = label.color, vertex.label.cex = 0.75, 
			main = main)
	
	if (!is.null(file)) {
		if (requireNamespace("Cairo", quietly = TRUE)) {
			Cairo::CairoPDF(file)
		} else {
			pdf(file)
		}

		plot(G, edge.width = edge.width, edge.color = edge.color, edge.lty = edge.lty, 
				vertex.color = vertex.color, vertex.frame.color = vertex.frame.color, 
				vertex.size = 20, vertex.label.color = label.color, vertex.label.cex = 0.5, 
				main = main)
		dev.off()
	}
	
	invisible(G$layout)
}


##### HELPER FUNCTIONS BELOW THIS LINE #########
.make.cormat <- function (dtf, zero.action) {
	# Helper function.  Make a correlation matrix from data frame.  
	if (!all(sapply(dtf, is.numeric))) dtf2 <- data.frame(sapply(dtf, as.numeric))
	else dtf2 <- dtf
	
	# Floating point imprecision may make 0-variance rows appear to have variance slightly higher than 0.
	zeros <- which(apply(dtf2, 1, var) <= 0.000000001)
	
	if (zero.action[1] == "drop" & (length(zeros) > 0)) {
		dtf2 <- dtf2[-zeros,]	
	}
	
	rv <- abs(cor(t(dtf2)))
	
	attributes(rv)$zeros <- zeros
	attributes(rv)$zero.action <- zero.action[1]
	attributes(rv)$dtf <- dtf2
	
	if ((zero.action[1] == "ownclass") & length(zeros) > 0) {
		rv[zeros,] <- 0
		rv[,zeros] <- 0
		rv[zeros,zeros] <- 1
	}
	
  	diag(rv) <- 0
	
	return (rv)
}


.filter.insignif <- function (corr, N.vars, pcutoff = 0.05) {
	# Helper function.
	# Filter signififances at p <= pcutoff (two-tailed).
	corr <- abs(corr)
	
	if (any(diag(corr) != 0))
		stop("Non-zero elements on the diagonal. diag(corr) <- 0 before running this function.")
	
	suppressWarnings(tvalues <- corr * sqrt ((N.vars-2) / (1 - corr^2)))
	if (any(is.infinite(tvalues))) {
		tvalues[is.infinite(tvalues)] <- 9999 # a very big number
	}
	cutoff <- abs(qt(pcutoff / 2, N.vars))
	
	isolates.pre <- sum(apply(corr, 1, sum) == 0)
	corr[tvalues < cutoff] <- 0
	isolates.post <- sum(apply(corr, 1, sum) == 0)
	
	if (isolates.post > isolates.pre) {
		warn1 <- paste ("Significance filtering left", isolates.post - isolates.pre, "rows with no non-zero ties. The CCA result will contain at least one small degenerate class.")
		warning(warn1)
	}
	
	return (corr)
}

.cormat.to.igraph <- function (corr, absolute.value = TRUE) { 
	# Helper function.
	# Make igraph object from correlation matrix.
	if (absolute.value)
		corr <- abs(corr)
	diag(corr) <- 0
	graph <- igraph::graph.adjacency(corr, mode="undirected", weighted = TRUE, diag = FALSE)
	return(graph)
}

.separate <- function (dtf, membership) {
	# separate a data frame and cormat according to the membership vector.
	ids <- sort(unique(membership))
	modules <- list()
	if (class(dtf) == "matrix") {
		rownames(dtf) <- NULL
		dtf <- data.frame(dtf)
	}
	
	for (i in 1:length(ids)) {
		curmod <- list()
		class(curmod) <- "cca.module"
		curmod$dtf <- dtf[membership == ids[i],]
		curmod$cormat <- cor(curmod$dtf)
		curmod$degenerate <- any(is.na(curmod$cormat))
		modules[[i]] <- curmod
	}
	
	return (modules)
}






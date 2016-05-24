#' @title Plots a flow chart with model comparisons
#'
#' @description
#' Plots a flow chart with model comparisons from a RSA object
#'
#' @details
#' The plot can be either requested within the \code{compare} function:
#' \code{compare(r1, plot=TRUE)}
#' Or it can be plotted from a cRSA object (= output from the \code{\link{compare}} function):
#' \code{c1 <- compare(r1)}
#' \code{plot(c1)}
#'
#' @export
#' @param x A cRSA object (= output from the \code{\link{compare}} function)
#' @param digits The number of digits to which numbers are rounded
#' @param sig Threshold for models to be marked as "not significant"
#' @param borderline Threshold for models to be marked as "borderline significant" (used for color of arrows)
#' @param ... Additional parameters (not used yet)
#'
#' @seealso \code{\link{RSA}}, \code{\link{compare}}
#'
#' @examples
#' \dontrun{
#' data(motcon)
#' r.m <- RSA(postVA~ePow*iPow, motcon)
#' c1 <- compare(r.m)
#' modeltree(c1)
#' }
modeltree <- function(x, digits=3, sig=.05, borderline=.10, ...) {
	if (!requireNamespace("qgraph", quietly = TRUE)) {
		stop('`qgraph` package needed for modeltrees. Please install with install.packages("qgraph")', call. = FALSE)
	}

	c1 <- x
	c1$fromto <- c("", paste0(c1$model[1:(nrow(c1)-1)], "_", c1$model[2:(nrow(c1))]))
	
	# define labels of boxes
	m <- c(
		"full", 	# 1
		"SRRR", 	# 2
		"IA", 		# 3
		"additive", # 4
		"diff", 	# 5
		"SRR", 		# 6
		"RR", 		# 7	
		"SQD", 	# 8
		"SRSQD", 	# 9
		"SSQD", 		#10
		"onlyx2", 	#11
		"onlyy2", 	#12
		"onlyx", 	#13
		"onlyy", 	#14
		"null"		#15
		)

	# get values from compare-object
	R2.adj <- c()
	CFI <- c()
	R2.p <- c()
	for (i in 1:length(m)) {
		R2.adj <- c(R2.adj, c1[c1$model==m[i], "R2.adj"][1])
		CFI <- c(CFI, c1[c1$model==m[i], "cfi"][1])
		R2.p <- c(R2.p, c1[c1$model==m[i], "R2.p"][1])
	}

	pos <- matrix(ncol=8, byrow=TRUE, c(
	16,  0,  0,  0,  1,  0,  0,  0,
	17,  0,  0,  0,  2,  0,  0,  0,	
	18,  0,  0,  3,  6,  9,  0,  0,
	19, 11,  0,  4,  7, 10,  0, 12,
	20, 13,  0,  5,  0,  8,  0, 14,
	21,  0,  0,  0, 15,  0,  0,  0
	))

	# define edgelist, without weight
	eL <- matrix(ncol=2, byrow=TRUE, data=c(
	1,2,
	1,3,
	1,11,
	1,12,
	3,4,
	4,5,
	5,15,
	2,6,
	6,7,
	7,8,
	8,15,
	2,9,
	9,10,
	10,8,
	6,10,
	11,13,
	12,14,
	13,15,
	14,15,
	
	16,16,
	17,17,
	18,18,
	19,19,
	20,20,
	21,21
	))

	# define weights of edges: will be translated to color
	w <- c()
	# 19 = number of edges defined in eL (upper block)
	for (i in 1:19) {
		#print(paste0(m[eL[i, 1]], "_", m[eL[i, 2]]))
	    w <- c(w, c1[c1$fromto == paste0(m[eL[i, 1]], "_", m[eL[i, 2]]), "Pr(>Chisq)"][1])
	}
	w[is.na(w)] <- 0

	# compute box labels
	lab <- c("full", "SRRR", "Interaction", "Additive", "Difference", "SRR", "RR", "Squared difference", "SRSQD", "SSQD", "x + x^2", "y + y^2", "x", "y", "Intercept only", "k = 5", "k = 4", "k = 3", "k = 2", "k = 1", "k = 0")	
	lab2 <- list()
	lab.color <- c()
	
	for (i in 1:length(lab)) {
		if (i <= length(m)) {
			lab2[[i]] <- paste(lab[i],
				paste("CFI = ", f2(CFI[i], digits)),
				paste("R^2[adj] = ", f2(R2.adj[i], digits)), sep="\n")
			lab.color <- c(lab.color, ifelse(R2.p < sig, "black", "grey"))	
		} else {
			lab2[[i]] <- lab[i]
			lab.color <- c(lab.color, "black")
		}
	}
	
	dev.new(width=11.5, height=6.5)
	p1 <- qgraph::qgraph(eL, 
		edgeList	= TRUE,
		nNodes		= length(lab),
		layout 		= pos, 
	
		# define edges
		edge.labels = c(paste0("p = ", f2(w, digits)), rep("", 6)), 
		edge.color	= c(pRamp(w, sig=sig, borderline=borderline), rep("#FFFFFF", 6)), 
		edge.label.cex = 0.6,
		edge.label.bg = "white",
		edge.label.position = 0.5,
		curve = c(0, -0.5, -1, 1, 0, 0, 0, 0,0,-0.5,0,0.5,0,0,0,0,0,-0.5,0.5, rep(0, 6)),	# hand-crafted edge curvature ...
		curveAll = TRUE,
	
		# define boxes
		labels		= lab2,
		label.color	= lab.color,
		label.prop	= 0.75,
		shape		= "rectangle",
		border.color = c(rep("black", 11), rep("white", 6), rep("black", 2)),
		asize		= 3,		# size of arrowhead
		vsize 		= 14,		# horizontal size of boxes
		vsize2 		= 7,		# vertical size of boxes
		border.width = R2.adj
	)
	invisible(p1)
}

#' @method plot cRSA
#' @export
plot.cRSA <- function(x, ...) {
	modeltree(x, ...)
}
#' Functions that create different example graphs
#' 
#' Functions that creates example graphs, e.g. graphs that represents a
#' Bonferroni-Holm adjustment, parallel gatekeeping or special procedures from
#' selected papers.
#' 
#' We are providing functions and not the resulting graphs directly because
#' this way you have additional examples: You can look at the function body
#' with \code{\link{body}} and see how the graph is built.
#' 
#' \describe{ \item{list("BonferroniHolm")}{ Returns a graph that represents a
#' Bonferroni-Holm adjustment.  The result is a complete graph, where all nodes
#' have the same weights and each edge weight is \eqn{\frac{1}{n-1}}{1/(n-1)}.
#' } \item{list("BretzEtAl2011")}{ Graph in figure 2 from Bretz et al. See
#' references (Bretz et al. 2011).  } \item{list("HommelEtAl2007")}{ Graph from
#' Hommel et al. See references (Hommel et al. 2007).  }
#' \item{list("parallelGatekeeping")}{ Graph for parallel gatekeeping. See
#' references (Dmitrienko et al. 2003).  }
#' \item{list("improvedParallelGatekeeping")}{ Graph for improved parallel
#' gatekeeping. See references (Hommel et al. 2007).  }
#' \item{list("HungEtWang2010")}{ Graph from Hung et Wang. See references (Hung
#' et Wang 2010).  } \item{list("MaurerEtAl1995")}{ Graph from Maurer et al.
#' See references (Maurer et al. 1995).  } \item{list("cycleGraph")}{ Cycle
#' graph. The weight \code{weights[i]} specifies the edge weight from node
#' \eqn{i}{i} to node \eqn{i+1}{i+1} for \eqn{i=1,\ldots,n-1}{i=1,...,n-1} and
#' \code{weight[n]} from node \eqn{n}{n} to node 1.  }
#' \item{list("improvedFallbackI")}{ Graph for the improved Fallback Procedure
#' by Wiens & Dmitrienko. See references (Wiens et Dmitrienko 2005).  }
#' \item{list("improvedFallbackII")}{ Graph for the improved Fallback Procedure
#' by Hommel & Bretz. See references (Hommel et Bretz 2008).  }
#' \item{list("Ferber2011")}{ Graph from Ferber et al. See references (Ferber
#' et al. 2011).  } \item{list("FerberTimeDose2011")}{ Graph from Ferber et al.
#' See references (Ferber et al. 2011).  } \item{list("Entangled1Maurer2012")}{
#' Entangled graph from Maurer et al. TODO: Add references as soon as they are
#' available.  } }
#' 
#' @name exampleGraphs
#' @aliases exampleGraphs BonferroniHolm BretzEtAl2011 parallelGatekeeping
#' improvedParallelGatekeeping HommelEtAl2007 HommelEtAl2007Simple
#' HungEtWang2010 MaurerEtAl1995 improvedFallbackI improvedFallbackII
#' cycleGraph fixedSequence generalSuccessive simpleSuccessiveI
#' simpleSuccessiveII truncatedHolm fallback HuqueAloshEtBhore2011
#' BauerEtAl2001 BretzEtAl2009a BretzEtAl2009b BretzEtAl2009c Ferber2011
#' Entangled1Maurer2012 Entangled2Maurer2012 FerberTimeDose2011 WangTing2014
#' @param n Number of hypotheses.
#' @param nodes Character vector of node names.
#' @param weights Numeric vector of node weights.
#' @param times Number of time points.
#' @param doses Number of dose levels.
#' @param w Further variable weight(s) in graph.
#' @param gamma An optional number in [0,1] specifying the value for variable gamma.
#' @param delta An optional number in [0,1] specifying the value for variable delta.
#' @param nu An optional number in [0,1] specifying the value for variable nu.
#' @param tau An optional number in [0,1] specifying the value for variable tau.
#' @param omega An optional number in [0,1] specifying the value for variable omega.
#' @return A graph of class \code{\link{graphMCP}} that represents a
#' sequentially rejective multiple test procedure.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @references Holm, S. (1979). A simple sequentally rejective multiple test
#' procedure. Scandinavian Journal of Statistics 6, 65-70.
#' 
#' Dmitrienko, A., Offen, W., Westfall, P.H. (2003). Gatekeeping strategies for
#' clinical trials that do not require all primary effects to be significant.
#' Statistics in Medicine. 22, 2387-2400.
#' 
#' Bretz, F., Maurer, W., Brannath, W., Posch, M.: A graphical approach to
#' sequentially rejective multiple test procedures. Statistics in Medicine 2009
#' vol. 28 issue 4 page 586-604.
#' \url{http://www.meduniwien.ac.at/fwf_adaptive/papers/bretz_2009_22.pdf}
#' 
#' Bretz, F., Maurer, W. and Hommel, G. (2011), Test and power considerations
#' for multiple endpoint analyses using sequentially rejective graphical
#' procedures. Statistics in Medicine, 30: 1489--1501.
#' 
#' Hommel, G., Bretz, F. und Maurer, W. (2007). Powerful short-cuts for
#' multiple testing procedures with special reference to gatekeeping
#' strategies. Statistics in Medicine, 26(22), 4063-4073.
#' 
#' Hommel, G., Bretz, F. (2008): Aesthetics and power considerations in
#' multiple testing - a contradiction? Biometrical Journal 50:657-666.
#' 
#' Hung H.M.J., Wang S.-J. (2010). Challenges to multiple testing in clinical
#' trials. Biometrical Journal 52, 747-756.
#' 
#' W. Maurer, L. Hothorn, W. Lehmacher: Multiple comparisons in drug clinical
#' trials and preclinical assays: a-priori ordered hypotheses. In Biometrie in
#' der chemisch-pharmazeutischen Industrie, Vollmar J (ed.). Fischer Verlag:
#' Stuttgart, 1995; 3-18.
#' 
#' Maurer, W., & Bretz, F. (2013). Memory and other properties of multiple test 
#' procedures generated by entangled graphs. Statistics in medicine, 32 (10), 1739-1753.
#' 
#' Wiens, B.L., Dmitrienko, A. (2005): The fallback procedure for evaluating a
#' single family of hypotheses. Journal of Biopharmaceutical Statistics
#' 15:929-942.
#' 
#' Wang, B., Ting, N. (2014). An Application of Graphical Approach to 
#' Construct Multiple Testing Procedures in a Hypothetical Phase III Design. 
#' Frontiers in public health, 1 (75).
#' 
#' Ferber, G. Staner, L. and Boeijinga, P. (2011): Structured multiplicity and
#' confirmatory statistical analyses in pharmacodynamic studies using the
#' quantitative electroencephalogram, Journal of neuroscience methods, Volume
#' 201, Issue 1, Pages 204-212.
#' @keywords misc graphs
#' @examples
#' 
#' g <- BonferroniHolm(5)
#' 
#' gMCP(g, pvalues=c(0.1, 0.2, 0.4, 0.4, 0.7))
#' 
#' HungEtWang2010()
#' HungEtWang2010(nu=1)
#' 
NULL

#' @rdname exampleGraphs
BonferroniHolm <- function(n, weights=rep(1/n, n)) {
	if (missing(n)) { stop("Please provide the number of hypotheses as parameter n.") }
	hnodes <- paste("H", 1:n, sep="")
	m <- matrix(1/(n-1), nrow=n, ncol=n)
	diag(m) <- 0
	rownames(m) <- colnames(m) <- hnodes
	BonferroniHolm <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- 100+(0:(n-1))*200
	nodeY <- rep(200, n)
	BonferroniHolm@nodeAttr$X <- nodeX
	BonferroniHolm@nodeAttr$Y <- nodeY		
	# Label settings
	for (i in 1:n) {
		n1 <- hnodes[i]
		for (j in (1:n)[-i]) {
			n2 <- hnodes[j]
			x <- ((i+j)*200-200)/2
			y <- 200 + ((i-j)*50)
			edgeAttr(BonferroniHolm, n1, n2, "labelX") <- x			
			edgeAttr(BonferroniHolm, n1, n2, "labelY") <- y
		}
	}
	attr(BonferroniHolm, "description") <- paste("Graph representing the (unweighted) Bonferroni-Holm-Procedure", 
			"",
			#"Most powerful test procedure (without further assumptions) that treats all hypotheses equally.",
			"The graph is a complete graph, where all nodes have the same weights and each edge weight is 1/(n-1).",
			"",
			"Literature: Holm, S. (1979). A simple sequentally rejective multiple test procedure. Scandinavian Journal of Statistics 6, 65-70.", sep="\n")
	return(BonferroniHolm)
}

#' @rdname exampleGraphs
BretzEtAl2011 <- function() {
	# M:
	m <- rbind(H11=c(0,   0.5, 0,   0.5, 0,   0  ),
			H21=c(1/3, 0,   1/3, 0,   1/3, 0  ),
			H31=c(0,   0.5, 0,   0,   0,   0.5),
			H12=c(0,   1,   0,   0,   0,   0  ),
			H22=c(0.5, 0,   0.5, 0,   0,   0  ),
			H32=c(0,   1,   0,   0,   0,   0  ))	
	# Graph creation
	weights <- c(1/3, 1/3, 1/3, 0, 0, 0)
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- rep(c(100, 300, 500), 2)
	nodeY <- rep(c(100, 300), each=3)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	# Label placement
	edgeAttr(graph, "H11", "H21", "labelX") <- 200
	edgeAttr(graph, "H11", "H21", "labelY") <- 80	
	edgeAttr(graph, "H31", "H21", "labelX") <- 400
	edgeAttr(graph, "H31", "H21", "labelY") <- 80	
	edgeAttr(graph, "H21", "H11", "labelX") <- 200
	edgeAttr(graph, "H21", "H11", "labelY") <- 120	
	edgeAttr(graph, "H21", "H31", "labelX") <- 400
	edgeAttr(graph, "H21", "H31", "labelY") <- 120	
	edgeAttr(graph, "H12", "H21", "labelX") <- 150
	edgeAttr(graph, "H12", "H21", "labelY") <- 250	
	edgeAttr(graph, "H22", "H11", "labelX") <- 250
	edgeAttr(graph, "H22", "H11", "labelY") <- 250	
	edgeAttr(graph, "H32", "H21", "labelX") <- 450
	edgeAttr(graph, "H32", "H21", "labelY") <- 250	
	edgeAttr(graph, "H22", "H31", "labelX") <- 350
	edgeAttr(graph, "H22", "H31", "labelY") <- 250	
	attr(graph, "description") <- paste("Graph representing the procedure from Bretz et al. (2011) - Figure 2", 
			"",
			"H11, H21 and H31 represent three primary hypotheses and H21, H22 and H23 the associated secondary hypotheses.",
			"",			
			"A secondary hypothesis is only tested if the associated primary hypotheses is rejected.",
			"",
			"Since in this example it is preferred to reject two adjacent hypotheses (like H11 and H21 instead of H11 and H31) there are only edges between adjacent nodes.",
			"",
			"Literature: Bretz, F., Maurer, W. and Hommel, G. (2011), Test and power considerations for multiple endpoint analyses using sequentially rejective graphical procedures. Statistics in Medicine, 30: 1489-1501.", sep="\n")
	return(graph)	
}

#' @rdname exampleGraphs
BauerEtAl2001 <- function() {
	# M:
	hnodes <- c("H_{E_1}","H_{E_2}","H_{E_3}","H_{S_1}","H_{S_2}","H_{S_3}")
	m <- rbind(c(0,    0,  0, 1, 0, 0),
			   c(0,    0,  0, 0, 1, 0),
			   c(0,    0,  0, 0, 0, 1),
			   c(0,  1/2,1/2, 0, 0, 0),
			   c(1/2,  0,1/2, 0, 0, 0),
			   c(1/2,1/2,  0, 0, 0, 0))
	rownames(m) <- colnames(m) <- hnodes
	# Graph creation
	weights <- c(1/3, 1/3, 1/3, 0, 0, 0)
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- rep(c(100, 300, 500), 2)
	nodeY <- rep(c(100, 300), each=3)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	# Label placement
	edgeAttr(graph, "H_{S_1}", "H_{E_3}", "labelX") <- 200
	edgeAttr(graph, "H_{S_1}", "H_{E_3}", "labelY") <- 100	
	edgeAttr(graph, "H_{S_3}", "H_{E_1}", "labelX") <- 400
	edgeAttr(graph, "H_{S_3}", "H_{E_1}", "labelY") <- 100	
	attr(graph, "description") <- paste("Graph representing the procedure from Bretz et al. (2011) - Figure 2", 
			"",
			"Literature: Bauer P., Brannath W., Posch M.: Multiple testing for identifying effective and safe treatments. Biometrical Journal 2001; 43:605-616.",
			"", 
			"Bretz F., Maurer W., Brannath W., Posch M.: A graphical approach to sequentially rejective multiple test procedures. Statistics in Medicine 2009; 28:586-604. Figure 8.", sep="\n")
	return(graph)	
}

#' @rdname exampleGraphs
BretzEtAl2009a <- function() {
	# M:
	hnodes <- c("H_{11}","H_{21}","H_{12}","H_{22}","H_{13}","H_{23}")
	m <- rbind(c(0,    0,  1, 0, 0, 0),
			c(0,    0,  0, 1, 0, 0),
			c(0,    0,  0, 0, 1, 0),
			c(0,    0,  0, 0, 0, 1),
			c(  0,  1,  0, 0, 0, 0),
			c(  1,  0,  0, 0, 0, 0))
	rownames(m) <- colnames(m) <- hnodes
	# Graph creation
	weights <- c(1/2, 1/2, 0, 0, 0, 0)
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings	
	graph <- placeNodes(graph, nrow=3, ncol=2)	
	# Label placement
	attr(graph, "description") <- paste("Graph from Bretz et al. (2009) - Figure 14", 
			"",
			"Literature: Bretz F., Maurer W., Brannath W., Posch M.: A graphical approach to sequentially rejective multiple test procedures. Statistics in Medicine 2009; 28:586-604. Figure 14.", sep="\n")
	return(graph)	
}

#' @rdname exampleGraphs
BretzEtAl2009b <- function() {
	# M:
	hnodes <- c("H_{11}","H_{21}","H_{12}","H_{22}","H_{13}","H_{23}")
	m <- rbind(c(0,    1/2,  1/2, 0, 0, 0),
			c(1/2,    0,  0, 1/2, 0, 0),
			c(0,    0,  0, 1/2, 1/2, 0),
			c(0,    0,  1/2, 0, 0, 1/2),
			c(  0,  1,  0, 0, 0, 0),
			c(  1,  0,  0, 0, 0, 0))
	rownames(m) <- colnames(m) <- hnodes
	# Graph creation
	weights <- c(1/2, 1/2, 0, 0, 0, 0)
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings	
	graph <- placeNodes(graph, nrow=3, ncol=2)	
	edgeAttr(graph, "H_{12}", "H_{22}", "labelX") <- 200
	edgeAttr(graph, "H_{12}", "H_{22}", "labelY") <- 260	
	edgeAttr(graph, "H_{22}", "H_{12}", "labelX") <- 200
	edgeAttr(graph, "H_{22}", "H_{12}", "labelY") <- 340	
	# Label placement
	attr(graph, "description") <- paste("Graph from Bretz et al. (2009) - Figure 14", 
			"",
			"Literature: Bretz F., Maurer W., Brannath W., Posch M.: A graphical approach to sequentially rejective multiple test procedures. Statistics in Medicine 2009; 28:586-604. Figure 14.", sep="\n")
	return(graph)	
}

#' @rdname exampleGraphs
BretzEtAl2009c <- function() {
	# M:
	hnodes <- c("H_{11}","H_{21}","H_{12}","H_{22}","H_{13}","H_{23}")
	m <- rbind(c(0,    1/2,  1/2, 0, 0, 0),
			   c(0,    0,  "1-\\epsilon", "\\epsilon", 0, 0),
			   c(0,    "1-\\epsilon",  0, 0, "\\epsilon", 0),
			   c(0,    0,  0, 0, "1-\\epsilon", "\\epsilon"),
			   c(  0,  0,  0, "1-\\epsilon", 0, "\\epsilon"),
			   c(  0,  0,  0, 0, 0, 0))
	rownames(m) <- colnames(m) <- hnodes
	# Graph creation
	weights <- c(1, 0, 0, 0, 0, 0)
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings	
	nodeX <- c(200, 300, 100, 300, 100, 200)
	nodeY <- c(100, 200, 200, 300, 300, 400)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	# Label placement
	attr(graph, "description") <- paste("Graph from Bretz et al. (2009) - Figure 15", 
			"",
			"Literature: Bretz F., Maurer W., Brannath W., Posch M.: A graphical approach to sequentially rejective multiple test procedures. Statistics in Medicine 2009; 28:586-604. Figure 15.", sep="\n")
	return(graph)	
}

#' @rdname exampleGraphs
HommelEtAl2007 <- function() {
	# Nodes:
	weights <- c(rep(1/3, 3), rep(0,4))	
	hnodes <- c("E1", "QoL", "E2", "D1", "D2", "D3", "D4")
	# Edges:
	m <- rbind(                 
		c("0",         "1", "0",         "0",                 "0",                 "0",                 "0"                ),
		c("0",         "0", "0",         "0.25",              "0.25",              "0.25",              "0.25"             ),
		c("0",         "1", "0",         "0",                 "0",                 "0",                 "0"                ),
		c("\\epsilon", "0", "\\epsilon", "0",                 "1/3-2/3*\\epsilon", "1/3-2/3*\\epsilon", "1/3-2/3*\\epsilon"),
		c("\\epsilon", "0", "\\epsilon", "1/3-2/3*\\epsilon", "0",                 "1/3-2/3*\\epsilon", "1/3-2/3*\\epsilon"),
		c("\\epsilon", "0", "\\epsilon", "1/3-2/3*\\epsilon", "1/3-2/3*\\epsilon", "0",                 "1/3-2/3*\\epsilon"),
		c("\\epsilon", "0", "\\epsilon", "1/3-2/3*\\epsilon", "1/3-2/3*\\epsilon", "1/3-2/3*\\epsilon", "0"))
	rownames(m) <- colnames(m) <- hnodes
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- c(200, 400, 600, 100, 300, 500, 700)
	nodeY <- c(100, 100, 100, 300, 300, 300, 300)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	for (i in 1:4) {
		n1 <- hnodes[3+i]
		for (j in (1:4)[-i]) {
			n2 <- hnodes[3+j]			
			x <- ((i+j)*200-200)/2+sign(i-j)*30
			y <- 300 + ((abs(i-j)-1)*60)+sign(i-j)*10+10
			edgeAttr(graph, n1, n2, "labelX") <- x
			edgeAttr(graph, n1, n2, "labelY") <- y
		}
	}
	attr(graph, "description") <- paste("Graph representing the procedure from Hommel et al. (2007)", 
			"",
			"In this clinical trial example three primary endpoints are investigated: QoL (Quality of Life score), E1 and E2.",
			"If QoL is rejected, four secondary hypotheses D1, D2, D3 and D4 are also be tested.",
			"",
			"Literature: Hommel, G., Bretz, F. und Maurer, W. (2007). Powerful short-cuts for multiple testing procedures with special reference to gatekeeping strategies. Statistics in Medicine, 26(22), 4063-4073.", sep="\n")
	attr(graph, "pvalues") <- c(0.097, 0.015, 0.005, 0.006, 0.004, 0.008, 0.04)
	return(graph)	
}

#' @rdname exampleGraphs
HommelEtAl2007Simple <- function() {
	# Nodes:
	weights <- c(rep(1/3, 3), rep(0,4))	
	# Edges:
	m <- structure(c(0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
					0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0.25, 
					0, 0, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0), .Dim = c(7L, 7L), .Dimnames = list(
					c("QoL", "E1", "E2", "D1", "D2", "D3", "D4"), c("QoL", "E1", 
							"E2", "D1", "D2", "D3", "D4"))) 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- c(200, 350, 450,  50, 150, 250, 350)
	nodeY <- c(150, 150, 150, 350, 350, 350, 350)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY
	
	edgeAttr(graph, "E2", "QoL", "labelX") <- 350
	edgeAttr(graph, "E2", "QoL", "labelY") <- 100
	
	attr(graph, "description") <- paste("Simplified graph representing the procedure from Hommel et al. (2007)", 
			"",
			"In this clinical trial example three primary endpoints are investigated: QoL (Quality of Life score), E1 and E2.",
			"If QoL is rejected, four secondary hypotheses D1, D2, D3 and D4 are also be tested.",
			"",
			"Literature: Hommel, G., Bretz, F. und Maurer, W. (2007). Powerful short-cuts for multiple testing procedures with special reference to gatekeeping strategies. Statistics in Medicine, 26(22), 4063-4073.", sep="\n")
	attr(graph, "pvalues") <- c(0.015, 0.097, 0.005, 0.006, 0.004, 0.008, 0.04)
	return(graph)	
}

#' @rdname exampleGraphs
parallelGatekeeping <- function() {
	# Nodes:
	weights <- rep(c(1/2,0), each=2)	
	hnodes <- paste("H", 1:4, sep="")
	# Edges:
	m <- rbind(
		c(0, 0, 0.5, 0.5),
		c(0, 0, 0.5, 0.5),
		c(0, 0, 0.0, 1.0),
		c(0, 0, 1.0, 0.0))
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- rep(c(100, 300), 2)
	nodeY <- rep(c(100, 300), each=2)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	# Label placement
	edgeAttr(graph, "H1", "H4", "labelX") <- 150
	edgeAttr(graph, "H1", "H4", "labelY") <- 150
	edgeAttr(graph, "H2", "H3", "labelX") <- 250
	edgeAttr(graph, "H2", "H3", "labelY") <- 150
	attr(graph, "description") <- paste("Graph representing a parallel gatekeeping procedure from Dmitrienko et al. (2003) Table I", 
			"",
			"Literature: Dmitrienko, A., Offen, W., Westfall, P.H. (2003). Gatekeeping strategies for clinical trials that do not require all primary effects to be significant. Statistics in Medicine. 22, 2387-2400.", sep="\n")
	return(graph)	
}

#' @rdname exampleGraphs
improvedParallelGatekeeping <- function() {
	graph <- parallelGatekeeping()
	graph <- setEdge("H3", "H1", graph, "\\epsilon")
	graph <- setEdge("H4", "H2", graph, "\\epsilon")
	graph <- setEdge("H3", "H4", graph, "1-\\epsilon")
	graph <- setEdge("H4", "H3", graph, "1-\\epsilon")
	edgeAttr(graph, "H1", "H3", "labelX") <- 100
	edgeAttr(graph, "H1", "H3", "labelY") <- 200
	edgeAttr(graph, "H2", "H4", "labelX") <- 300
	edgeAttr(graph, "H2", "H4", "labelY") <- 200
	edgeAttr(graph, "H3", "H1", "labelX") <-  70
	edgeAttr(graph, "H3", "H1", "labelY") <- 200
	edgeAttr(graph, "H4", "H2", "labelX") <- 330
	edgeAttr(graph, "H4", "H2", "labelY") <- 200
	attr(graph, "description") <- paste("Graph representing an improved parallel gatekeeping procedure", 
			"",
			"Literature: Bretz, F., Maurer, W., Brannath, W., Posch, M.: A graphical approach to sequentially rejective multiple test procedures. Statistics in Medicine 2009 vol. 28 issue 4 page 586-604. URL: http://www.meduniwien.ac.at/fwf_adaptive/papers/bretz_2009_22.pdf .", sep="\n")
	return(graph)	
}

#' @rdname exampleGraphs
fallback <- function(weights) {
	if (missing(weights)) { stop("Please provide weights.") }
	n <- length(weights)
	hnodes <- paste("H", 1:n, sep="")
	m <- matrix(0, nrow=n, ncol=n)
	for (i in 2:n) {
		m[i-1,i] <- 1
	}
	rownames(m) <- colnames(m) <- hnodes
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- 50+(0:(n-1))*150
	nodeY <- rep(100, n)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY		
	attr(graph, "description") <- paste("Graph representing the fallback (a fixed sequence Bonferroni) procedure", 
			"",
			"Literature: Wiens B.L.: A fixed sequence Bonferroni procedure for testing multiple endpoints. Pharmaceutical Statistics 2003; 2: 211-215.",
			"",
			"Bretz F., Maurer W., Brannath W., Posch M.: A graphical approach to sequentially rejective multiple test procedures. Statistics in Medicine 2009; 28:586-604.", sep="\n")
	return(graph)
}

#' @rdname exampleGraphs
fixedSequence <- function(n) {
	if (missing(n)) { stop("Please provide the number of hypotheses as parameter n.") }
	weights <- c(1, rep(0, n-1))
	graph <- fallback(weights)		
	attr(graph, "description") <- paste("Graph representing the fixed sequence test", 
			"",
			"Literature:  Maurer W., Hothorn L., Lehmacher W.: Multiple comparisons in drug clinical trials and preclinical assays: a-priori ordered hypotheses. In Biometrie in der chemisch-pharmazeutischen Industrie, Vollmar J (ed.). Fischer Verlag: Stuttgart, 1995; 3-18.",
            "",
			"Westfall P.H., Krishen A.: Optimally weighted, fixed sequence, and gatekeeping multiple testing procedures. Journal of Statistical Planning and Inference 2001; 99:25-40.", sep="\n")
	return(graph)
}

#' @rdname exampleGraphs
simpleSuccessiveI <- function() {
	graph <- generalSuccessive()
	graph <- replaceVariables(graph, variables=list("gamma"=0, "delta"=0))
	attr(graph, "description") <- paste("Simple successive graph from Maurer et al. (2011)", 
			"",
			"Literature: Maurer W., Glimm E., Bretz F.: Multiple and repeated testing of primary, co-primary and secondary hypotheses. Statistics in Biopharmaceutical Reserach 2011; (in press).",
			"",
			"F. Bretz, M. Posch, E. Glimm, F. Klinglmueller, W. Maurer, K. Rohmeyer (2011), Graphical approaches for multiple comparison procedures using weighted Bonferroni, Simes or parametric tests. To be published. Figure 12.", sep="\n")
	return(graph)
}

#' @rdname exampleGraphs
simpleSuccessiveII <- function() {
	graph <- generalSuccessive()
	graph <- replaceVariables(graph, variables=list("gamma"=1/2, "delta"=1/2))
	attr(graph, "description") <- paste("Simple successive graph from Maurer et al. (2011)", 
			"",
			"Literature: Maurer W., Glimm E., Bretz F.: Multiple and repeated testing of primary, co-primary and secondary hypotheses. Statistics in Biopharmaceutical Reserach 2011; (in press).", sep="\n")
	return(graph)
}

#' @rdname exampleGraphs
truncatedHolm <- function(gamma) {
	# Nodes:
	weights <- c(1/2, 1/2, 0, 0)
	hnodes <- paste("H", 1:4, sep="")
	# Edges:
	m <- rbind(
			c("0",      "\\gamma","(1-\\gamma)/2", "(1-\\gamma)/2"),
			c("\\gamma","0",      "(1-\\gamma)/2", "(1-\\gamma)/2"),
			c("0",      "0",      "0",             "1"),       
			c("0",      "0",      "1",             "0"))     
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- rep(c(100, 300), 2)
	nodeY <- rep(c(100, 300), each=2)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	attr(graph, "description") <- paste("Example of the Truncated Holm Procedure", 
			"",
			"Literature: Dmitrienko A, Tamhane A, Wiens B. General multi-stage gatekeeping procedures. Biometrical Journal 2008; 50:667-677.",
			"",
			"F. Bretz, M. Posch, E. Glimm, F. Klinglmueller, W. Maurer, K. Rohmeyer (2011), Graphical approaches for multiple comparison procedures using weighted Bonferroni, Simes or parametric tests. To be published. Figure 5.", sep="\n")
	if (!missing(gamma)) {
	  graph <- replaceVariables(graph, variables=list("gamma"=gamma))
	}
	return(graph)
}

# From Maurer, Glimm and Bretz 2011:
# A graph generates a successive procedure if
# - initially has weights 0 on all secondary hypotheses
# - the only edges with positive weight leading into a secondary
#   hypothesis are those originating at its parent primary hypotheses
# - and there are no edges leading from a secondary hypothesis to
#   another secondary hypothesis that has not the same parents.

#' @rdname exampleGraphs
generalSuccessive <- function(weights=c(1/2,1/2), gamma, delta) {
	if (length(weights)!=2) stop("Please specify the weights for H1 and H2 and only these.")
	# Nodes:
	weights <- c(weights, 0, 0)
	hnodes <- paste("H", 1:4, sep="")
	# Edges:
	m <- rbind(
		c("0",      "\\gamma","1-\\gamma","0"),
		c("\\delta","0",      "0",        "1-\\delta"),
		c("0",      "1",      "0",        "0"),       
		c("1",      "0",      "0",        "0"))     
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- rep(c(100, 300), 2)
	nodeY <- rep(c(100, 300), each=2)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	attr(graph, "description") <- paste("General successive graph from Bretz et al. (2011), Figure 6", 
			"",
			"Literature: Bretz, F., Maurer, W. and Hommel, G. (2011), Test and power considerations for multiple endpoint analyses using sequentially rejective graphical procedures. Statistics in Medicine, 30: 1489-1501.", sep="\n")
	variables <- list()
	if (!missing(gamma)) variables[["gamma"]] <- gamma
	if (!missing(delta)) variables[["delta"]] <- delta	
	if (length(variables)>0) {
	  graph <- replaceVariables(graph, variables=variables, partial=TRUE)
	}
	return(graph)		
}

#' @rdname exampleGraphs
HuqueAloshEtBhore2011 <- function() {
	graph <- HungEtWang2010()
	graph <- replaceVariables(graph, variables=list("nu"=1/2, "omega"=1/2, "tau"=0))
	# TODO: Create function for renaming of graph nodes in R
	rownames(graph@m) <- colnames(graph@m) <- paste("H", 1:4, sep="")
	names(graph@nodeAttr$rejected) <- rownames(graph@m)
	names(graph@weights) <- rownames(graph@m)
	# END of Todo
	graph@m["H4","H2"] <- 1
	attr(graph, "description") <- paste("Graph representing the procedure from Huque, Alosh and Bhore (2011)", 
			"",
			"Literature: Huque, M.F. and Alosh, M. and Bhore, R. (2011), Addressing Multiplicity Issues of a Composite Endpoint and Its Components in Clinical Trials. Journal of Biopharmaceutical Statistics, 21: 610-634.", sep="\n")
	return(graph)	
}

#' @rdname exampleGraphs
HungEtWang2010 <- function(nu, tau, omega) {
	# Nodes:
	weights <- c(1,0,0,0)	
	hnodes <- c("H_{1,NI}","H_{1,S}","H_{2,NI}","H_{2,S}")
	# Edges:
	m <- rbind(  
		c("0",     "\\nu",    "1-\\nu", "0"),        
		c("0",     "0",       "\\tau",  "1-\\tau"),  
		c("0",     "\\omega", "0",      "1-\\omega"),
		c("0",     "0",       "0",      "0"))
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	nodeX <- rep(c(100, 300), 2)
	nodeY <- rep(c(100, 300), each=2)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	attr(graph, "description") <- paste("Graph representing the procedure from Hung and Wang (2010)",
			"",
			"$H_{1,NI}$ : Non-inferiority of the primary endpoint",
			"$H_{1,S}$  : Superiority of the primary endpoint",
			"$H_{2,NI}$ : Non-inferiority of the secondary endpoint",
			"$H_{2,S}$  : Superiority of the secondary endpoint",
			"",
			"Literature: Hung H.M.J., Wang S.-J. (2010). Challenges to multiple testing in clinical trials. Biometrical Journal 52, 747-756.", sep="\n")
  variables <- list()
  if (!missing(nu)) variables[["nu"]] <- nu
	if (!missing(omega)) variables[["omega"]] <- omega
	if (!missing(tau)) variables[["tau"]] <- tau
  if (length(variables)>0) {
	  graph <- replaceVariables(graph, variables=variables, partial=TRUE)
  }
	return(graph)
}

#' @rdname exampleGraphs
MaurerEtAl1995 <- function() {
	# Nodes:
	weights <- c(1,0,0,0,0)	
	hnodes <- paste("H", 1:5, sep="")
	# Edges:
	m <- rbind( 
		c(0, 1, 0, 0.0, 0.0),
		c(0, 0, 1, 0.0, 0.0),
		c(0, 0, 0, 0.5, 0.5),
		c(0, 0, 0, 0.0, 0.0),
		c(0, 0, 0, 0.0, 0.0))
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	# Visualization settings
	nodeX <- c(100, 200, 300, 400, 400)
	nodeY <- c(100, 100, 100, 50, 150)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY	
	attr(graph, "description") <- paste("Graph representing a procedure in drug clinical trials (from Maurer et al. 1995, Scenario 1)",
			"",
			"In a univariate one-way design a drug A is compared against placebo and two positive control drugs B and C.",
			"",
			"The order of importance is that first the sensitivity has to be shown, i.e. that drug B and C are better than placebo. Than the efficacy of A vs. placebo is tested and if this can be shown, it is tested (with Bonferroni correction) whether A is superior to drug B and/or C.",
			"",
			"These hypotheses are represented in the graph as follows:",
			"H1: drug B better than placebo",
			"H2: drug C better than placebo",
			"H3: drug A better than placebo",
			"H4: drug A better than drug B",
			"H5: drug A better than drug C",
			"",
			"(Maurer et al. apply the intersection-union principle to H1 and H2 to test sensitivity, so sensitivity is shown if and only if H1 and H2 are both rejected.)",
			"",
			"Note that you could improve the test procedure by using a Bonferroni-Holm correction instead of the Bonferroni correction in the last step by adding an edge from H4 to H5 with weight 1 and vice versa.",
			"",
			"Literature:",
			"W. Maurer, L. Hothorn, W. Lehmacher: Multiple comparisons in drug clinical trials and preclinical assays: a-priori ordered hypotheses. In Biometrie in der chemisch-pharmazeutischen Industrie, Vollmar J (ed.). Fischer Verlag: Stuttgart, 1995; 3-18.", sep="\n")	
	return(graph)	
}

#' @rdname exampleGraphs
cycleGraph <- function(nodes, weights) {
	# Edges:
	n <- length(nodes)
	m <- diag(n)
	m <- rbind(m[2:n,],m[1,])
	rownames(m) <- colnames(m) <- nodes 	
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	return(graph)
}

#TODO @rdname exampleGraphs
gatekeeping <- function(n, type=c("serial", "parallel", "imporved parallel"), weights=rep(1/n, n)) {
	# Nodes:
	hnodes <- paste("H", 1:(2*n), sep="")
	# Edges:
	edges <- vector("list", length=4)
	for (i in 1:n) {
		
	}
	names(edges)<-hnodes
	# Graph creation
	graph <- new("graphMCP", nodes=hnodes, edgeL=edges, weights=c(weights, rep(0, n)))
	# Visualization settings
	nodeX <- c(100, 200, 300, 400, 400)
	nodeY <- c(100, 100, 100, 50, 150)
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY
	attr(graph, "description") <- paste("Graph representing ...",
			"",
			"Literature:", sep="\n")	
	return(graph)
}

#' @rdname exampleGraphs
improvedFallbackI <- function(weights=rep(1/3, 3)) {
	# Nodes:
	hnodes <- paste("H", 1:3, sep="")
	# Edges:
	m <- rbind( 
		c(0.0, 1.0,  0),
		c(0.0, 0.0,  1),
		c(0.5, 0.5,  0))
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	graph <- placeNodes(graph, nrow=1, ncol=3)
	attr(graph, "description") <- paste("Improved Fallback Method I by Wiens & Dmitrienko",
			"",
			"Literature: B.L. Wiens, A. Dmitrienko (2005): The fallback procedure for evaluating a single family of hypotheses. Journal of Biopharmaceutical Statistics 15:929-942.",
			"",
			"Bretz, F., Maurer, W. and Hommel, G. (2011), Test and power considerations for multiple endpoint analyses using sequentially rejective graphical procedures. Statistics in Medicine, 30: 1489-1501.", sep="\n")
	edgeAttr(graph, "H3", "H1", "labelX") <- 300
	edgeAttr(graph, "H3", "H1", "labelY") <- 200
	edgeAttr(graph, "H2", "H3", "labelX") <- 400
	edgeAttr(graph, "H2", "H3", "labelY") <- 100
	edgeAttr(graph, "H3", "H2", "labelX") <- 400
	edgeAttr(graph, "H3", "H2", "labelY") <- 135	
	return(graph)
} 

#' @rdname exampleGraphs
improvedFallbackII <- function(weights=rep(1/3, 3)) {
	# Nodes:
	hnodes <- paste("H", 1:3, sep="")
	# Edges:
	m <- rbind(         
		c("0",           "1", "0"        ),
		c("1-\\epsilon", "0", "\\epsilon"),
		c("1",           "0", "0"        ))
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=weights)
	graph <- placeNodes(graph, nrow=1, ncol=3)
	attr(graph, "description") <- paste("Improved Fallback Method II by Hommel & Bretz",
			"",
			"Literature: G. Hommel, F. Bretz (2008): Aesthetics and power considerations in multiple testing - a contradiction? Biometrical Journal 50:657-666.",
			"",
			"Bretz, F., Maurer, W. and Hommel, G. (2011), Test and power considerations for multiple endpoint analyses using sequentially rejective graphical procedures. Statistics in Medicine, 30: 1489-1501.", sep="\n")
	edgeAttr(graph, "H3", "H1", "labelX") <- 300
	edgeAttr(graph, "H3", "H1", "labelY") <- 200
	edgeAttr(graph, "H1", "H2", "labelX") <- 200
	edgeAttr(graph, "H1", "H2", "labelY") <- 100
	edgeAttr(graph, "H2", "H1", "labelX") <- 200
	edgeAttr(graph, "H2", "H1", "labelY") <- 135	
	return(graph)
} 

#' @rdname exampleGraphs
FerberTimeDose2011 <- function(times, doses, w="\\nu") {
	# Nodes:
	hnodes <- paste(rep(paste("T", 1:times, sep=""), each=doses),"D",1:doses, sep="")
	if (times<2) stop("times has to be an integer > 1")
	if (doses<2) stop("doses has to be an integer > 1")
	# Edges:
	m <- matrix(0, times*doses, times*doses)
	w2 <- paste("1-", w, sep="")
	for (i in 2:(doses)) {
		m[i,i-1] <- w	
		m[doses*(times-1)+i, doses*(times-1)+i-1] <- w2
		for (j in 2:(times)) {			
			m[doses*(j-1)+(i-1),doses*(j-2)+i] <- w		
			m[doses*(j-2)+i, doses*(j-1)+(i-1)] <- w2
		}
	}
	for (j in 2:(times)) {
		m[doses*(j-1)+1,doses*(j-2)+1] <- w2
		m[doses*(j), doses*(j-1)] <- w
	}
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=c(rep(0, times*doses-1), 1))
	attr(graph, "description") <- paste("Second graph from Ferber et al. 2011",
			"",
			"Literature: G. Ferber, L. Staner and P. Boeijinga (2011): Structured multiplicity and confirmatory statistical analyses in pharmacodynamic studies using the quantitative electroencephalogram, Journal of neuroscience methods, Volume 201, Issue 1, Pages 204-212.", sep="\n")
	# Placing nodes and converting to numeric
	graph <- placeNodes(parse2numeric(graph), times, doses)
	return(graph)
}

#' @rdname exampleGraphs
Ferber2011 <- function(w) {
	# Nodes:
	hnodes <- c("\\delta", "\\theta", "\\beta", "\\alpha", "\\alpha_1", 
			"\\alpha_2", "\\beta_1", "\\beta_2", "\\beta_3")
	# Edges:
	m <- matrix(c("0", "1/3", "1/6", "1/6", "0", "0", "0", "0", "0", 
					"1/3", "0", "1/6", "1/6", "0", "0", "0", "0", "0", 
					"1/3", "1/3", "0", "1/6", "0", "0", "1-w", "1-w", "1-w", 
					"1/3", "1/3", "1/6", "0", "1-w", "1-w", "0", "0", "0", 
					"0", "0", "0", "0.25", "0", "w", "0", "0", "0", 
					"0", "0", "0", "0.25", "w", "0", "0", "0", "0", 
					"0", "0", "1/6", "0", "0", "0", "0", "w/2", "w/2", 
					"0", "0", "1/6", "0", "0", "0", "w/2", "0", "w/2", 
					"0", "0", "1/6", "0", "0", "0", "w/2", "w/2", "0"), nrow=9)
	rownames(m) <- colnames(m) <- hnodes 
	# Graph creation
	graph <- new("graphMCP", m=m, weights=c(0.2,0.2,0.3,0.3,0,0,0,0,0))
	attr(graph, "description") <- paste("Graph from Ferber et al. 2011",
			"",
			"Literature: G. Ferber, L. Staner and P. Boeijinga (2011): Structured multiplicity and confirmatory statistical analyses in pharmacodynamic studies using the quantitative electroencephalogram, Journal of neuroscience methods, Volume 201, Issue 1, Pages 204-212.", sep="\n")
	
	graph@nodeAttr$X <- c(80, 350, 350, 80, 0, 150, 230, 480, 350)
	graph@nodeAttr$Y <- c(30, 30, 180, 180, 300, 300, 300, 300, 430)
	
	graph@edgeAttr$labelX <- structure(c(NA, 242, 296, 62, NA, NA, NA, NA, NA, 196, NA, 333, 
					130, NA, NA, NA, NA, NA, 145, 366, NA, 172, NA, NA, 265, 443, 
					323, 94, 303, 268, NA, 10, 146, NA, NA, NA, NA, NA, NA, 64, NA, 
					78, NA, NA, NA, NA, NA, NA, 96, 72, NA, NA, NA, NA, NA, NA, 314, 
					NA, NA, NA, NA, 354, 263, NA, NA, 392, NA, NA, NA, 348, NA, 404, 
					NA, NA, 377, NA, NA, NA, 291, 439, NA), .Dim = c(9L, 9L), .Dimnames = list(
					c("\\delta", "\\theta", "\\beta", "\\alpha", "\\alpha_1", 
							"\\alpha_2", "\\beta_1", "\\beta_2", "\\beta_3"), c("\\delta", 
							"\\theta", "\\beta", "\\alpha", "\\alpha_1", "\\alpha_2", 
							"\\beta_1", "\\beta_2", "\\beta_3")))
	graph@edgeAttr$labelY <- structure(c(NA, 39, 131, 120, NA, NA, NA, NA, NA, 10, NA, 121, 
					134, NA, NA, NA, NA, NA, 77, 89, NA, 171, NA, NA, 231, 233, 301, 
					89, 70, 192, NA, 252, 248, NA, NA, NA, NA, NA, NA, 233, NA, 319, 
					NA, NA, NA, NA, NA, NA, 233, 281, NA, NA, NA, NA, NA, NA, 234, 
					NA, NA, NA, NA, 339, 393, NA, NA, 233, NA, NA, NA, 265, NA, 370, 
					NA, NA, 299, NA, NA, NA, 365, 384, NA), .Dim = c(9L, 9L), .Dimnames = list(
					c("\\delta", "\\theta", "\\beta", "\\alpha", "\\alpha_1", 
							"\\alpha_2", "\\beta_1", "\\beta_2", "\\beta_3"), c("\\delta", 
							"\\theta", "\\beta", "\\alpha", "\\alpha_1", "\\alpha_2", 
							"\\beta_1", "\\beta_2", "\\beta_3")))
	
	if (!missing(w)) {
	  graph <- replaceVariables(graph, variables=list("w"=w))
	}
  
	return(graph)
}

#' @rdname exampleGraphs
Entangled1Maurer2012 <- function() {
	m <- rbind(H1=c(0, 0, 1, 0, 0),
			H2=c(0, 0, 1, 0, 0),
			H3=c(0, 0, 0, "1-\\epsilon", "\\epsilon"),
			H4=c(0, 1, 0, 0, 0),
			H5=c(0, 0, 0, 0, 0))
	
	weights <- c(1, 0, 0, 0, 0)
	
	graph1 <- new("graphMCP", m=m, weights=weights)
	
	graph1@nodeAttr$X <- c(100, 300, 200, 100, 300)
	graph1@nodeAttr$Y <- c(100, 100, 200, 300, 300)
	
	edgeAttr(graph1, "H4", "H2", "labelX") <- 50
	edgeAttr(graph1, "H4", "H2", "labelY") <- 50
	
	edgeAttr(graph1, "H1", "H3", "labelX") <- 170
	edgeAttr(graph1, "H1", "H3", "labelY") <- 130
	
	edgeAttr(graph1, "H2", "H3", "labelX") <- 270
	edgeAttr(graph1, "H2", "H3", "labelY") <- 170
	
	edgeAttr(graph1, "H3", "H4", "labelX") <- 170
	edgeAttr(graph1, "H3", "H4", "labelY") <- 270
	
	edgeAttr(graph1, "H3", "H5", "labelX") <- 270
	edgeAttr(graph1, "H3", "H5", "labelY") <- 230
	
	m <- rbind(H1=c(0, 0, 1, 0, 0),
			H2=c(0, 0, 1, 0, 0),
			H3=c(0, 0, 0, "\\epsilon", "1-\\epsilon"),
			H4=c(0, 0, 0, 0, 0),
			H5=c(1, 0, 0, 0, 0))
	
	weights <- c(0, 1, 0, 0, 0)
	
	graph2 <- new("graphMCP", m=m, weights=weights)
	
	edgeAttr(graph2, "H5", "H1", "labelX") <- 350
	edgeAttr(graph2, "H5", "H1", "labelY") <- 50
	
	edgeAttr(graph2, "H1", "H3", "labelX") <- 130
	edgeAttr(graph2, "H1", "H3", "labelY") <- 170
	
	edgeAttr(graph2, "H2", "H3", "labelX") <- 230
	edgeAttr(graph2, "H2", "H3", "labelY") <- 130
	
	edgeAttr(graph2, "H3", "H4", "labelX") <- 130
	edgeAttr(graph2, "H3", "H4", "labelY") <- 230
	
	edgeAttr(graph2, "H3", "H5", "labelX") <- 230
	edgeAttr(graph2, "H3", "H5", "labelY") <- 270
	
	graph <- new("entangledMCP", subgraphs=list(graph1,graph2), weights=c(0.5,0.5))

	attr(graph, "description") <- paste("Graph from Maurer and Bretz 2012",
			"",
			"Literature: Maurer, W., & Bretz, F. (2013). Memory and other properties of multiple test procedures generated by entangled graphs. Statistics in medicine, 32(10), 1739-1753.", sep="\n")
	
	return(graph)
}

#' @rdname exampleGraphs
Entangled2Maurer2012 <- function() {
	m <- rbind(H1=c(0, 0, 1, 0, 0),
			H2=c(0, 0, 0, 0, 1),
			H3=c(0, 0, 0, 1, 0),
			H4=c(0, 1, 0, 0, 0),
			H5=c(0, 0, 0, 0, 0))
	
	weights <- c(1, 0, 0, 0, 0)
	
	graph1 <- new("graphMCP", m=m, weights=weights)
	
	graph1@nodeAttr$X <- c(100, 300, 200, 100, 300)
	graph1@nodeAttr$Y <- c(100, 100, 200, 300, 300)
	
	edgeAttr(graph1, "H4", "H2", "labelX") <- 50
	edgeAttr(graph1, "H4", "H2", "labelY") <- 50
	
	m <- rbind(H1=c(0, 0, 0, 1, 0),
			H2=c(0, 0, 1, 0, 0),
			H3=c(0, 0, 0, 0, 1),
			H4=c(0, 0, 0, 0, 0),
			H5=c(1, 0, 0, 0, 0))
	
	weights <- c(0, 1, 0, 0, 0)
	
	graph2 <- new("graphMCP", m=m, weights=weights)
	
	edgeAttr(graph2, "H5", "H1", "labelX") <- 350
	edgeAttr(graph2, "H5", "H1", "labelY") <- 50
	
	graph <- new("entangledMCP", subgraphs=list(graph1,graph2), weights=c(0.5,0.5))
	
	attr(graph, "description") <- paste("Graph from Maurer and Bretz 2012",
			"",
			"Literature: Maurer, W., & Bretz, F. (2013). Memory and other properties of multiple test procedures generated by entangled graphs. Statistics in medicine, 32(10), 1739-1753.", sep="\n")
	
	return(graph)
}

#' @rdname exampleGraphs
WangTing2014 <- function(nu, tau) {
  m <- rbind(H1=c("0", "1-\\nu", "\\nu/2", "\\nu/2", "0", "0"),
             H2=c("1-\\nu", "0", "0", "0", "\\nu/2", "\\nu/2"),
             H3=c("0", "\\tau", "0", "1-\\tau", "0", "0"),
             H4=c("0", "\\tau", "1-\\tau", "0", "0", "0"),
             H5=c("\\tau", "0", "0", "0", "0", "1-\\tau"),
             H6=c("\\tau", "0", "0", "0", "1-\\tau", "0"))
  
  weights <- c(1/2, 1/2, 0, 0, 0, 0)
  
  graph <- new("graphMCP", m=m, weights=weights)
  
  graph@nodeAttr$X <- c(350, 350, 100, 250, 450, 600)-50
  graph@nodeAttr$Y <- c(100, 400, 250, 250, 250, 250)-50
  
  edgeAttr(graph, "H1", "H3", "labelX") <- 125
  edgeAttr(graph, "H1", "H3", "labelY") <- 125
  edgeAttr(graph, "H1", "H4", "labelX") <- 225
  edgeAttr(graph, "H1", "H4", "labelY") <- 125
  edgeAttr(graph, "H3", "H2", "labelX") <- 125
  edgeAttr(graph, "H3", "H2", "labelY") <- 275
  edgeAttr(graph, "H4", "H2", "labelX") <- 225
  edgeAttr(graph, "H4", "H2", "labelY") <- 275
  
  edgeAttr(graph, "H6", "H1", "labelX") <- 550-125+50
  edgeAttr(graph, "H6", "H1", "labelY") <- 125
  edgeAttr(graph, "H5", "H1", "labelX") <- 550-225+50
  edgeAttr(graph, "H5", "H1", "labelY") <- 125
  edgeAttr(graph, "H2", "H6", "labelX") <- 550-125+50
  edgeAttr(graph, "H2", "H6", "labelY") <- 275
  edgeAttr(graph, "H2", "H5", "labelX") <- 550-225+50
  edgeAttr(graph, "H2", "H5", "labelY") <- 275
  
  attr(graph, "description") <- paste("Graph from Wang and Ting 2014",
                                      "With \\nu=1, \\tau=0.5 graph from figure 3,",
                                      "with \\nu=1, \\tau=\\epsilon graph from figure 4 and",
                                      "with \\nu=\\epsilon, \\tau=\\epsilon graph from figure 5.",
                                      "",
                                      "Literature: Wang, B., Ting, N. (2014). An Application of Graphical Approach to Construct Multiple Testing Procedures in a Hypothetical Phase III Design. Frontiers in public health, 1 (75).",
                                      "URL: http://journal.frontiersin.org/Journal/10.3389/fpubh.2013.00075/full", sep="\n")
  
  variables <- list()
  if (!missing(nu)) variables[["nu"]] <- nu
  if (!missing(tau)) variables[["tau"]] <- tau
  if (length(variables)>0) {
    graph <- replaceVariables(graph, variables=variables, partial=TRUE)
  }
  
  return(graph)
}
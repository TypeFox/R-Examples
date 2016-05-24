updateGraphToNewClassDefinition <- function(object, askToInstallGraph=TRUE) {
	if ("entangledMCP" %in% class(object)) return(object)
	if (!("try-error"%in% class(try(object@m, silent=TRUE)))) return(object)
	if (!requireLibrary("graph")) {
		stop(paste("Library graph is required for converting but was not installed.",
						"Please try to install package 'graph' manually.",
						"See http://www.bioconductor.org/packages/release/bioc/html/graph.html"
				, sep="\n"))
	}	
	nodes <- object@nodes
	edges <- object@edgeL
	eData <- object@edgeData@data
	nData <- object@nodeData@data
	renderInfo <- object@renderInfo
	m <- matrix(0, nrow=length(nodes), ncol=length(nodes))
	edgeAttr <- list()
	nodeAttr <- list()
	w <- c()
	rejected <- c()
	for (n in nodes) {
		w <- c(w, nData[[n]]$nodeWeight)
		rejected <- c(rejected, (!is.null(nData[[n]]$rejected)) && nData[[n]]$rejected)
	}
	names (w) <- names(rejected) <- nodes
	nodeAttr$rejected <- rejected
	nodeAttr$X <- renderInfo@nodes$nodeX
	names(nodeAttr$X) <- nodes
	nodeAttr$Y <- renderInfo@nodes$nodeY
	names(nodeAttr$Y) <- nodes
	for (i in 1:length(edges)) {
		indices <- edges[i][[1]]$edges	
		for(j in indices) {
			d <- eData[paste(nodes[i], nodes[j], sep="|")][[1]]
			m[i,j] <- d$weight			
			if (is.nan(d$weight)) {				
				if (!is.null(d$variableWeight)) {
					m[i,j] <- d$variableWeight
				}
			}
			if (!is.null(d$epsilon)) {
				m[i,j] <- makeEpsilonString(d$weight, d$epsilon)	
			}
			if (is.null(edgeAttr[[nodes[i]]])) edgeAttr[[nodes[i]]] <- list()
			if (!any(is.null(c(d$labelX, d$labelY)))) {
				if (is.null(edgeAttr[[nodes[i]]][[nodes[j]]])) edgeAttr[[nodes[i]]][[nodes[j]]] <- list()
				edgeAttr[[nodes[i]]][[nodes[j]]]$coordinates <- c(d$labelX, d$labelY)
			}
		}		
	}
	rownames(m) <- colnames(m) <- nodes
	return(new("graphMCP", m=m, weights=w, edgeAttr=edgeAttr, nodeAttr=nodeAttr))
}

# The function 'makeEpsilonString' takes the real part of the weight and 
# a coefficent vector 'p' for different powers of epsilon and returns the combined weight as a character vector.
# Examples:
#  makeEpsilonString(1/2,c(1,-1,1/2)) == "1/2+\\epsilon-\\epsilon^2+1/2*\\epsilon^3"
#  makeEpsilonString(0,c(1,-1,1/2)) == "\\epsilon-\\epsilon^2+1/2*\\epsilon^3"
#  makeEpsilonString(1/2,c()) == "1/2"

makeEpsilonString <- function(weight, p) {	
	frac <- function(x) {as.character(fractions(x))}
	e <- "\\epsilon"	
	pStr <- ""
	for (i in 1:length(p)) {
		if (length(p)!=0 && !isTRUE(all.equal(p[i], 0))) {
			if (p[i]>=0) {
				pStr <- paste(pStr, "+", sep="")
			}
			if (!isTRUE(all.equal(p[i], 1))) {
				if (!isTRUE(all.equal(p[i], -1))) {
					pStr <- paste(pStr, frac(p[i]), "*", e, sep="")
				} else {
					pStr <- paste(pStr, "-", e, sep="")
				}
			} else {
				pStr <- paste(pStr, e, sep="")
			}
			if (i>1) {
				pStr <- paste(pStr, "^", i, sep="")
			}
		}
	}
	if (weight==0 && pStr!="") { # Remove the first "+" and just return the epsilon part:
		return(substring(pStr, 2))
	}
	return(paste(frac(weight), pStr, sep=""))	
}
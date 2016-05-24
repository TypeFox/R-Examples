getEdges <- function(graph){
	fromL <- c()
	toL <- c()
	weightL <- numeric(0)
	labelx <- numeric(0)
	labely <- numeric(0)
	curveL <- logical(0)
	weightStrL <- character(0)
	for (from in getNodes(graph)) {
		for (to in getNodes(graph)) {
			options(warn=-1)
			if (is.na(as.numeric(graph@m[from, to])) || as.numeric(graph@m[from, to])!=0) {	
				options(warn=0)			
				x <- try(unlist(edgeAttr(graph, from, to, "labelX")), silent = TRUE)
				if (class(x)!="try-error" && !is.null(x) && !is.na(x)) {
					labelx <- c(labelx, x)
				} else {
					labelx <- c(labelx, -100)
				}
				y <- try(unlist(edgeAttr(graph, from, to, "labelY")), silent = TRUE)
				if (class(y)!="try-error" && !is.null(y) && !is.na(y)) {
					labely <- c(labely, y)
				} else {
					labely <- c(labely, -100)
				}
				
				fromL <- c(fromL, from)
				toL <- c(toL, to)
				options(warn=-1)
				weightL <- c(weightL, as.numeric(graph@m[from, to]))				
				curve <- (is.na(as.numeric(graph@m[to, from])) || as.numeric(graph@m[to, from])!=0)
				options(warn=0)
				curveL <- c(curveL, curve)
				weightStrL <- c(weightStrL, graph@m[from, to])
			}
		}
	}
	return(list(from=fromL, to=toL, weight=weightL, labelx=labelx, labely=labely, curve=curveL, weightStr=as.character(weightStrL)))
}

checkPSD <- function(m) {
	if (!all(!is.na(m))) {
		return("Eigen values couldn't be checked, since matrix contains NAs.")	
	}
	ev <- eigen(m, symmetric=TRUE, only.values=TRUE)
	# We use the same check as mvtnorm to minimize problems:
	if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
		return(paste("Matrix is not positive semidefinite, lowest eigen value is:",min(ev$values)))
	}
	return("")
}

# Checks the following properties:
# Values must be between -1 and 1.
# Diagonal must be equal to 1.
# Matrix must be symmetric.
checkCorrelation <- function(m, returnMessage=FALSE, na.allowed=TRUE) {
  if (!na.allowed && any(is.na(m))) {
    if (returnMessage) return("Matrix can not contain NAs.")
    return(FALSE)
  }
  if (!is.numeric(m) || !is.matrix(m)) {
    if (returnMessage) return("Matrix must be a numeric matirx.")
    return(FALSE)
  }
  if (!isTRUE(all.equal(1, max(1, max(abs(m)[!is.na(m)]))))) {
    if (returnMessage) return("Values must be between -1 and 1.")
    return(FALSE)
  }
  if (!isTRUE(all.equal(diag(m), rep(1, dim(m)[1])))) {
    if (returnMessage) return("Diagonal must be equal to 1.")
    return(FALSE)
  }  
  if (!isSymmetric(unname(m))) {
    if (returnMessage) return("Matrix must be symmetric.")
    return(FALSE)
  }
  return(TRUE)
}

checkQuadraticMatrix <- function(m, n=NULL) {
  if (!is.numeric(m) || !is.matrix(m)) return(FALSE)
  if (!is.null(n) && dim(m)[1]!=n) return(FALSE)
  return(dim(m)[1]==dim(m)[2])
}

forceCorrelation <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(cov2cor(m))
}

#' Placement of graph nodes
#' 
#' Places the nodes of a graph according to a specified layout.
#' 
#' If one of \code{nrow} or \code{ncol} is not given, an attempt is made to
#' infer it from the number of nodes of the \code{graph} and the other
#' parameter.  If neither is given, the graph is placed as a circle.
#' 
#' @param graph A graph of class \code{\link{graphMCP}} or class
#' \code{\link{entangledMCP}}.
#' @param nrow The desired number of rows.
#' @param ncol The desired number of columns.
#' @param byrow Logical whether the graph is filled by rows (otherwise by
#' columns).
#' @param topdown Logical whether the rows are filled top-down or bottom-up.
#' @param force Logical whether a graph that has already a layout should be
#' given the specified new layout.
#' @return The graph with nodes placed according to the specified layout.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{graphMCP}}, \code{\link{entangledMCP}}
#' @keywords graphs
#' @examples
#' 
#' 
#' g <- matrix2graph(matrix(0, nrow=6, ncol=6))
#' 
#' g <- placeNodes(g, nrow=2, force=TRUE)
#' 
#' \dontrun{
#' graphGUI(g)
#' 
#' }
#' 
#' 
#' @export placeNodes

placeNodes <- function(graph, nrow, ncol, byrow = TRUE, topdown = TRUE, force = FALSE) {
	entangled <- NULL
	# If the graph is entangled only place the nodes of the first graph
	if ("entangledMCP" %in% class(graph)) {
		entangled <- graph
		graph <- entangled@subgraphs[[1]]
	}
	# Only place nodes if  no placement data exists or parameter force is set to TRUE
	if (is.null(graph@nodeAttr$X) || force) {
		n <- length(getNodes(graph))
		if (missing(nrow) && missing(ncol)) {		
			v <- (1:n)/n*2*pi
			nodeX <- 300 + 250*sin(v)
			nodeY <- 300 + 250*cos(v)			
		} else {
			if (missing(nrow)) {
				nrow <- ceiling(length(getNodes(graph))/ncol)
			}
			if (missing(ncol)) {
				ncol <- ceiling(length(getNodes(graph))/nrow)
			}
			# if nrow*ncol<n increase the number of rows
			if (nrow*ncol<n) {
				nrow <- ceiling(length(getNodes(graph))/ncol)
			}
			if (byrow) {
				nodeX <- rep(((1:ncol)-1)*200+100, nrow)
				nodeY <- rep(((1:nrow)-1)*200+100, each = ncol)
			} else {
				nodeX <- rep(((1:ncol)-1)*200+100, each = nrow)
				nodeY <- rep(((1:nrow)-1)*200+100, ncol)
			}
			if (!topdown) nodeY <- max(nodeY) - nodeY + 100
		}
		graph@nodeAttr$X <- nodeX[1:n]
		graph@nodeAttr$Y <- nodeY[1:n]
		for (i in getNodes(graph)) {
			for (j in getNodes(graph)) {
				if (graph@m[i,j]!=0) {
					edgeAttr(graph, i, j, "labelX") <- -100
					edgeAttr(graph, i, j, "labelY") <- -100
				}
			}
		}		
	}	
	if (!is.null(entangled)) {
		entangled@subgraphs[[1]] <- graph
		return(entangled)
	}
	return(graph)	
}

# I guess I simply don't understand how the graph package is supposed to be used.
# Or they have a bug. I have contacted them but got no response. Therefore I still use this stupid work-around:
stupidWorkAround <- function(graph) {
	if (length(graph@edgeAttr@data)>0) {
		for (i in 1:length(graph@edgeAttr@data)) {
			if (length(graph@edgeAttr@data[[i]])>0) {
				for (j in 1:length(graph@edgeAttr@data[[i]])){
					graph@edgeAttr@data[[i]][[j]] <- unname(graph@edgeAttr@data[[i]][[j]])
				}
			}
		}
	}
	return(graph)
}

getAllCorrelationMatrices <- function(envir=globalenv(), n="all") {
  objects <- ls(envir)
  matrices <- c()
  for (obj in objects) {
    candidate <- get(obj, envir=envir)
    if (is.matrix(candidate) && dim(candidate)[1] == dim(candidate)[2] && checkCorrelation(candidate)) {
      if (n=="all" || dim(candidate)[1]==n) {
        matrices <- c(matrices, obj)
      }
    }
  }
  if (length(matrices)==0) {
    if (n=="all") {
      return("No quadratic matrices found.")
    } else {
      return(paste("No ",n,"x",n,"-matrices found.", sep=""))
    }
  }
  return(matrices)
}

getAllQuadraticMatrices <- function(envir=globalenv(), n="all") {
	objects <- ls(envir)
	matrices <- c()
	for (obj in objects) {
		candidate <- get(obj, envir=envir)
		if (is.matrix(candidate) && dim(candidate)[1] == dim(candidate)[2]) {
			if (n=="all" || dim(candidate)[1]==n) {
				matrices <- c(matrices, obj)
			}
		}
	}
	if (length(matrices)==0) {
		if (n=="all") {
			return("No quadratic matrices found.")
		} else {
			return(paste("No ",n,"x",n,"-matrices found.", sep=""))
		}
	}
	return(matrices)
}

getAllGraphs <- function(envir=globalenv()) {
	objects <- ls(envir)
	graphs <- c()
	for (obj in objects) {
		candidate <- get(obj, envir=envir)
		if ("graphMCP" %in% class(candidate) || "entangledMCP"  %in% class(candidate)) {
			graphs <- c(graphs, obj)
		}
	}
	if (length(graphs)==0) return("No graphMCP objects found.")
	return(graphs)
}

getObjectInfo <- function(object) {
  s <- paste(capture.output(print(object)), collapse="\n")
  descr <- attr(object, "description")
  if (!is.null(descr)) {
    s <- paste(s, "\nDescription:", descr, sep="\n")
  }
	return(s)
}

gMCPVersion <- function() {
	x <- try(as.character(packageVersion("gMCP")), silent=TRUE)
	if (class(x)!="try-error") {
		return(x)
	} else {
		return("unknown")
	}
}
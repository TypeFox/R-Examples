## Graph representation in gMCP
setClass("graphMCP",	
		representation(m="matrix", 
				weights="numeric", 
				nodeAttr="list", 
				edgeAttr="list"),
		validity=function(object) validWeightedGraph(object))

setMethod("initialize", "graphMCP",
		function(.Object, m, weights, nodeAttr=list(), edgeAttr=list()) {			
			if (length(weights)) {			
				checkValidWeights(weights)
			}			
      if (is.null(rownames(m))) {
        rownames(m) <- paste("H", 1:dim(m)[1], sep="")
      }
			colnames(m) <- rownames(m)
			.Object@m <- m
			names(weights) <- rownames(m)
			.Object@weights <- weights
			.Object@nodeAttr <- nodeAttr
			.Object@edgeAttr <- edgeAttr
			if(is.null(.Object@nodeAttr$rejected)) {
				.Object@nodeAttr$rejected <- rep(FALSE, dim(m)[1])
				names(.Object@nodeAttr$rejected) <- rownames(m)
			}
			validObject(.Object)
			return(.Object)
		})

validWeightedGraph <- function(object) {
	# if (sum(object@weights)>1)
	return(TRUE)
}

setClass("gMCPResult",		
		representation(graphs="list",
				pvalues="numeric",
				alpha="numeric",
				rejected="logical",
				adjPValues="numeric")
)

setMethod("print", "gMCPResult",
		function(x, ...) {
			callNextMethod(x, ...)
			#for (node in getNodes(x)) {
			#	cat(paste(node, " (",ifelse(unlist(nodeAttr(x, node, "rejected")),"rejected","not rejected"),", alpha=",format(unlist(nodeAttr(x, node, "nodeWeight")), digits=4 ,drop0trailing=TRUE),")\n", sep=""))	
			#}
			#cat(paste("alpha=",paste(format(getWeights(x), digits=4 ,drop0trailing=TRUE),collapse="+"),"=",sum(getWeights(x)),"\n", sep=""))			
		})

setMethod("show", "gMCPResult",
		function(object) {
			# callNextMethod(x, ...)
			cat("gMCP-Result\n")			
			cat("\nInitial graph:\n")
			print(object@graphs[[1]])
			cat("\nP-values:\n")
			print(object@pvalues)						
			if (length(object@adjPValues)>0) {
				cat("\nAdjusted p-values:\n")
				print(object@adjPValues)
			}
			cat(paste("\nAlpha:",object@alpha,"\n"))
			if (all(!object@rejected)) {
				cat("\nNo hypotheses could be rejected.\n")				
			} else {
				cat("\nHypothesis rejected:\n")
				print(object@rejected)
			}
			if (length(object@graphs)>1) {
				cat("\nFinal graph after", length(object@graphs)-1 ,"steps:\n")
				print(object@graphs[[length(object@graphs)]])
			}			
		})

setMethod("plot", "gMCPResult",
		function(x, y, ...) {
			# TODO Show visualization of graph			
		})

setGeneric("getNodes", function(object, ...) standardGeneric("getNodes"))

setMethod("getNodes", c("graphMCP"),
		function(object, ...) {			
			return(rownames(object@m))
		})

setGeneric("getMatrix", function(object, ...) standardGeneric("getMatrix"))

setMethod("getMatrix", c("graphMCP"),
		function(object, ...) {
			m <- object@m
			return(m)
		})

setGeneric("getWeights", function(object, node, ...) standardGeneric("getWeights"))

setMethod("getWeights", c("graphMCP"),
		function(object, node, ...) {
			weights <- object@weights
			names(weights) <- getNodes(object)
			if (!missing(node)) {
				return(weights[node])
			}
			return(weights)
		})

setGeneric("setWeights", function(object, weights, node, ...) standardGeneric("setWeights"))

setMethod("setWeights", c("graphMCP"),
		function(object, weights, node, ...) {
			if (missing(node)) {
				node <- getNodes(object)
			}
			object@weights[node] <- weights			
			return(object)
		})

setMethod("getWeights", c("gMCPResult"),
		function(object, node, ...) {
			graph <- object@graphs[[length(object@graphs)]]			
			return(getWeights(graph, node))
		})

setGeneric("setEdge", function(from, to, graph, weights) standardGeneric("setEdge"))

setMethod("setEdge", signature=signature(from="character", to="character",
				graph="graphMCP", weights="character"),
		function(from, to, graph, weights) {
			graph@m[from, to] <- weights
			graph
		})

setMethod("setEdge", signature=signature(from="character", to="character",
				graph="graphMCP", weights="numeric"),
		function(from, to, graph, weights) {
			graph@m[from, to] <- weights
			graph
		})

setGeneric("edgeAttr", function(self, from, to, attr) standardGeneric("edgeAttr"))
setGeneric("edgeAttr<-", function(self, from, to, attr, value) standardGeneric("edgeAttr<-"))

setMethod("edgeAttr", signature(self="graphMCP", from="character", to="character",
				attr="character"),
		function(self, from, to, attr) {
			self@edgeAttr[[attr]][from, to]
		})

setReplaceMethod("edgeAttr",
		signature(self="graphMCP", from="character", to="character", attr="character", value="ANY"),
		function(self, from, to, attr, value) {
			if (is.null(self@edgeAttr[[attr]])) self@edgeAttr[[attr]] <- matrix(NA, nrow=dim(self@m)[1], ncol=dim(self@m)[2])			
			rownames(self@edgeAttr[[attr]]) <- colnames(self@edgeAttr[[attr]]) <- getNodes(self)
			self@edgeAttr[[attr]][from, to] <- value		
			self
		})

setGeneric("nodeAttr", function(self, n, attr) standardGeneric("nodeAttr"))
setGeneric("nodeAttr<-", function(self, n, attr, value) standardGeneric("nodeAttr<-"))

setMethod("nodeAttr", signature(self="graphMCP", n="character", attr="character"),
		function(self, n, attr) {
			self@nodeAttr[[attr]][n]
		})

setReplaceMethod("nodeAttr",
		signature(self="graphMCP", n="character", attr="character", value="ANY"),
		function(self, n, attr, value) {
			if (is.null(self@nodeAttr[[attr]])) self@nodeAttr[[attr]] <- logical(length=length(getNodes(self)))
			self@nodeAttr[[attr]][n] <- value			
			self
		})

setGeneric("getRejected", function(object, node, ...) standardGeneric("getRejected"))

setMethod("getRejected", c("graphMCP"), function(object, node, ...) {
			rejected <- object@nodeAttr$rejected
			if (!missing(node)) {
				return(rejected[node])
			}
			return(rejected)
		})

setMethod("getRejected", c("gMCPResult"), function(object, node, ...) {			
			rejected <- object@rejected
			if (!missing(node)) {
				return(rejected[node])
			}
			return(rejected)
		})

setGeneric("setRejected", function(object, node, value) standardGeneric("setRejected"))
setGeneric("setRejected<-", function(object, node, value) standardGeneric("setRejected<-"))

setMethod("setRejected", c("graphMCP"),
		function(object, node, value) {
			if (missing(node)) {
				node <- getNodes(object)
			}
			object@nodeAttr$rejected[node] <- value			
			return(object)
		})

setReplaceMethod("setRejected", c("graphMCP"),
		function(object, node, value) {
			if (missing(node)) {
				node <- getNodes(object)
			}
			object@nodeAttr$rejected[node] <- value			
			return(object)
		})

setGeneric("getXCoordinates", function(graph, node) standardGeneric("getXCoordinates"))

setMethod("getXCoordinates", c("graphMCP"), function(graph, node) {
			x <- graph@nodeAttr$X
			if (is.null(x)) return(x)
			names(x) <- getNodes(graph)
			if (!missing(node)) {
				return(x[node])
			}
			return(x)
		})

setGeneric("getYCoordinates", function(graph, node) standardGeneric("getYCoordinates"))

setMethod("getYCoordinates", c("graphMCP"), function(graph, node) {
			y <- graph@nodeAttr$Y
			if (is.null(y)) return(y)
			names(y) <- getNodes(graph)
			if (!missing(node)) {
				return(y[node])
			}
			return(y)
		})

canBeRejected <- function(graph, node, alpha, pvalues) {	
	return(getWeights(graph)[[node]]*alpha>pvalues[[node]] | (all.equal(getWeights(graph)[[node]]*alpha,pvalues[[node]])[1]==TRUE));
}

setMethod("print", "graphMCP",
		function(x, ...) {
			callNextMethod(x, ...)
			#for (node in getNodes(x)) {
			#	cat(paste(node, " (",ifelse(unlist(nodeAttr(x, node, "rejected")),"rejected","not rejected"),", alpha=",format(unlist(nodeAttr(x, node, "nodeWeight")), digits=4 ,drop0trailing=TRUE),")\n", sep=""))	
			#}
			#cat(paste("alpha=",paste(format(getWeights(x), digits=4 ,drop0trailing=TRUE),collapse="+"),"=",sum(getWeights(x)),"\n", sep=""))			
		})

setMethod("show", "graphMCP",
		function(object) {
			#callNextMethod(object)
			nn <- getNodes(object)
			cat("A graphMCP graph\n")
			if (!isTRUE(all.equal(sum(getWeights(object)),1))) {
				cat(paste("Sum of weight: ",sum(getWeights(object)),"\n", sep=""))
			}			
			for (node in getNodes(object)) {
				cat(paste(node, " (",ifelse(nodeAttr(object, node, "rejected"),"rejected, ",""),"weight=",format(object@weights[node], digits=4, drop0trailing=TRUE),")\n", sep=""))
			}
			printEdge <- FALSE;
			for (i in getNodes(object)) {
				for (j in getNodes(object)) {
					if (object@m[i,j]!=0) {
						if (!printEdge) {
							cat("Edges:\n")
							printEdge <- TRUE
						}
						cat(paste(i, " -(", object@m[i,j], ")-> ", j, "\n"))
					}
				}
			}
			if (!printEdge) cat("No edges.\n")
			#if (!is.null(attr(object, "pvalues"))) { # TODO Do we want to have more output here?
			#  cat("\nAttached p-values: ", attr(object, "pvalues"), "\n")
			#}
			#cat(paste("\nalpha=",paste(format(getWeights(object), digits=4 ,drop0trailing=TRUE),collapse="+"),"=",sum(getWeights(object)),"\n", sep=""))
			cat("\n")
		}
)

getWeightStr <- function(graph, from, to, LaTeX=FALSE) {
	weight <- graph@m[from,to]	
	if (LaTeX) {
		if (is.numeric(weight)) {
			return(getLaTeXFraction(weight))
		} else {
      asNr <- try(eval(parse(text=weight), envir=NULL, enclos=NULL), silent=TRUE)
		  if (!is.na(asNr) && !("try-error"%in%class(asNr))) {
		    return(getLaTeXFraction(asNr))
		  }
		}
		# TODO / Bonus: Parse fractions in polynomials
		return(weight)
	}	
	if (is.numeric(weight)) {
		return(getFractionString(weight))
	}
	return(as.character(weight))	
}

getFractionString <- function(x, eps=1e-07) {
	xStr <- as.character(fractions(x))
	xStr <- ifelse(abs(sapply(xStr, function(x) {eval(parse(text=x))})-x)>eps, as.character(x), xStr)
	return(xStr)
}

setMethod("plot", "graphMCP",
		function(x, y, ...) {
			# TODO Show visualization of graph			
		})

setGeneric("simConfint", function(object, pvalues, confint, alternative=c("less", "greater"), estimates, df, alpha=0.05, mu=0) standardGeneric("simConfint"))

setMethod("simConfint", c("graphMCP"), function(object, pvalues, confint, alternative=c("less", "greater"), estimates, df, alpha=0.05, mu=0) {
			result <- gMCP(object, pvalues, alpha=alpha)
			if (all(getRejected(result))) {
				alpha <- getWeights(object)*alpha				
			} else {
				alpha <- getWeights(result)*alpha				
			}
			if (class(confint)=="function") {
				f <- function(node, alpha, rejected) {
					if (rejected && alternative=="less") return(c(-Inf, mu))
					if (rejected && alternative=="greater") return(c(mu, Inf))
					return(confint(node, alpha))
				}
				m <- mapply(f, getNodes(object), alpha, getRejected(result))	
				m <- rbind(m[1,], estimates, m[2,])
				rownames(m) <- c("lower bound", "estimate", "upper bound")
				return(t(m))
			} else if (confint=="t") {
				dist <- function(x) {qt(p=x, df=df)}
			} else if (confint=="normal") {
				dist <- qnorm
			} else {
				stop('Parameter confint has to be a function or "t" or "normal".')
			}
			if (alternative=="greater") {			
				stderr <- abs(estimates/dist(1-pvalues))
				lb <- estimates+dist(alpha)*stderr				
				lb <- ifelse(getRejected(result), max(0,lb), lb) 
				ub <- rep(Inf,length(lb))
			} else if (alternative=="less") {			
				stderr <- abs(estimates/dist(pvalues))								 
				ub <- estimates+dist(1-alpha)*stderr				
				ub <- ifelse(getRejected(result), min(0,ub), ub)
				lb <- rep(-Inf,length(ub))
			} else {
				stop("Specify alternative as \"less\" or \"greater\".")
			}
			m <- matrix(c(lb, estimates, ub), ncol=3)
			colnames(m) <- c("lower bound", "estimate", "upper bound")
			return(m)
		})


setClass("gPADInterim",
		representation(Aj="matrix",
				BJ="numeric",
				z1="numeric",
				v="numeric",
				preplanned="graphMCP",
				alpha="numeric"),
		validity=function(object) validPartialCEs(object))

setMethod("print", "gPADInterim",
		function(x, ...) {
			callNextMethod(x, ...)
			
		})


setMethod("show","gPADInterim",
		function(object) {
			cat("Pre-planned graphical MCP at level:",object@alpha,"\n")
			show(object@preplanned)
			n <- length(object@z1)
			cat("Proportion of pre-planned measurements\n collected up to interim:\n")
			v <- object@v
                        if(length(v) == 1) v <- rep(v,n)
			names(v) <- paste('H',1:n,sep='')
                        print(v)
			cat("Z-scores computed at interim\n")
			z1 <- object@z1
			names(z1) <- paste('H',1:n,sep='')
			print(z1)
			cat("\n Interim PCE's by intersection\n")
			tab <- round(cbind(object@Aj,object@BJ),3)
			rownames(tab) <- to.intersection(1:nrow(tab))
			colnames(tab) <- c(paste('A(',1:n,')',sep=''),'BJ')
			print(tab)
		}
)

############################## Entangled graphs #################################

## Entangled graph representation in gMCP
setClass("entangledMCP",	
		representation(subgraphs="list", 
				weights="numeric",
				graphAttr="list"),
		validity=function(object) validEntangledGraph(object))


validEntangledGraph <- function(graph) {
	if (!all("graphMCP" == lapply(graph@subgraphs, class))) stop("Subgraphs need to be of class 'graphMCP'.")
	return(TRUE)
}

setGeneric("getMatrices", function(object, ...) standardGeneric("getMatrices"))

setMethod("getMatrices", c("entangledMCP"),
		function(object, ...) {
			result <- list()
			for (g in object@subgraphs) {
				result[[length(result)+1]] <- g@m
			}
			return(result)
		})

setMethod("getWeights", c("entangledMCP"),
		function(object, node, ...) {
			result <- c()
			for (g in object@subgraphs) {
				result <- rbind(result, getWeights(g, node, ...))
			}
			return(result)
		})

setMethod("getNodes", c("entangledMCP"),
		function(object, ...) {			
			return(getNodes(object@subgraphs[[1]]))
		})

setMethod("getXCoordinates", c("entangledMCP"), function(graph, node) {			
			return(getXCoordinates(graph@subgraphs[[1]], node))
		})

setMethod("getYCoordinates", c("entangledMCP"), function(graph, node) {
			return(getYCoordinates(graph@subgraphs[[1]], node))
		})

setMethod("getRejected", c("entangledMCP"), function(object, node, ...) {
			return(getRejected(object@subgraphs[[1]], node))
		})


gMCPenv <- new.env()

checkValidWeights <- function(weights) {
	if(!is.numeric(weights)) {
		stop("Weights have to be numeric!")
	}
	if(any(is.na(weights) | is.infinite(weights))) {
		warning("Some of the weights are not real numbers. NA, NaN, Inf and -Inf are not supported.")
	}
	if(any(0 > weights | weights > 1 + .Machine$double.eps ^ 0.25)) {
		warning("Invalid weights: weights must be between 0 and 1")
	}
	if(sum(weights) > 1 + .Machine$double.eps ^ 0.25) {
		warning("Invalid weights: the sum of all weights must be less than 1")
	}
}

# Converts a string like "5+3*e+5*e^2" to the tupel representation c(5,3,5) 
parseEpsPolynom <- function(s) {
	env <- new.env(parent = baseenv())
	assign("epsilon", polynom(), envir=env)
	p <- try(eval(parse(text=s), envir=env), silent = TRUE)
	if (class(p)=="try-error") {
		stop("String does not represent a polynom in epsilon.")
	}
	if(is.numeric(p)) {
		return(p)
	} else {
    	return(coef(p))
	}
}

getDebugInfo <- function() {
	graphs <- ls(pattern="\\.InitialGraph*", all.names=TRUE, envir=gMCPenv)
	graphs <- c(graphs, ls(pattern="\\.tmpGraph*", all.names=TRUE, envir=gMCPenv))
	graphs <- c(graphs, ls(pattern="\\.ResetGraph*", all.names=TRUE, envir=gMCPenv))
	graphInfo <- c()
	for (graph in graphs) {
		.DebugGraph <- get(graph, envir=gMCPenv)
		graphTXT <- paste(capture.output(print(.DebugGraph)), collapse="\n")
		matrixTXT <- paste("m <-",paste(capture.output(dput(graph2matrix(.DebugGraph))), collapse="\n"),"\n")
		weightsTXT <- paste("w <-",paste(capture.output(dput(getWeights(.DebugGraph))), collapse="\n"),"\n")
		createTXT <- paste("graph <- matrix2graph(m)", "setWeights(graph, w)", sep="\n")
		graphInfo <- c(graphInfo, paste(graphTXT, matrixTXT, weightsTXT, createTXT, sep="\n"))
	}
	if (length(graphInfo)!=0) {
		return(paste(graphInfo, collapse="\n\n"))
	}
	return("Graph not available.")
}

#' Create a Block Diagonal Matrix with NA outside the diagonal
#' 
#' Build a block diagonal matrix with NA values outside the diagonal given
#' several building block matrices.
#' 
#' This function is usefull to build the correlation matrices, when only
#' partial knowledge of the correlation exists.
#' 
#' @param ...  individual matrices or a \code{list} of matrices.
#' @return A block diagonal matrix with NA values outside the diagonal.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{gMCP}}
#' @examples
#' 
#' 
#' bdiagNA(diag(3), matrix(1/2,nr=3,nc=3), diag(2))
#' 
#' 
#' @export bdiagNA
bdiagNA <- function(...) {	
	if (nargs() == 0) 
		return(matrix(nrow=0, ncol=0))
	if (nargs() == 1 && !is.list(...)) 
		return(as.matrix(...))
	asList <- if (nargs() == 1 && is.list(...)) list(...)[[1]] else list(...)
	if (length(asList) == 1) 
		return(as.matrix(asList[[1]]))
	n <- 0
	for (m in asList) {
		if (!is.matrix(m)) {
			stop("Only matrices are allowed as arguments.")
		}
		if (dim(m)[1]!=dim(m)[2]) {
			stop("Only quadratic matrices are allowed.")
		}
		n <- n + dim(m)[1]	
	}
	M <- matrix(NA, nrow=n, ncol=n)
	k <- 0
	for (m in asList) {
		for (i in 1:dim(m)[1]) {
			for (j in 1:dim(m)[1]) {
				M[i+k,j+k] <- m[i,j]
			}
		}
		k <- k + dim(m)[1]	
	}	
	return(M)
}

requireLibrary <- function(package) {
	if(!require(package, character.only=TRUE)) {
		answer <- readline(paste("Package ",package," is required - should we install it (y/n)? ", sep=""))
		if (substr(answer, 1, 1) %in% c("y","Y")) {
			if (package %in% c("graph", "Rgraphviz")) {	
				source("http://www.bioconductor.org/biocLite.R")
				biocLite <- get("biocLite", envir=globalenv())
				biocLite(package)
			} else {
				install.packages(package)				
			}
			return(require(package, character.only=TRUE))
		} else {
			stop(paste("Required package ",package," should not be installed.", sep=""))
		}
	} else {
		return(TRUE)
	}
}

triangle <- function(min, peak, max) {
	
}

# Depending of the number of hypotheses it is calculated 
# which standard designs could be a possibility and how many 
# groups they would have.
getAvailableStandardDesigns <- function(n) {
	designs <- c()
	numberOfGroups <- c()
	possibleDesigns = c("Dunnett", "Tukey", "Sequen", "AVE", "Changepoint", "Williams", "Marcus", "McDermott", "UmbrellaWilliams", "GrandMean")
	for (design in possibleDesigns) {
		for (i in 2:(2*n)) {
			m <- try(contrMat(n=rep(10, i), type=design), silent = TRUE)
			if (!("try-error" %in% class(m)) && dim(m)[1]==n) {
				designs <- c(designs, design)
				numberOfGroups <- c(numberOfGroups, i)
			}
		}
	}	
	return(list(designs, numberOfGroups))
}

getCorrMat <- function(n, type) {
	m <- contrMat(n, type)	
	m <- m%*%diag(1/sqrt(n))%*%t(m)
	s <- diag(1/sqrt(diag(m)))
	r <- s%*%m%*%s
	return(as.numeric(r))
}

# Given a character string x only the letters are returned (upper and lower case)
removeSymbols <- function(x, numbers=TRUE) {
	pattern <-ifelse(numbers, "[a-z]|[0-9]", "[a-z]")
	m <- gregexpr(pattern, x, ignore.case=TRUE)
	return(unlist(lapply(regmatches(x, m), paste, collapse="")))
}

nextAvailableName <- function(x) {
	if (exists(x)) {
		i <- 1
		while (exists(paste(x, i, sep=""))) {
			i <- i + 1
		}
		x <- paste(x, i, sep="")
	}
	return(x)
}

layers <- function(graph) {
	if ("graphMCP" %in% class(graph)) {
		return(as.integer(1))
	} else if ("entangledMCP" %in% class(graph)){
		return(as.integer(length(graph@subgraphs)))
	} else {
		stop("This function should only be used for objects of class graphMCP or entangledMCP.")
	}
}

#' Get Memory and Runtime Info from JVM
#' 
#' Get Memory and Runtime Info from JVM
#' 
#' @param memory Logical whether to include memory information + number of available cores
#' @param filesystem Logical whether to include filesystem information (Total, free and usable space)
#' @param runtime Logical whether to include runtime information (Class Path, Library Path, Input Arguments)
#' @return character vector of length 1 containing the memory and runtime info.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @examples
#' 
#' cat(getJavaInfo())
#' 
getJavaInfo <- function(memory=TRUE, filesystem=TRUE, runtime=TRUE) {
  return(.jcall("org.af.commons.logging.SystemInfo", "S", "getSystemInfo", memory, filesystem, runtime))  
}

getNewestRVersion <- function() {
  # What about something like setInternet2(true)?
  warn <- getOption("warn")
  options(warn=-1)
  line <- try(grep("R-[0-9.]+-win", readLines("http://cran.r-project.org/bin/windows/base/", warn=FALSE), value=TRUE), silent=TRUE)
  if("try-error"==class(line) || length(line)==0) {
    return("Unknown")
  }
  options(warn=warn)
  return(strsplit(regmatches(line[1], regexpr("R-[0-9.]+-win", line[1])), split="-")[[1]][2])
}
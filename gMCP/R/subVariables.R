#' Substitute Epsilon
#' 
#' Substitute Epsilon with a given value.
#' 
#' For details see the given references.
#' 
#' @param graph A graph of class \code{\link{graphMCP}} or class
#' \code{\link{entangledMCP}}.
#' @param eps A numeric scalar specifying a value for epsilon edges.
#' @return A graph where all epsilons have been replaced with the given value.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{graphMCP}}, \code{\link{entangledMCP}}
#' @keywords print graphs
#' @examples
#' 
#' 
#' graph <- improvedParallelGatekeeping()
#' graph
#' substituteEps(graph, eps=0.01)
#' 
#' 
#' @export substituteEps
substituteEps <- function(graph, eps=10^(-3)) {
	# Call this function recursivly for entangled graphs.
	if ("entangledMCP" %in% class(graph)) {
		for(i in 1:length(graph@subgraphs)) {
			graph@subgraphs[[i]] <- substituteEps(graph@subgraphs[[i]], eps)
		}
		return(graph)
	}
	# Real function:
	if (is.numeric(graph@m)) return(graph)
	m <- matrix(gsub("\\\\epsilon", eps, graph@m), nrow=length(getNodes(graph)))
	options(warn=-1)
	m2 <- matrix(sapply(m, function(x) {
						result <- try(eval(parse(text=x)), silent=TRUE);
						ifelse(class(result)=="try-error",NA,result)
					}), nrow=length(getNodes(graph)))
	options(warn=0)
	if (all(is.na(m)==is.na(m2))) m <- m2
	rownames(m) <- colnames(m) <- getNodes(graph)
	graph@m <- m
	return(graph)
}

#' Replaces variables in a general graph with specified numeric values
#' 
#' Given a list of variables and real values a general graph is processed and
#' each variable replaced with the specified numeric value.
#' 
#' 
#' @param graph A graph of class \code{\link{graphMCP}} or class
#' \code{\link{entangledMCP}}.
#' @param variables A named list with one or more specified real values, for example
#' \code{list(a=0.5, b=0.8, "tau"=0.5)} or \code{list(a=c(0.5, 0.8), b=0.8, "tau"=0.5)}.  
#' If \code{ask=TRUE} and this list is
#' missing at all or single variables are missing from the list, the user is
#' asked for the values (if the session is not interactive an error is thrown).
#' For interactively entered values only single numbers are supported.
#' @param ask If \code{FALSE} all variables that are not specified are not
#' replaced.
#' @param partial IF \code{TRUE} only specified variables are replaced and 
#' parameter \code{ask} is ignored.
#' @param expand Used internally. Don't use yourself.
#' @param list If \code{TRUE} the result will always be a list, even if only one
#' graph is returned in this list.
#' @return A graph or a matrix with variables replaced by the specified numeric
#' values. Or a list of theses graphs and matrices if a variable had more than one value.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{graphMCP}}, \code{\link{entangledMCP}}
#' @keywords print graphs
#' @examples
#' 
#' 
#' graph <- HungEtWang2010()
#' \dontrun{
#' replaceVariables(graph)
#' }
#' replaceVariables(graph, list("tau"=0.5,"omega"=0.5, "nu"=0.5))
#' replaceVariables(graph, list("tau"=c(0.1, 0.5, 0.9),"omega"=c(0.2, 0.8), "nu"=0.4))
#' 
#' @export replaceVariables
replaceVariables <-function(graph, variables=list(), ask=TRUE, partial=FALSE, expand=TRUE, list=FALSE) {
  if (expand) variables <- varcombs(variables)
  if (is.list(variables[[1]])) {
    result <- list()
    for (v in variables) {
      r2 <- replaceVariables(graph, v, ask, partial, expand=FALSE)
      attr(r2, "label") <- attr(v, "label")
      attr(r2, "variables") <- unlist(v)
      result <- c(result, list(r2))
    }
    if (length(result)==1 && !list) return(result[[1]])
    return(result)
  }
	# Call this function recursivly for entangled graphs.
	if ("entangledMCP" %in% class(graph)) {
		for(i in 1:length(graph@subgraphs)) {
			graph@subgraphs[[i]] <- replaceVariables(graph@subgraphs[[i]], variables, ask)
		}
		return(graph)
	}
	# Real function:
	greek <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", 
			"theta", "iota", "kappa", "lambda", "mu", "nu", "xi", 
			"omicron", "pi", "rho", "sigma", "tau", "nu", "phi",
			"chi", "psi", "omega")
	if (is.matrix(graph)) { m <- graph } else {m <- graph@m}	
	for (g in c(greek,  letters)) {
		if (length(grep(g, m))!=0) {
			if (is.null(answer <- variables[[g]])) {
				if (!partial && ask) {
					if(interactive()) {
						answer <- readline(paste("Value for variable ",g,"? ", sep=""))
					} else {
						stop(paste("Value for variable",g,"not specified."))
					}
				}
			}
			if(!is.null(answer)) {
				m <- gsub(paste(ifelse(nchar(g)==1,"","\\\\"), g, sep=""), answer, m)
			}
		}
	}
	if (is.matrix(graph)) return(parse2numeric(m))
	graph@m <- m
	return(parse2numeric(graph))
}

# Parses matrices of graphs (simple and entangled)
# 
# Parses matrices of graphs (simple and entangled) when values are of type character, e.g. "1/3".
# 
# @param g Graph of class \code{\link{graphMCP}} or \code{\link{entangledMCP}}.
# @param force Logical whether conversion to numeric should be forced or not. 
# If forced all values that could not be parsed will be \code{NA}. 
# Otherwise the original unchanged graph will be returned.
# @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
# @keywords Converted graph (if all values could be parsed or \code{force=TRUE}) or original graph.
# @examples
#
# # Nothing changes:
# gMCP:::parse2numeric(HungEtWang2010())
# # Note that other methods like printing don't handle NAs well:
# gMCP:::parse2numeric(HungEtWang2010(), force=TRUE)
# 
parse2numeric <- function(graph, force=FALSE) {
	# Call this function recursivly for entangled graphs.
	if ("entangledMCP" %in% class(graph)) {
		for(i in 1:length(graph@subgraphs)) {
			graph@subgraphs[[i]] <- parse2numeric(graph@subgraphs[[i]])
		}
		return(graph)
	}
	# Real function:
	if (is.matrix(graph)) { m <- graph } else {m <- graph@m}
	names <- rownames(m)
	m <- matrix(sapply(m, function(x) {
						result <- try(eval(parse(text=x)), silent=TRUE);
						ifelse(class(result)=="try-error",NA,result)
					}), nrow=dim(m)[1])
	if (!force && any(is.na(m))) return(graph)
	rownames(m) <- colnames(m) <- names
	if (is.matrix(graph)) return(m)
	graph@m <- m
	return(graph)
}

isEpsilon <- function(w) {
	x <- try(eval(parse(text = gsub("\\\\epsilon", 0, w)), envir = emptyenv()), silent=TRUE)
	if ("try-error" %in% class(x)) return(FALSE)
	return(x==0)
}

# For testing purposes: variables <- list(a=c(1,2), b=(3), x=c(2,3,4), d=c(1,2))
varcombs <- function(variables) {
  combs <- list()
  m <- do.call(expand.grid, lapply(variables, function(x){1:length(x)}))
  for (i in 1:dim(m)[1]) {
    variablesII <- rep(0, length(variables))
    for(k in 1:length(variables)) {
      variablesII[k] <- variables[[k]][m[i,k]]
    }
    names(variablesII) <- names(variables)
    x <- as.list(variablesII)
    attr(x, "label") <- paste(paste(names(variables),"=",variablesII,sep=""), collapse=", ")
    combs <- c(combs, list(x))
  }
  # GII <- replaceVariables(G, as.list(variablesII))
  # additionalLabel <- paste(",", paste(paste(names(variables),"=",variablesII,sep=""), collapse=", "))
  return(combs)
}
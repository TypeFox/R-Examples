#' Comparing and testing network objects
#' 
#' Compare or test network objects for (near) equality.
#' 
#' Arguments \code{target} and \code{current} can be network objects of one of
#' the supported classes. They do not have to be of the same class though.
#' 
#' The function does a series of comparisons between \code{target} and
#' \code{current}:
#'
#' \enumerate{
#' \item The network structure comparison is made based on adjacency matrices
#' (mind this when using for huge networks).
#' 
#' \item Network/edge/vertex attributes are checked for presence in both
#' objects.
#' 
#' \item Common network/edge/vertex attribures are checked for equality.
#' }

#' All the results are collected in a list of class \code{netcompare} with an
#' associated \code{print} method.
#' 
#' If \code{test} is TRUE then instead of the detailed test results the function
#' returns TRUE or FALSE depending on some of the checks resulted positively.
#' Currently attribute checks are ignored, i.e. what is taken into account is:
#'
#' \enumerate{
#' \item Equivalence of adjacency matrices
#' \item Directed / undirected character of the network
#' \item Edge set size
#' \item Vertex set size
#' }
#' 
#' @param target,current network objects, currently \code{network} and
#' \code{igraph} classes are supported
#'
#' @param test logical, whether to perform the test or return comparison data,
#' see Details
#'
#' @param \dots other arguments, currently ignored
#'
#' @return Depending on the value of \code{test} either an object of class
#' \code{netcompare} containing the results of all the tests (if
#' \code{test=FALSE}) or (if \code{test=TRUE}) a logical whether or not the
#' networks are (nearly) the same.
#'
#' @seealso \code{\link{all.equal}}, \code{\link{identical}}
#'
#' @export
#'
#' @examples
#'
#'netcompare( asIgraph(exNetwork), exNetwork)
#'netcompare( asIgraph(exNetwork), exNetwork, test=TRUE)
#'
netcompare <- function(target, current, test=FALSE, ...)
{
  # trivial checks
  rval <- list()
  # class
  rval$class <- c(target=class(target), current=class(current))
  # number of vertices
  rval$vcount <- c(target=igVcount(target), current=igVcount(current))
  # number of edges
  rval$ecount <- c( target=igEcount(target), current=igEcount(current))
  # directedness
  rval$directed <- c( target=igDirected(target), current=igDirected(current))
  # compare adjacency matrices
  rval$identical_am <- compareEdges(target, current)
  # compare attributes
  targeta <- dumpAttr(target)
  currenta <- dumpAttr(current)
  rval$network <- compareAttributes( targeta$network, currenta$network)
  rval$vertex <- compareAttributes( targeta$vertex, currenta$vertex)
  rval$edge <- compareAttributes( targeta$edge, currenta$edge)
  rval <- structure(rval, class="netcompare")
  if(test)
    compareTest(rval)
  else return(rval)
}







print.netcompare <- function(x, ...)
{
  cat("\n")
  cat("Identical adjacency matrices:\n")
  cat( paste(x$identical_am, collapse=", "), "\n", fill=TRUE, labels="   ")
  cat("Network-level features:\n")
  m <- do.call("rbind", lapply(x[c("vcount", "ecount", "directed")], format))
  print(m, quote=FALSE)
  cat("\n")
  cat("Presence of network-level attributes\n")
  print(x$network)
  cat("\n")
  cat("Presence of vertex-level attributes\n")
  print(x$vertex)
  cat("\n")
  cat("Presence of edge-level attributes\n")
  print(x$edge)
}

# comparing and testing for (near) equality of networks
#
# Result of comparison:
#
# 1. computed network-level comparisons
#
# 2. built-in network-level comparisons
#
# 3. attributes
#
# 3a. attribute presence
#
# 3b. identical common attributes
#
# NOTE: Makes use of non-exported generic functions


compareTest <- function(object)
{
  stopifnot(inherits(object, "netcompare"))
  rval <- logical()
  rval["adjacency"] <- object$identical_am
  rval["vcount"] <- object$vcount[1] == object$vcount[2]
  rval["ecount"] <- object$ecount[1] == object$ecount[2]
  rval["directed"] <- object$directed[1] == object$directed[2]
  all(rval)
}


compareEdges <- function(target, current, use.names=FALSE)
{
  tr <- try(utils::getS3method("as.matrix", class=class(target)), silent=TRUE)
  if(inherits(tr, "try-error"))
    stop("cannot find 'as.matrix' method for class ", dQuote(class(target)))
  tr <- try(utils::getS3method("as.matrix", class=class(current)), silent=TRUE)
  if(inherits(tr, "try-error"))
    stop("cannot find 'as.matrix' method for class ", dQuote(class(current)))
  mtar <- as.matrix(target, "adjacency")
  mcur <- as.matrix(current, "adjacency")
  # compare matrices (no dimnames)
  if(use.names)
    all.equal(mtar, mcur)
  else all.equal( structure(mtar, dimnames=NULL), structure(mcur,
            dimnames=NULL) )
}


# compare common components of a list (by name)
# return a list of all.equal results
compareAlist <- function(target, current)
{
  # common components
  nams <- intersect(names(target), names(current))
  if( length(nams) == 0 )
    return(as.character(NA))
  rval <- lapply(nams, function(n) all.equal( target[[n]], current[[n]]) )
  names(rval) <- nams
  rval
}



# Compare lists of attributes (as returned by 'dumpAttr')

compareAttributes <- function(target, current)
{
  # compare number of attributes
  rval <- list()
  rval$n <- c(target=length(target), current=length(current))
  # compare names
  pre <- list()
  u <- union(names(target), names(current))
  r <- t(sapply(u, function(a)
      c( a %in% names(target),
          a %in% names(current) )
      ))
  pre <- c(pre, list(r))
  pre <- do.call("rbind", pre)
  dimnames(pre) <- list(rownames(pre), c("target", "current"))
  rval$presence <- pre
  rval$bycomp <- compareAlist(target, current)
  structure(rval, class="netcomparea")
}


# Print method for the result of 'compareAttributes'
print.netcomparea <- function(x, ...)
{
  m <- do.call("rbind", lapply( x[c("n", "presence")], format))
  print(m, quote=FALSE)
  cat("Common attributes comparison (TRUE=identical)\n")
  if( identical( x$bycomp, as.character(NA)) )
  {
    cat("   No common attributes\n")
  } else
  {
    l <- sapply(x$bycomp, paste, collapse=", ")
    for(i in seq(along=l))
      cat(names(l)[i], ":", l[i], fill=TRUE, labels="   ")
  }
}

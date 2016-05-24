#' Analyzes Clickstreams Based on Markov Chains
#' 
#' This package allows modeling clickstreams with Markov chains. It supports to
#' model clickstreams as zero-order, first-order or higher-order Markov chains.
#' 
#' \tabular{ll}{ Package: \tab clickstream\cr Type: \tab Package\cr Version:
#' \tab 1.1.8\cr Date: \tab 2016-04-26\cr License: \tab GPL-2\cr Depends: \tab
#' R (>= 3.0), methods\cr }
#' 
#' @name clickstream-package
#' @aliases clickstream-package clickstream
#' @docType package
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @references For modeling clickstreams with Markov chains, see Ching, W.-K.
#' et al.: \emph{Markov Chains -- Models, Algorithms and Applications}, 2nd
#' edition, Springer, 2013.
#' @import methods Rsolnp arules data.table plyr linprog
#' @importFrom stats runif rpois kmeans cor
#' @importFrom utils read.table write.table count.fields
#' @importFrom igraph E graph.adjacency
#' @keywords click stream Markov chain
#' @examples
#' 
#' # fitting a simple Markov chain and predicting the next click
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf <- tempfile()
#' writeLines(clickstreams, csf)
#' cls <- readClickstreams(csf, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' startPattern <- new("Pattern", sequence = c("h", "c"))
#' predict(mc, startPattern)
#' plot(mc)
#' 
NULL

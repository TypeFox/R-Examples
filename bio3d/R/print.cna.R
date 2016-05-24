print.cna <- function(x, ...) {

  ## Check for presence of igraph package
  oops <- requireNamespace("igraph", quietly = TRUE)
  if (!oops) {
     stop("igraph package missing: Please install, see: ?install.packages")
  }
  
  ## y <- summary.cna(x, ...)

  l1 <- paste( "\n - NETWORK NODES#:  ", x$communities$vcount,
              "\tEDGES#:", igraph::ecount(x$network))
 
  l2 <- paste( "\n - COMMUNITY NODES#:", max(x$communities$membership),
              "\tEDGES#:", igraph::ecount(x$community.network))

  cat("\nCall:\n  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nStructure:",l1,l2,"\n\n ")
  
  i <- paste( attributes(x)$names, collapse=", ")
  cat(strwrap(paste(" + attr:",i,"\n"),width=60, exdent=8), sep="\n")
  #print.igraph(x$network)
  #print.igraph(x$community.network)
}

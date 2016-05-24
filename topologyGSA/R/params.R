.procParams <- function(y1, y2, dag) {
  if (!is.matrix(y1) && !is.data.frame(y1))
    stop("y1 is not a matrix or a data frame.")
  else if (!is.matrix(y2) && !is.data.frame(y2))
    stop("y2 is not a matrix or a data frame.")
  else if (ncol(y1) != ncol(y2))
    stop("y1 and y2 differ in the number of columns (genes)")
  else if (any(colnames(y1) != colnames(y2)))
    stop("y1 and y2 differ in the column names (gene names)")
  else if (nrow(y1) < 3)
    stop("y1 should have at least 3 rows (samples)")
  else if (nrow(y2) < 3)
    stop("y2 should have at least 3 rows (samples)")

  common <- intersect(colnames(y1), nodes(dag))
  if (length(common) < 3)
    stop("need at least 3 genes in common between expression and dag")

  y1  <- y1[,common,drop=FALSE]
  y2  <- y2[,common,drop=FALSE]
  dag   <- subGraph(common, dag)
  graph <- .processGraph(dag)

  list(y1=y1, y2=y2, graph=graph)
}

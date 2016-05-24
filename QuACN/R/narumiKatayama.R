narumiKatayama <- function(g, deg=NULL) {
  if (is.null(deg)) {

    if (class(g)[1] != "graphNEL")
      stop("'g' has to be a 'graphNEL' object")
    stopifnot(.validateGraph(g))
    
    deg = graph::degree(g)
  }

  prod(deg)
}

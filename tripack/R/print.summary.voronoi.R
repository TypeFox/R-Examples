print.summary.voronoi<-function(x,...)
  {
    cat("voronoi mosaic\n")
    cat("Call:", deparse(x$call),"\n")
    cat(x$nn, "nodes\n")
    cat(x$nd, "dummy nodes\n")
  }

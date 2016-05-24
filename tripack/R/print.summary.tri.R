print.summary.tri<-function(x,...)
  {
    cat("triangulation:\n")
    cat("Call:", deparse(x$call),"\n")
    cat("number of nodes:",x$n,"\n")
    cat("number of arcs:",x$na,"\n")
    cat("number of boundary nodes:",x$nb,"\n")
    cat("number of triangles:",x$nt,"\n")
    cat("number of constraints:",x$nc,"\n")
  }

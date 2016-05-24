##
## plot (gRain)
##

plot.grain <- function(x, type, ...){
  #cat("plot.grain; type:", type, "\n")
   #' if (!require("Rgraphviz")){
   #'   cat("The Rgraphviz package (from Bioconductor) must be installed to display the models\n")
   #'   return()
   #' }
    
   if (!requireNamespace("Rgraphviz", quietly = TRUE)) {
       cat("The Rgraphviz package (from Bioconductor) must be installed to display the models\n")
       return()
   }

   if (missing(type)){
       if (x$isCompiled){
           Rgraphviz::plot(x$ug)
       } else {
           if ("pot-grain" %in% class(x)){
               Rgraphviz::plot(x$ug)
           } else {
               Rgraphviz::plot(x$dag)
           }
       }
   } else {
       if (type=="jt") ## For backward compatibility; p 57 in GMwR-book
           type="rip"
       zz <- x[[type]]
       if (!is.null(zz))
           Rgraphviz::plot(zz)
       else
           cat("Slot", type, "does not exist \n")
   }
}

iplot.grain <- function(x,type, ...){
  #.primiplot(x$dag)

  if (missing(type)){
    if (x$isCompiled){
      .primiplot(x$ug)
    } else {
      if ("pot-grain" %in% class(x)){
        .primiplot(x$ug)
      } else {
        .primiplot(x$dag)
      }
    }
  } else {
    zz <- x[[type]]
    if (!is.null(zz))
      .primiplot(zz)
    else
      cat("Slot", type, "does not exist \n")
  }
}

.primiplot <- function(grp){
  ig <- igraph.from.graphNEL(grp)
  V(ig)$label <- V(ig)$name
  V(ig)$size  <- 40
  ig$cex   <-  4
  ig$layout   <- layout.graphopt
  plot(ig)
}




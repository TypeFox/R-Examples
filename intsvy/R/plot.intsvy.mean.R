plot.intsvy.mean <- function(x, se=TRUE, sort=FALSE, ...) {
  # it is assumed that last three columns are:  "Freq"       "Mean" "Std.err."  
  vars <- setdiff(colnames(x), c("Freq", "Mean", "s.e.", "SD", "s.e"))
  nvar <- length(vars)
  x$MeanL <- x$Mean - x$s.e.
  x$MeanH <- x$Mean + x$s.e.
  pl <- NA
  if (isTRUE(sort)) {
    x[,1] <- reorder(factor(x[,1]), x[,"Mean"], mean, na.rm=TRUE)
  }
  if (nvar == 1) {
    pl <- ggplot(data=x, aes_string(y = "Mean", x=vars[1])) + 
      geom_point(size=5) + 
      theme_bw() + coord_flip() + 
      theme(legend.position="top")
  } else {
    if (nvar == 2) {
      pl <- ggplot(data=x, aes_string(y = "Mean", x=vars[1], col=vars[2])) + 
        geom_point(size=5) + 
        theme_bw() + coord_flip() + 
        theme(legend.position="top")
    } else {
      pl <- ggplot(data=x, aes_string(y = "Mean", x=vars[1], col=vars[2], shape= vars[3])) + 
        geom_point(size=5) + 
        theme_bw() + coord_flip() + 
        facet_wrap(as.formula(paste0("~", vars[3])))+ 
        theme(legend.position="top")
    }
  }
  
  if (se) {
    pl <- pl + geom_errorbar(aes_string(ymin="MeanL", ymax="MeanH"), width=.5) 
  } 
  pl
}


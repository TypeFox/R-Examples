"plot.kver" <- function(x, which = names(x), colpol=rainbow(length(which)),
                        colborder=rep("black", length(which)), lwd = 2,
                        add=FALSE, ...)
{
    ## Verifications
    if (!inherits(x, "kver"))
        stop("x should be of class kver")
    x <- x[which]

    ## A first plot if the plot is not to be added to another plot
    if (!add) {
      xc <- unlist(lapply(x, function(y) y[,2]))
      yc <- unlist(lapply(x, function(y) y[,3]))
      plot(xc, yc, asp=1, ty="n", ...)
    }

    ## The plot
    lapply(1:length(x),
           function(i) plot.area(x[[i]],
                                 colpol = rep(colpol[i], nlevels(x[[i]][,1])),
                                 colborder = rep(colborder[i],
                                   nlevels(x[[i]][,1])),
                                 lwd = lwd, add = TRUE))
    invisible(NULL)
  }


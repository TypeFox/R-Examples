
images <- function(matrices, names=NULL, ...) {

  no.rows <- nrow(matrices[[1]])
  no.cols <- ncol(matrices[[1]])
  no.mods <- length(matrices)

  if (no.mods==0) {
    stop("No matrices to plot")
  }

  alldata <- structure(unlist(matrices), dim=c(no.rows, no.cols, no.mods),
                       dimnames=list(NULL, NULL, names))

  rx <- range(alldata, finite=TRUE)
  nn <- 100
  n0 <- min(nn, max(0, round((0-rx[1])/(rx[2]-rx[1])*nn)))
  col.regions <- c(colorRampPalette(c("blue3", "gray80"))(n0),
                   colorRampPalette(c("gray75", "red3"))(nn-n0))

  levelplot(alldata, scale=list(y=list(draw=FALSE), x=list(draw=FALSE)),
            col.regions=col.regions, between=list(x=0.5,y=0.5),
            ...)
} 

setMethod("plotModules", signature(modules="list"),
          function(modules, ...)
          plotModules.default(modules, ...))

plotModules.default <- function(modules, to.plot=seq_len(ncol(modules$rows)),
                                data, binary=TRUE, names=NULL,
                                xlab="", ylab="",
                                ...) {

  no.mods <- length(to.plot)
  no.rows <- nrow(modules$rows)
  no.cols <- nrow(modules$columns)
  
  if (no.mods==0) {
    stop("No modules")
  }

  if (!missing(data)) {
    if (nrow(data) != no.rows || ncol(data) != no.cols) {
      stop("`data' has a different dimension than `modules'")
    }
  }
  
  M <- lapply(to.plot, function(x) {
    if (binary) {
      outer(modules$rows[,x]!=0, modules$columns[,x]!=0)
    } else {
      outer(modules$rows[,x], modules$columns[,x])
    }
  })  
  
  if (!missing(data)) {
    no.mods <- no.mods + 1
    M <- c(list(data), M)
  }

  images(M, names=names, xlab=xlab, ylab=ylab, ...)
}

                       

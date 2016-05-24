# ++++++++++++++++++++++++++++++++++++++++++
# Generation of a landscape composed by regular polygones
# AB: 29/09/2009
# +++++++++++++++++++++++++++++++++++++++++
generPoly <- function(step=5, np=10, file="data", plot=TRUE ) {
  # -----------------------------------------------
  # FUNCTION
  #  Generate a regular grid of polygons
  # If file is not NULL:
  #  Generate an input file for CaliFlopp in format 2
  #  with  tabulate as separator
  # INPUT
  # step: step of the grid; >1
  # np: number of polygons on each axis
  # file: name of the file to build, or NULL
  # plot=T, when plotting is requested
  # OUTPUT
  # an object of class 'listpoly'
  # -----------------------------------------------
  if (step<=1)
    stop("step should be greater than 1")
  debut <- 1
  y <- c(debut, step, step,debut)
  xdeb <- c(debut,debut, step, step)
  npoly <- np*np
  if (!is.null(file)) {
    fic <- file(file, open="w")
    write(npoly, file=fic)
  }
  nvertices <- 4
  retour <- list()
#  centres<- list();  centres$x <- rep(NA, npoly);  centres$y <- rep(NA, npoly)
  ipoly <- 1
  for (i in 1:np) {
    x <- xdeb
    for (j in 1:np) {
      if (!is.null(file)) {
    tete <- c(ipoly, ipoly, nvertices)
    write(tete, file=fic, sep="\t", append=TRUE)
    write(x,  file=fic, sep="\t", append=TRUE)
    write(y,  file=fic, sep="\t", append=TRUE)
#    centres$x[ipoly] <- min(x) + (max(x)-min(x))/2-debut
#    centres$y[ipoly] <- min(y) + (max(y)-min(y))/2-debut
  } # end !is.null(file)
    a <- as.poly(x,y)
    retour[[ipoly]] <- a
    ipoly <- ipoly+1
    x <- x + step-1
  } #end j
    y <- y+step-1
  } # end i
  if (!is.null(file))
    close(fic)
  
    class(retour) <- "listpoly"
  if (plot==TRUE) {
    plot(retour)
#    pointmap(as.points(centres), add=T, axes=F, pch=as.character(1:npoly))
  }
  return(retour)
} # end generPoly

    

## functions starting with the dot are considered internal
## for this file and the one specific function that is exported from it.

## initialize the rgl window, clear etc.
.pca2d_plot_init <- function(pca.coords, components=c(1,2), show.plane=FALSE, axe.titles=NULL, bg="white", new=TRUE, biplot=FALSE) {

  # open new window if requested or not opened yet
  if(new) {
    cat("Creating new device\n")
    #open3d(antialias= 3)
  }

  prange <- list(x= c(-1, 1), y= c(-1, 1))

  if(biplot) {
    prange$x <- c(-1, 1) * max(abs(range(pca.coords[,1])))
    prange$y <- c(-1, 1) * max(abs(range(pca.coords[,2])))
  } else {
    prange$x <- range(pca.coords[,1])
    prange$y <- range(pca.coords[,2])
  }


  if(is.null(axe.titles)) {
    xlab <- paste("PC", components[1]) 
    ylab <- paste("PC", components[2])
  } else {
    xlab <- axe.titles[1] 
    ylab <- axe.titles[2] ;
  }

  plot(NULL, type= "n", 
    xlim= prange$x,
    ylim= prange$y,
    xlab= xlab,
    ylab= ylab,
    bty= "none"
    ) 


  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col= bg, border=NA)

  if(show.plane) {
    abline(h= 0, col= "grey")
    abline(v= 0, col= "grey")
  }

  usr
}

## Create a framed text
.text.frames <- function(coords, text, col= "black", bg= "#66666633") {

  sw <- strwidth(text) * 1.5
  sh <- strheight(text) * 1.7

  rect(
        coords[,1] - sw / 2, 
        coords[,2] - sh  / 2,
        coords[,1] + sw / 2, 
        coords[,2] + sh  / 2,
        col= bg,
        border= F
        
       )
  text(coords[,1], coords[,2], text, col= col)

}

## Add a biplot (variable loading plot)
.biplot.2D <- function(biplot.coords, biplot.vars) {

  r <- c(-1.1, 1.1) * max(abs(range(biplot.coords[,1:2])))
  par(usr= rep(r, 2))
  biplot.coords <- biplot.coords[ biplot.vars, ]

  nb <- nrow(biplot.coords)

  arrows(rep(0, nb), rep(0, nb), biplot.coords[,1], biplot.coords[,2], col= "red", lwd= 2, angle= 20, length= 0.1)
  labels <- rownames(biplot.coords)

  .text.frames(biplot.coords, labels)

}

## Plot ellipses surrounding each group
.group.ellipses.2D <- function(coords, group, group.col, ellipse.ci) {
  centr.coords <- calc.centroids(coords, group)
  lg <- levels(group)
  for(i in 1:length(lg)) {
    g   <- lg[i]
    sel <- group == g
    s   <- cov(coords[sel,,drop=FALSE])
    cc  <- centr.coords[i,]
    lines(ellipse(s, centre=cc, level=ellipse.ci), col=group.col[i])
  }
}

## Add group labels to the plot at the coordinates of centroids
.group.labels.2D <- function(coords, group, col= "black", ...) {
   centr.coords <- calc.centroids(coords, group)
   .text.frames(centr.coords, levels(group), col= col)
}

## Add centroids to the plot
.centroids.2D <- function(coords, group, group.shape, group.col, radius, col) {
  centr.coords <- calc.centroids(coords, group)
  points(centr.coords, col= group.col, pch= shape2pch(group.shape), cex = 2 * radius)
  segments(coords[,1], coords[,2], centr.coords[ group, 1], centr.coords[ group, 2 ], col= col)
}

## select coordinates from the pca object
.prepare.coords.2D <- function(pca, components) {
  pca.coords <- get.pca.coords(pca, components, 2)
  if(ncol(pca.coords) < 2) 
    stop(sprintf("Not enough components: found %d, need at least 2", ncol(pca.coords)))
  pca.coords <- pca.coords[ , components ]
  pca.coords
}

#' @rdname pca3d-package
#' @export
pca2d <-
function(pca, components=1:2, col=NULL, title=NULL, new=FALSE,
  axes.color ="black",
  bg="white",
  radius=1,
  group=NULL,
  shape=NULL,
  palette=NULL,
  fancy=FALSE,
  biplot=FALSE,
  biplot.vars=5,
  legend=NULL,
  show.scale=FALSE,
  show.labels=FALSE,
  labels.col="black",
  show.axes=TRUE,
  show.axe.titles=TRUE,
  axe.titles=NULL,
  show.plane=TRUE,
  show.shadows=FALSE,
  show.centroids=FALSE,
  show.group.labels=FALSE,
  show.ellipses=FALSE,
  ellipse.ci=0.95,
  ...
 ) {

  if(fancy) {
    show.labels       <- TRUE
    show.shadows      <- TRUE
    show.centroids    <- TRUE
    show.group.labels <- TRUE
  }

  if(missing(palette)) palette <- defaultPalettePCA3D()

  # prepare pca coordinates to plot
  pca.coords <- get.pca.coords(pca, 2, components)
  n.p <- nrow(pca.coords)

  if(!missing(biplot)) {
    biplot.coords <- get.biplot.coords(pca, biplot)
    biplot.coords <- biplot.coords[ , components]
    biplot <- TRUE
    biplot.vars <- get.biplot.vars(biplot.coords, biplot.vars)
  }

  # ----------------------------------
  # selecting shapes and colors
  if(missing(group)) 
    group <- rep("default", n.p)

  group <- factor(group)

  col <- autocol(n.p, col, group, palette, "col")
  group.col <- getGroupVals(col, group)

  shape <- autocol(n.p, shape, group, all.shapes, "shape")
  shape <- matchShapes(shape)
  group.shape <- getGroupVals(shape, group)

  ret <- print.legend(group, group.col, group.shape)
  # ----------------------------------

  # initialize the plot
  oldpar <- .pca2d_plot_init(pca.coords, components=components, 
    show.plane=show.plane, axe.titles=axe.titles, bg=bg, new=new,
    biplot=biplot) 
  on.exit(par(oldpar))

  # draw shadows ("lollipop" plot)
  if(show.shadows) 
    segments(pca.coords[,1], pca.coords[,2], pca.coords[,1], 0, col= "grey")

  # draw the points
  pch <- shape2pch(shape) ## "cube", "terahaedron" etc. are allowed
  points(pca.coords[,1], pca.coords[,2], col= col, pch= pch, cex= radius, ...)

  # add graphical extras
  if(show.centroids)    .centroids.2D(pca.coords, group, group.shape, group.col, 2 * radius, col)
  if(show.group.labels) .group.labels.2D(pca.coords, group, col= group.col)
  if(show.ellipses)     .group.ellipses.2D(pca.coords, group, group.col, ellipse.ci)

  # add legends for all points
  if(! missing(show.labels)) {
    if(class(show.labels) == "logical" & length(show.labels) == 1) {
      show.labels <- rownames(pca.coords)
    }

    text(pca.coords[,1], pca.coords[,2], show.labels, col= labels.col)
  }

  # add a biplot to the picture
  if(biplot) .biplot.2D(biplot.coords, biplot.vars)

  # show legend if it makes sense
  if(!is.null(legend) && !is.null(ret)) {
    legend(legend, 
      ret$groups,
      col=ret$colors,
      pch=shape2pch(ret$shapes),
      bg="white",
      box.lty=0
      )
  }


  return(invisible(ret))
}

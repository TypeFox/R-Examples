## functions starting with the dot are, by my convention, considered internal
## for this file and the one specific function that is exported from it.


## radius depends on the scale in each direction
.calc.radius <- function(radius, pca.coords) {
  radius <- rep(radius, 3)

  for(i in 1:3) {
    r <- 2*max(abs(pca.coords[ , i ]))
    radius[i] <- radius[i] * r / 100
  }
  print(radius)
  radius
}



## for drawing circles in 3D, precalculate some values
.sin.t <- sin(seq(0, 2 * pi, len= 10))
.cos.t <- cos(seq(0, 2 * pi, len= 10))


## draw labels for each centroid
.group.labels <- function(coords, group, ...) {
   centr.coords <- calc.centroids(coords, group)
   texts3d(centr.coords[,1], centr.coords[,2], centr.coords[,3], levels(group), col= "black", cex= 2, ...)
}

## add titles of axes
.show.axe.titles <- function(r, show.scale, axes.color, components, axe.titles) {
      if(missing(axe.titles) | is.null(axe.titles)) {
        axe.titles = paste("PC", components)
      }

      cex <- 1.2
      
      texts3d(
        c(cex *r[ 2,1 ], 0, 0),
        c(0, cex *r[ 2,2 ], 0),
        c(0, 0, cex *r[ 2,3 ]), col= "black", axe.titles) 

    
    if(show.scale) {
      axis3d('x',pos=c(NA, 0, 0), col= axes.color)
      axis3d('y',pos=c(0, NA, 0), col= axes.color)
      axis3d('z',pos=c(0, 0, NA), col= axes.color)
    }
}

## add the xy plane
.show.plane <- function(r) {
   quads3d(
     rbind(
       c(r[1,1],0,r[1,3]),
       c(r[1,1],0,r[2,3]),
       c(r[2,1],0,r[2,3]),
       c(r[2,1],0,r[1,3]) 
      ),
     col= "grey",
     alpha= 0.2
  ) 
}

## add the axes
.show.axes <- function(axes.color, ranges) {
  axes <- rbind(
    c(ranges[1,1], 0, 0),
    c(ranges[2,1], 0, 0),
    c(0, ranges[1,2], 0),
    c(0, ranges[2,2], 0),
    c(0, 0, ranges[1,3]),
    c(0, 0, ranges[2,3])
   )
  segments3d(axes, col= axes.color)

  radius <- 10

  scale <- c(1, 1, 1)
  if(! missing(ranges)) {
    scale <- ranges[2,]
    radius <- max(scale) / 50
  }

  arrows3d(axes, radius= radius, scale= scale,  col= axes.color)
}

## Plot ellipses surrounding each group
.group.ellipses.3D <- function(coords, group, group.col, ellipse.ci=0.95) {
  centr.coords <- calc.centroids(coords, group)
  lg <- levels(group)
  for(i in 1:length(lg)) {
    g   <- lg[i]
    sel <- group == g
    s   <- cov(coords[sel,,drop=FALSE])
    cc  <- centr.coords[i,]
    #lines(ellipse(s, centre=cc), col=group.col[i])
    shade3d(ellipse3d(s, centre=cc, level=ellipse.ci), col=group.col[i], alpha=0.2)
  }
}

## add group centroids
.centroids.3D <- function(coords, group, radius, group.shape, group.col, shape.functions) {

  centr.coords <- calc.centroids(coords, group)

  for(l in levels(group)) {
    shape.functions[[ group.shape[l] ]](centr.coords[l,,drop= F], col= group.col[l], radius= radius, alpha= 0.7)
    n.l <- length(which(group == l))
    tmp <- coords[ group == l, ]

    tmp <- rbind(tmp, t(sapply(1:nrow(tmp), function(x) centr.coords[ l, ])))
    tmp <- tmp[ as.vector(rbind(1:n.l, (n.l + 1):(2*n.l))), ]
    segments3d(tmp, col= group.col[ l ])
  }
}

## drop shadows on the xy plane
.draw.shadows.3D <- function(pca.coords, xz.functions, shape, radius) {

  seg <- NULL
  for(i in 1:nrow(pca.coords)) {
    tmp <- pca.coords[ i, ]
    seg <- rbind(seg, tmp)
    tmp[ 2 ] <- 0
    seg <- rbind(seg, tmp)
    #squarexz(tmp[1], tmp[2] + 0.001, tmp[3], radius= radius, col= "grey")
    xz.functions[[ shape[i] ]](tmp[1], 0 + 0.001 * radius[2], tmp[3], radius= radius, col= "grey")
    xz.functions[[ shape[i] ]](tmp[1], 0 - 0.001 * radius[2], tmp[3], radius= radius, col= "grey")
  }
  segments3d(seg, col= "grey", alpha= 0.5)
}

.squarexz <- function(x, y, z, radius= c(1, 1, 1), col= "grey", alpha= 0.5) {
  tmp <- NULL
  r <- 2*radius/3
  for(i in 1:length(x)) {
    tmp <- rbind(tmp,
      c(x[i] - r[1], y[i], z[i] - r[3]), 
      c(x[i] - r[1], y[i], z[i] + r[3]), 
      c(x[i] + r[1], y[i], z[i] - r[3]), 
      c(x[i] + r[1], y[i], z[i] + r[3]), 
      c(x[i] - r[1], y[i], z[i] + r[3]), 
      c(x[i] + r[1], y[i], z[i] - r[3]))
  }
  triangles3d(tmp, col= col, alpha= alpha)
}

.trianglexz <- function(x, y, z, radius= c(1, 1, 1), col= "grey", alpha= 0.5) {
  tmp <- NULL
  for(i in 1:length(x)) {
    tmp <- rbind(tmp, 
      c(x[i]                            , y[i], z[i] + radius[3]    ),
      c(x[i] + sqrt(3) * radius[1] / 2, y[i], z[i] - radius[3] / 2),
      c(x[i] - sqrt(3) * radius[1] / 2, y[i], z[i] - radius[3] / 2)
     )
  }
  triangles3d(tmp, col= col, alpha= alpha)

}

.circlexz <- function(x, y, z, radius= c(1, 1, 1), col= "grey", alpha= 0.5) {
  
  r <- radius * 2 / 3 
  n <- length(.sin.t)
  xv <- x + r[1] * .sin.t
  yv <- rep(y, n)
  zv <- z + r[3] * .cos.t

  tmp <- NULL
  for(i in 1:(n-1)) {
    tmp <- rbind(tmp, 
      c(x, y, z),
      c(xv[i],   yv[i],   zv[i]  ),
      c(xv[i+1], yv[i+1], zv[i+1]))
  }

  triangles3d(tmp, col= col, alpha= alpha)
}

.myspheres3d <- function(coords, radius= c(1, 1, 1), ...) {
    spheres3d(coords, radius= mean(radius), ...)
}

## create a 3D biplot
.biplot.3D <- function(biplot.coords, biplot.vars, ranges= NULL) {

  nvar <- length(biplot.vars)

  coords <- matrix(0, ncol= 3, nrow= nvar * 2)
  tips <- biplot.coords[ biplot.vars, ]
  coords[ seq(2, nvar * 2, by= 2), ] <- tips

  mx <- max(abs(coords))

  coords <- coords / mx 
  scale <- c(1, 1, 1)
  if(! missing(ranges)) {
    scale <- ranges[2,]
    coords <- t(apply(coords, 1, function(x) x * scale))
  }

  arrows3d(coords, col= "red", scale= scale)

  if(length(rownames(biplot.coords)) > 0) 
    names <- rownames(biplot.coords)[biplot.vars]
  else
    names <- biplot.vars

  coords <- coords[ seq(2, nvar * 2, by = 2), ]
  texts3d(coords[,1], coords[,2], coords[,3], names, col= "black")
}

.rgl.window.init <- function(bg) {
  cat("Creating new device\n")
  rgl.open()
  #open3d(antialias= 3)
  par3d(windowRect= 50 + c(0, 0, 640, 640))
  rgl.viewpoint(theta= 45, phi= 30, fov= 60, zoom= 1)
  bg3d(bg)
  rgl.pop("lights")
  light3d(theta=-45, phi=15, specular= "black")
  light3d(theta=55, phi=55,  specular= "black")
  #light3d(theta=15, phi=-15,  specular= "white")
}

## initialize the rgl window, clear etc.
.pca_plot_init3d <- function(new, show.plane, show.axes, show.axe.titles, show.scale, 
  axes.color, axe.titles, main.ranges, components, bg) {

  # open new window if requested or not opened yet
  if(new || rgl.cur() == 0) 
    .rgl.window.init(bg)

  par3d(skipRedraw= TRUE) 
  rgl.clear()

  #r <- apply(pca.coords, 2, function(x) c(-max(abs(range(x))), max(abs(range(x)))))
  if(show.plane) .show.plane(main.ranges)
  if(show.axes) .show.axes(axes.color, main.ranges * 1.1)
  if(show.axe.titles) .show.axe.titles(main.ranges, show.scale, axes.color, components, axe.titles)
  if(show.scale) {
    axis3d('x',pos=c(NA, 0, 0))
    axis3d('y',pos=c(0, NA, 0))
    axis3d('z',pos=c(0, 0, NA))
  }
}

## finalize the plot
.pca_plot_finish_3d <- function() {
  aspect3d(c(1, 1, 1))
  #aspect3d("iso")
  par3d(skipRedraw= F)
  rgl.bringtotop()
}


.plot.points <- function(shape, shape.functions, radius, col, pca.coords) {

  for(s in unique(shape)) {
    for(c in unique(col)) {
      sel <- which(shape == s & col == c)
      if(length(sel) > 0) {
        shape.functions[[s]](pca.coords[ sel,, drop=F ], 
          specular="#111111",
          col= as.character(c), alpha= 1, radius= radius, shininess=0 ) 
      }
    }
  }


}

.labels.3D <- function(show.labels, pca.coords, labels.col, n.p) {
  if(class(show.labels) == "logical" && length(show.labels) == 1) {
    show.labels <- rownames(pca.coords)
  }
  if(is.null(show.labels)) show.labels <- 1:n.p
  texts3d(pca.coords[,1:3], texts=show.labels, col=labels.col, adj=c(1, 1))
}


#' @rdname pca3d-package
#' @export
pca3d <- function(pca, components=1:3, col=NULL, title=NULL, new=FALSE,
  axes.color ="grey",
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
  show.shapes=TRUE,
  show.ellipses=FALSE,
  ellipse.ci=0.95
 ) {

  if(fancy) {
    show.labels       <- TRUE
    show.shadows      <- TRUE
    show.centroids    <- TRUE
    show.group.labels <- TRUE
  }

  if(missing(palette)) 
    palette <- defaultPalettePCA3D()

  # prepare pca coordinates to plot
  pca.coords <- get.pca.coords(pca, 3, components)
  n.p <- nrow(pca.coords)
  main.ranges <- apply(pca.coords, 2, function(x) c(-max(abs(range(x))), max(abs(range(x)))))

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

  # other things we need for drawing in 3D
  radius <- .calc.radius(radius, pca.coords)
  shape.functions <- list(sphere= .myspheres3d, tetrahedron= tetrahedrons3d, cube= cubes3d, octahedron=octahedrons3d)
  xz.functions    <- list(sphere= .circlexz,    tetrahedron= .trianglexz,     cube= .squarexz, octahedron= .squarexz)

  # initialize the plot
  .pca_plot_init3d(new, 
    show.plane, show.axes, show.axe.titles, show.scale, 
    axes.color, axe.titles, main.ranges, components, bg)

  # graphical add-ons
  if(show.shadows)   .draw.shadows.3D(pca.coords, xz.functions, shape, radius) # draw shadows ("lollipop" plot)
  if(show.centroids) .centroids.3D(pca.coords, group, 2 * radius, group.shape, group.col, shape.functions)
  if(show.group.labels) .group.labels(pca.coords, group)
  if(show.ellipses) .group.ellipses.3D(pca.coords, group, group.col, ellipse.ci)

  # draw the samples
  if(show.shapes) {
    .plot.points(shape, shape.functions, radius, col, pca.coords) 
  }

  # add the biplot
  if(biplot) 
    .biplot.3D(biplot.coords, biplot.vars, main.ranges)

  # label the samples
  if(!missing(show.labels)) 
    .labels.3D(show.labels, pca.coords, labels.col, n.p) 

  # finalize the plot
  .pca_plot_finish_3d() 
  
  # draw legend
  if(!is.null(legend) && !is.null(ret)) {
    legend3d(legend, 
      ret$groups,
      col=ret$colors,
      pch=ret$pch,
      box.lty=0
      )
  }

  return(invisible(ret))
}


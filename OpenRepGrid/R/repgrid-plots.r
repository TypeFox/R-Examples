
###############################################################################


#' Calculate coordinates for biplot.
#'
#' @param x             \code{repgrid} object.
#' @param g             Power of the singular value matrix assigned to the left singular 
#'                      vectors, i.e. the constructs.
#' @param h             Power of the singular value matrix assigned to the right singular 
#'                      vectors, i.e. the elements.
#' @param col.active    Columns (elements) that are no supplementary points, i.e. they are used
#'                      in the SVD to find principal components. default is to use all elements.
#' @param col.passive   Columns (elements) that are supplementary points, i.e. they are NOT used
#'                      in the SVD but projecte into the component space afterwards. They do not 
#'                      determine the solution. Default is \code{NA}, i.e. no elements are set 
#'                      supplementary.
#' @param ...           Parameters to be passed on to \code{center()} and \code{normalize}.       
#' @return              a \code{list}.
#'
#' @author        Mark Heckmann
#' @keywords internal
#' @export
#' 
calcBiplotCoords <- function(x, g=0, h=1-g, 
                             col.active=NA, 
                             col.passive=NA, 
                             ... ){
  # definition of active and passive (supplementary points) 
  if (!identical(col.active, NA) & !identical(col.passive, NA))
    stop("active OR passive columns must be defined")
  ne <- getNoOfElements(x)
  if (identical(col.active, NA)){                   # if no active points defined
    col.active <- seq_len(ne)                       # the rest is set active
    col.active <- setdiff(col.active, col.passive)
  } else if (identical(col.passive, NA)){           # if no passive points defined
    col.passive <- seq_len(ne)                      # the is set passive
    col.passive <- setdiff(col.passive, col.active)
  }
  
  X <- center(x, ...)           # center grid
  X <- normalize(X, ...)        # normalize grid
  
  X.active <- X[ , col.active]  # X with active columns (elements) only. Used for SVD. 
                                # The other supplementary elements are projected afterwards.

  dec <- svd(X.active)    # make SVD for reduced set of active points
  U <- dec$u              # left singular vector matrix
  D <- dec$d              # matrix of singular values
  V <- dec$v              # right singular vector matrix      
  
  # constructs coords
  C <- U %*% diag(D^g)    # standard form       
  # C <- X[, col.active] %*% V %*% (D^h)^-1 
  # C <- X[, col.active] %*% V %*% (D^(1-g))^-1 
  
  # element coords
  # E <- V %*% diag(D^h)                  # not used as supplementary points need to be calculated
  # t(X) %*% U %*% (D^g)^-1               # only works when g + h =1, thus: 
  E <- t(X) %*% U %*% diag((D^(1-h))^-1)  # only dependent on h not g                                            

  rownames(C) <- getConstructNames(x)[ ,2]  # names of direction into which vector points
  rownames(E) <- getElementNames(x)
    
  x@calcs$biplot <- list(X=X, element.coords=E, construct.coords=C,
                         D=D,U=U, V=V, col.passive=col.active, 
                         col.passive=col.passive)
  x
}


#' Map arbitrary numeric vector to a given range of values. 
#'
#' From a given numeric vector \code{z} the range is determined and 
#' the values are linearly mapped onto the interval 
#' given by \code{val.range}. This 
#' funtion can be used in order to map arbitrary vectors to a given
#' range of values.
#'
#' @param  z          numeric vector.
#' @param  val.range  numeric vector of lengths two (default \code{c(.5, 1)}).
#' @return numeric vector
#'
#' @author  Mark Heckmann
#' @keywords internal
#' @export
#'
mapCoordinatesToValue <- function(z, val.range=c(.5, 1)) {
	z.range <- c(min(z, na.rm=T), max(z, na.rm=T))
	slope <-  diff(val.range) / diff(z.range)
	int <- val.range[1] - z.range[1] * slope
	vals <- int + slope * z
  round(vals, 10)   # round at 10th digit to prevent values like 1.00000000001
}


#' Determine color values according to a given range of values. 
#'
#' From a given numeric vector z the range is determined and the values are 
#' linearly mapped onto the interval given by \code{val.range}. Then 
#' a color ramp using the colors given by \code{color} is created and 
#' the mapped values are transformed into hex color values. 
#'
#' @param  z          numeric vector.
#' @param  color      vector of length two giving color values \code{c("white", "black")}.
#' @param  val.range  numeric vector of lengths two (default \code{c(.2, .8)}).
#' @return numeric vector
#'
#' @author  Mark Heckmann
#' @keywords internal
#' @export
#'
mapCoordinatesToColor <- function(z, colors=c("white", "black"), val.range=c(.2,.8)){
  colorRamp <- makeStandardRangeColorRamp(colors)
  vals <- mapCoordinatesToValue(z, val.range)
  colorRamp(unlist(vals))     # unlist in case z comes as a data frame column
}


#' Coordinates of a sourrounding rectangle in direction of a given vector. 
#'
#' An arbitrary numeric vector in 2D is to be extented so it will 
#' end on the borders of a sourrounding rectangle of a given size.
#' Currently the vector is supposed to start in the origin \code{c(0,0)}.
#'
#' @param x      numeric vector of x coordinates x coordinates. 
#' @param y      numeric vector of y coordinates x coordinates. 
#' @param xmax   maximal x value for sourrounding rectangle (default is \code{1}).
#' @param ymax   maximal y value for sourrounding rectangle (default is \code{1}).
#' @param cx     center of retangle in x direction (not yet supported).
#' @param cy     center of retangle in x direction (not yet supported).
#' 
#' @return       a \code{dataframe} containing the x and y coordinates for the 
#'               extended vectors.
#'
#' @author    Mark Heckmann
#' @keywords internal
#' @export
#'
#' @examples \dontrun{
#'   calcCoordsBorders(1:10, 10:1)
#' 
#'   x <- c(-100:0, 0:100, 100:0, 0:-100)/10
#'   y <- c(0:100, 100:0, -(0:100), -(100:0))/10
#'   xy1 <- calcCoordsBorders(x, y)
#'   xy2 <- calcCoordsBorders(x, y, xm=1.2, ym=1.2)
#'   plot(xy2[,1], xy2[,2], type="n")
#'   segments(xy1[,1],xy1[,2],xy2[,1], xy2[,2])
#' }
#'
calcCoordsBorders <- function(x, y, xmax=1, ymax=1, cx=0, cy=0)
{
  is.lr.part <- abs(x*ymax/xmax) >= abs(y)        # which are left and right parts

  # left and right part             
  sign.x <- sign(x)             # positive or negative value    
  sign.x[sign.x == 0] <- 1      # zeros in posistive direction
  a.lr <- xmax * sign(x)        # x is fix on the left and right side
  b.lr <- y/x * a.lr            
  
  # upper and lower part
  sign.y <- sign(y)
  sign.y[sign.y == 0] <- 1
  b.ul <- ymax * sign(y)
  a.ul <- x/y * b.ul
  
  a.lr <- unlist(a.lr)
  b.lr <- unlist(b.lr)
  a.ul <- unlist(a.ul)
  b.ul <- unlist(b.ul)
  
  a.lr[is.nan(a.lr)] <- 0       # in case one of x or y is zero Inf results ans subsequently NaN
  b.lr[is.nan(b.lr)] <- 0
  a.ul[is.nan(a.ul)] <- 0
  b.ul[is.nan(b.ul)] <- 0  
  
  # join both parts
  b <- (b.ul * !is.lr.part) + (b.lr * is.lr.part)
  a <- (a.ul * !is.lr.part) + (a.lr * is.lr.part)
  a[abs(a) > xmax] <- (xmax * sign(a))[abs(a) > xmax]
  b[abs(b) > ymax] <- (ymax * sign(b))[abs(b) > ymax]
  
  data.frame(x=a, y=b)
}


# calculate the coords for the label rectangle
#
# TODO: supply x.ext in mm and convert to usr coords
#
# @param  xy        \code{dataframe} with x and y coords.
# @param  labels    vector of strings.
# @param  cex       vector of cex values (default is \code{.7}).
# @param  x.ext     scalar giving the horizontal margin 
#                   of the rectangle in NDC coordinates
#                   (default is \code{.02}).
# @param  y.ext     scalar giving the vertical margin 
#                   of the rectangle in NDC coordinates
#                   (default is \code{.02}).
# @return \code{dataframe} with coordinates for the lower left and 
#         upper right rectangle borders (\code{x0, y0, x1, y1}).
#
calcRectanglesCoordsForLabels <- function(xy, labels, cex=.7, 
                                          x.ext=.02, y.ext=.02){
  if (length(cex) == 1)
    cex <- rep(cex, dim(xy)[1])
  
  heights <- vector()
  widths <- vector()
  
  for (i in 1:dim(xy)[1]){
    heights[i] <- strheight(labels[i], cex=cex[i])   # determine necessary height for text 
    widths[i] <- strwidth(labels[i], cex=cex[i])     # determine necessary width for text 
  }
  # make adj adjustements
  leftSide <- xy[, 1] < 0
  labelsBorders <- data.frame(x0= xy[, 1] - (widths * leftSide), 
                              y0= xy[, 2] - heights/2, 
                              x1= xy[, 1] + (widths * !leftSide), 
                              y1= xy[, 2] + heights/2)
  # extend borders for neater look
  labelsBorders$x0 <- labelsBorders$x0 - x.ext
  labelsBorders$y0 <- labelsBorders$y0 - y.ext
  labelsBorders$x1 <- labelsBorders$x1 + x.ext
  labelsBorders$y1 <- labelsBorders$y1 + y.ext
  
  labelsBorders
}


#' Detect if two rectangles overlap. 
#'
#' The overlap is assessed in x AND y.
#'
#' @param  a   vector with four coordinates \code{c(x0,y0,x1,y1)}.
#' @param  b   vector with four coordinates \code{c(x0,y0,x1,y1)}.
#' @return     \code{logical}. TRUE if rectangles overlap.
#'
#' @keywords internal
#' @export
#'
#' @examples \dontrun{
#'   #overlap in x and y
#'   a <- c(0,0,2,2)
#'   b <- c(1,1,4,3)
#'   plot(c(a,b), c(a,b), type="n")
#'   rect(a[1], a[2], a[3], a[4])
#'   rect(b[1], b[2], b[3], b[4])
#'   doRectanglesOverlap(a,b)
#' 
#'   # b contained in a vertically
#'   a <- c(5,0,20,20)
#'   b <- c(0, 5,15,15)
#'   plot(c(a,b), c(a,b), type="n")
#'   rect(a[1], a[2], a[3], a[4])
#'   rect(b[1], b[2], b[3], b[4])
#'   doRectanglesOverlap(a,b)
#' 
#'   # overlap only in y
#'   a <- c(0,0,2,2)
#'   b <- c(2.1,1,4,3)
#'   plot(c(a,b), c(a,b), type="n")
#'   rect(a[1], a[2], a[3], a[4])
#'   rect(b[1], b[2], b[3], b[4])
#'   doRectanglesOverlap(a,b)
#' }
#'
doRectanglesOverlap <- function(a, b, margin=0){
  overlap1D <- function(a0, a1, b0, b1){ # overlap if one of four conditions is satisfied
    (a0 <= b1 & b1 <= a1) |     # b overlaps at bottom
    (a0 <= b0 & b0 <= a1) |     # b overlaps at top
    (a0 >= b0 & a1 <= b1) |     # b overlaps at bottom and top
    (a0 <= b0 & a1 >= b1)       # b contained within a
  }
  overlap.x <- overlap1D(a[1], a[3], b[1], b[3])  # overlap in x ?
  overlap.y <- overlap1D(a[2], a[4], b[2], b[4])  # overlap in y ?
  as.logical(overlap.x & overlap.y)               # overlap in x and y, strip off vector names ?
}


# calculate angle between vector and x-y plane
# a   vector
# n   plane normal vector
degreesBetweenVectorAndPlane <- function(a, n){
    rad <- asin( abs(n %*% a) / 
                (sum(n^2)^.5 * sum(a^2)^.5))
    rad * 180/pi                # convert from radians to degrees
}


#' A graphically unsophisticated version of a biplot.
#'
#' It will draw elements and constructs vectors using similar
#' arguments as \code{\link{biplot2d}}. It is a version for quick 
#' exploration used during development.
#' 
#' @param x             \code{repgrid} object.
#' @param dim           Dimensions (i.e. principal components) to be used for biplot 
#'                      (default is \code{c(1,2)}).
#' @param  center		    Numeric. The type of centering to be performed. 
#'                      0= no centering, 1= row mean centering (construct), 
#'                      2= column mean centering (elements), 3= double-centering (construct and element means),
#'                      4= midpoint centering of rows (constructs).
#'                      The default is \code{1} (row centering).
#' @param normalize     A numeric value indicating along what direction (rows, columns)
#'                      to normalize by standard deviations. \code{0 = none, 1= rows, 2 = columns}
#'                      (default is \code{0}).
#' @param g             Power of the singular value matrix assigned to the left singular 
#'                      vectors, i.e. the constructs.
#' @param h             Power of the singular value matrix assigned to the right singular 
#'                      vectors, i.e. the elements.
#' @param col.active    Columns (elements) that are no supplementary points, i.e. they are used
#'                      in the SVD to find principal components. default is to use all elements.
#' @param col.passive   Columns (elements) that are supplementary points, i.e. they are NOT used
#'                      in the SVD but projecte into the component space afterwards. They do not 
#'                      determine the solution. Default is \code{NA}, i.e. no elements are set 
#'                      supplementary.
#' @param unity         Scale elements and constructs coordinates to unit scale in 2D (maximum of 1)
#'                      so they are printed more neatly (default \code{TRUE}).
#' @param zoom          Scaling factor for all vectors. Can be used to zoom
#'                      the plot in and out (default \code{1}).
#' @param scale.e       Scaling factor for element vectors. Will cause element points to move a bit more
#'                      to the center. This argument is for visual appeal only.
#' @param e.point.col   Color of the element symbols (default is \code{"black"}.
#' @param e.point.cex   Size of the element symbol (default is \code{1}.
#' @param e.label.col   Color of the element labels (default is \code{"black"}.
#' @param e.label.cex   Size of the element labels (default is \code{.7}.
#' @param c.point.col   Color of the construct lines (default is \code{grey(.6)}.
#' @param c.label.col   Color of the construct labels (default is \code{grey(.6)}.
#' @param c.label.cex   Size of the costruct labels (default is \code{.6}.
#' @param ...           Parameters to be passed on to \code{center()} and \code{normalize}.       
#' @return              \code{repgrid} object.
#'
#' @author Mark Heckmann
#' @export
#' 
#' @seealso   Unsophisticated biplot: \code{\link{biplotSimple}}; \cr
#'            2D biplots:
#'            \code{\link{biplot2d}},
#'            \code{\link{biplotEsa2d}},
#'            \code{\link{biplotSlater2d}};\cr
#'            Pseudo 3D biplots:
#'            \code{\link{biplotPseudo3d}},  
#'            \code{\link{biplotEsaPseudo3d}},
#'            \code{\link{biplotSlaterPseudo3d}};\cr
#'            Interactive 3D biplots:
#'            \code{\link{biplot3d}},
#'            \code{\link{biplotEsa3d}},
#'            \code{\link{biplotSlater3d}};\cr
#'            Function to set view in 3D:
#'            \code{\link{home}}.
#'
#' @examples \dontrun{
#'    
#'    biplotSimple(boeker)
#'    biplotSimple(boeker, unity=F)
#'
#'    biplotSimple(boeker, g=1, h=1)              # INGRID biplot
#'    biplotSimple(boeker, g=1, h=1, center=4)    # ESA biplot
#'
#'    biplotSimple(boeker, zoom=.9)               # zooming out
#'    biplotSimple(boeker, scale.e=.6)            # scale element vectors
#'
#'    biplotSimple(boeker, e.point.col="brown")   # change colors
#'    biplotSimple(boeker, e.point.col="brown",
#'                 c.label.col="darkblue")
#' }
#'
biplotSimple <- function(x, dim=1:2, center=1, normalize=0, 
                          g=0, h=1-g, unity=T,
                          col.active=NA, 
                          col.passive=NA, 
                          scale.e=.9, zoom=1, 
                          e.point.col="black", 
                          e.point.cex=1,
                          e.label.col="black",
                          e.label.cex=.7, 
                          c.point.col=grey(.6),
                          c.label.col=grey(.6),
                          c.label.cex=.6,
                          ...){
  par(mar=c(1,1,1,1))
  d1 <- dim[1]
  d2 <- dim[2]
  
  x <- calcBiplotCoords(x, g=g, h=h, center=center, 
                        normalize=normalize, 
                        col.active=col.active, 
                        col.passive=col.passive, ...)
  cnames <- getConstructNames(x)
  E <- x@calcs$biplot$el
  C <- x@calcs$biplot$con
  X <- x@calcs$biplot$X

  max.e <- max(abs(E[ ,dim]))
  max.c <- max(abs(C[ ,dim]))
  mv <- max(max.e, max.c)
  if (unity){
    max.e <- max(apply(E[ ,dim[1:2]]^2, 1, sum)^.5)   # maximal length of element vectors
    max.c <- max(apply(C[ ,dim[1:2]]^2, 1, sum)^.5)   # maximal length of construct vectors
    se <- 1/max.e  * scale.e   # scale to unity to make E and C same size
    sc <- 1/max.c
  } else {
    se <- 1
    sc <- 1
  }
  Cu <- C * sc
  Eu <- E * se
  
  mv <- max(abs(rbind(Cu, Eu)))
  Cu <- Cu * zoom
  Eu <- Eu * zoom
  
  # make biplot
  plot(0, xlim=c(-mv, mv), ylim=c(-mv, mv), type="n", asp=1,
        xaxt="n", yaxt="n", xaxs="i", yaxs="i")
  abline(v=0, h=0, col="grey")
  
  # plot constructs and labels
  arrows(0,0, -Cu[ ,d1], -Cu[ ,d2], length=.05, 
         col=c.point.col, lty=1)              # plot left poles
  text(-Cu[ ,d1], -Cu[ ,d2], cnames[,1], pos=1, 
        cex=c.label.cex, col=c.label.col)
  arrows(0,0, Cu[ ,d1], Cu[ ,d2], length=.05, 
         col=c.point.col, lty=3)               # plot right poles
  text(Cu[ ,d1], Cu[ ,d2], cnames[,2], pos=1, 
       cex=c.label.cex, col=c.label.col)
  
  # plot elements and labels
  points(Eu[, dim], pch=15, col=e.point.col, cex=e.point.cex)           # plot elements
  text(Eu[, dim], labels=rownames(Eu), cex=e.label.cex,
        col=e.label.col, pos=2)                                      # label elements  
  invisible(x)
}

# library(xtable)
# x <- biplotSimple(raeithel, dim=1:2, g=1, h=1, col.act=c(1,2,3,5,10,12))
# ssq.table <- ssq(x)
# #ssq.table[ssq.table < 10] <- NA
# res <- xtable(round(ssq.table, 1), digits=1, 
#               align=c("l", rep("r", ncol(ssq.table))), caption="Percentage of element's SSQ explained")
# print(res, table.placement="H", hline.after=c(-1,0,nrow(ssq.table)-1, nrow(ssq.table)))



#' Prepare dataframe passed to drawing functions for biplots.
#'
#' Data frame contains the variables \code{type, show, x, y, 
#'  z, labels, color, cex}.
#'
#' @param x             \code{repgrid} object.
#' @param dim           Dimensions to be used for biplot (default is \code{c(1,2)}).
#' @param map.dim       Third dimension used to map aesthetic attributes (depth)
#'                      (default is \code{3}).
#' @param e.point.col   Color(s) of the element symbols. Two values can be entered that will
#'                      create a color ramp. The values of \code{map.dim} are mapped onto the ramp.
#'                      The default is \code{c("white", "black")}. If only one color color value
#'                      is supplied (e.g. \code{"black"}) no mapping occurs and all elements
#'                      will have the same color irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param e.point.cex   Size of the element symbols. Two values can be entered that will
#'                      represents the lower and upper size of a range of cex the values of \code{map.dim} 
#'                      are mapped onto. The default is \code{c(.4, .8)}. If only one cex value
#'                      is supplied (e.g. \code{.7}) no mapping occurs and all elements
#'                      will have the same size irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param e.label.col   Color(s) of the element labels. Two values can be entered that will
#'                      create a color ramp. The values of \code{map.dim} are mapped onto the ramp.
#'                      The default is \code{c("white", "black")}. If only one color color value
#'                      is supplied (e.g. \code{"black"}) no mapping occurs and all element labels
#'                      will have the same color irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param e.label.cex   Size of the element labels. Two values can be entered that will
#'                      represents the lower and upper size of a range of cex the values of \code{map.dim} 
#'                      are mapped onto. The default is \code{c(.4, .8)}. If only one cex value
#'                      is supplied (e.g. \code{.7}) no mapping occurs and all element labels
#'                      will have the same size irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param e.color.map   Value range to determine what range of the color ramp defined in 
#'                      \code{e.color} will be used for mapping the colors. 
#'                      Default is \code{c(.4, ,1)}. Usually not important for the user. 
#' @param c.point.col   Color(s) of the construct symbols. Two values can be entered that will
#'                      create a color ramp. The values of \code{map.dim} are mapped onto the ramp.
#'                      The default is \code{c("white", "darkred")}. If only one color color value
#'                      is supplied (e.g. \code{"black"}) no mapping occurs and all elements
#'                      will have the same color irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param c.point.cex   Size of the construct symbols. Two values can be entered that will
#'                      represents the lower and upper size of a range of cex the values of \code{map.dim} 
#'                      are mapped onto. The default is \code{c(.4, .8)}. If only one cex value
#'                      is supplied (e.g. \code{.7}) no mapping occurs and all elements
#'                      will have the same size irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param c.label.col   Color(s) of the construct labels. Two values can be entered that will
#'                      create a color ramp. The values of \code{map.dim} are mapped onto the ramp.
#'                      The default is \code{c("white", "black")}. If only one color color value
#'                      is supplied (e.g. \code{"black"}) no mapping occurs and all construct labels
#'                      will have the same color irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param c.label.cex   Size of the construct labels. Two values can be entered that will
#'                      represents the lower and upper size of a range of cex the values of \code{map.dim} 
#'                      are mapped onto. The default is \code{c(.4, .8)}. If only one cex value
#'                      is supplied (e.g. \code{.7}) no mapping occurs and all construct labels
#'                      will have the same size irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param c.color.map   Value range to determine what range of the color ramp defined in 
#'                      \code{c.color} will be used for mapping. Default is \code{c(.4, ,1)}.
#'                      Usually not important for the user. 
#' @param devangle      The deviation angle from the x-y plane in degrees. These can only be calculated
#'                      if a third dimension \code{map.dim} is specified. Only the constructs 
#'                      vectors that do not depart 
#'                      more than the specified degrees from the shown x-y plane will be printed.
#'                      This facilitates the visual interpretation, as only vectors represented in
#'                      the current plane are shown. Set the value to \code{91} (default) 
#'                      to show all vectors.
#' @param unity         Scale elements and constructs coordinates to unit scale in 2D (maximum of 1)
#'                      so they are printed more neatly (default \code{TRUE}).
#' @param unity3d       Scale elements and constructs coordinates to unit scale in 3D (maximum of 1)
#'                      so they are printed more neatly (default \code{TRUE}).
#' @param scale.e       Scaling factor for element vectors. Will cause element points to move a bit more
#'                      to the center. This argument is for visual appeal only.
#' @param ...           Not evaluated. 
#'
#' @return              \code{dataframe} containing the variables \code{type, show, x, y, 
#'                      z, labels, color, cex}. Usually not of interest to the user.
#' @note                TODO:  How to omit \code{map.dim}?
#'
#' @author  Mark Heckmann
#' @keywords internal
#' @export
#'
prepareBiplotData <- function(x, dim=c(1,2), map.dim=3, 
                               #e.color=c("white", "black"), 
                               #c.color=c("white", "darkred"),
                               e.label.cex=.8,
                               c.label.cex=.6,
                               e.label.col="black",
                               c.label.col=grey(.8),
                               e.point.cex=.7,
                               c.point.cex=.8,
                               e.point.col="black",
                               c.point.col="darkred",
                               #e.cex.map=c(.6, .8),
                               #c.cex.map=c(.6, .8),
                               e.color.map=c(.4, 1),
                               c.color.map=c(.4, 1),
                               c.points.devangle=90,
                               c.labels.devangle=90,
                               c.points.show=TRUE,
                               c.labels.show=TRUE,
                               e.points.show=TRUE,
                               e.labels.show=TRUE, 
                               unity=TRUE, 
                               unity3d=FALSE, 
                               scale.e=.9, 
                               ...)
{
  dim <- c(dim, map.dim)

  # make vector of length two if only one color/cex is specified
  #if (length(e.color) == 1)         # element color
  #  e.color <- rep(e.color, 2)          
  #if (length(c.color) == 1)         # construct color
  #  c.color <- rep(c.color, 2)
  #if (length(e.cex.map) == 1)       # element cex for pseudo 3d dimension
  #  e.cex.map <- rep(e.cex.map, 2)
  #if (length(c.cex.map) == 1)       # construct cex for pseudo 3d dimension
  #  c.cex.map <- rep(c.cex.map, 2)
  
  if (length(e.label.col) == 1)     # label color(s) for elements 
    e.label.col <- rep(e.label.col, 2)
  if (length(c.label.col) == 1)     # label color(s) for constructs 
    c.label.col <- rep(c.label.col, 2)  
  if (length(e.point.col) == 1)     # point color(s) for elements 
    e.point.col <- rep(e.point.col, 2)
  if (length(c.point.col) == 1)     # point color(s) for constructs 
    c.point.col <- rep(c.point.col, 2)
  
  if (length(e.label.cex) == 1)     # label cex(s) for elements 
    e.label.cex <- rep(e.label.cex, 2)
  if (length(c.label.cex) == 1)     # label cex(s) for constructs 
    c.label.cex <- rep(c.label.cex, 2)  
  if (length(e.point.cex) == 1)     # point cex(s) for elements 
    e.point.cex <- rep(e.point.cex, 2)
  if (length(c.point.cex) == 1)     # point cex(s) for constructs 
    c.point.cex <- rep(c.point.cex, 2)
  
    if (length(e.color.map) == 1)     # element color for pseudo 3d dimension
      e.color.map <- rep(e.color.map, 2)
    if (length(c.color.map) == 1)     # construct color for pseudo 3d dimension
      c.color.map <- rep(c.color.map, 2)
      
  # construct data frame containing all information needed for different plotting functions
  # (e.g. rgl and biplot functions)
  labels.e <- getElementNames(x)
  labels.cl <- getConstructNames(x)[,1] 
  labels.cr <- getConstructNames(x)[,2]
  labels.all <- c(labels.e, labels.cr, labels.cl)         # join all labels
  type <- factor(c(rep("e", getNoOfElements(x)),          # make factor specifying if row is element or construct
                   rep(c("cl", "cr"), each=getNoOfConstructs(x))))
  df <- data.frame(type=type, label=labels.all, stringsAsFactors=FALSE)
  df$cex <- .7            # default cex
  df$showpoint <- T       # default value for show point
  df$showlabel <- T       # default value for show label
  df$color <- grey(0)     # default color
  df$label.col <- "darkgreen"     # default label color
  df$point.col <- "purple"    # default point color
  df$label.cex <- .7          # default label size
  df$point.cex <- .7          # default point size
  
  # calculate and add coordinates
  #x <- calcBiplotCoords(x, ...)
   
  E <- x@calcs$biplot$el
  C <- x@calcs$biplot$con
  X <- x@calcs$biplot$X

  # scale to unity to make E and C same size.
  # Two types of unity, for 2D and 3D
  
  if (unity){
    max.e <- max(apply(E[ ,dim[1:2]]^2, 1, sum)^.5)   # maximal length of element vectors
    max.c <- max(apply(C[ ,dim[1:2]]^2, 1, sum)^.5)   # maximal length of construct vectors
    se <- 1/max.e  * scale.e   # scale to unity to make E and C same size
    sc <- 1/max.c   
  } 
  if (unity3d){
    #max.e <- max(abs(E[ ,dim[1:3]]), na.rm=T)    
    #max.c <- max(abs(C[ ,dim[1:3]]), na.rm=T)    
    max.e <- max(apply(E[ ,dim[1:3]]^2, 1, sum)^.5)   # maximal length of element vectors
    max.c <- max(apply(C[ ,dim[1:3]]^2, 1, sum)^.5)   # maximal length of construct vectors
    se <- 1/max.e * scale.e     # scale to unity to make E and C same size
    sc <- 1/max.c     
  } 
  if (!unity & !unity3d){
    se <- 1
    sc <- 1
  }
  Cu <- C * sc
  Eu <- E * se
   
  coords <- rbind(Eu[, dim], Cu[ ,dim], -Cu[ ,dim])
  colnames(coords) <- c("x", "y", "z")
  rownames(coords) <- NULL            # otherwise warning in cbind occurs
  df <- cbind(df, coords) #, check.rows=F)
  if (is.na(dim[3]))   # if no 3rd dimension in specified, all values are set to zero i.e. neutral
    df$z <- 0
    
  # plot coords for all points
  z <- subset(df, type=="e", sel=z)                     # z scores for elements
  #cex.e <- mapCoordinatesToValue(z, e.cex.map)
  cex.label.e <- mapCoordinatesToValue(z, e.label.cex)
  cex.point.e <- mapCoordinatesToValue(z, e.point.cex)
  #color.e <- mapCoordinatesToColor(z, color=e.color, val.range=e.color.map)
  color.label.e <- mapCoordinatesToColor(z, colors=e.label.col, val.range=e.color.map)
  color.point.e <- mapCoordinatesToColor(z, colors=e.point.col, val.range=e.color.map)
  
  z <- subset(df, type=="cl", sel=z)
  #cex.cl <- mapCoordinatesToValue(z, c.cex.map)
  cex.label.cl <- mapCoordinatesToValue(z, c.label.cex)
  cex.point.cl <- mapCoordinatesToValue(z, c.point.cex)  
  #color.cl <- mapCoordinatesToColor(z, color=c.color, val.range=c.color.map)
  color.label.cl <- mapCoordinatesToColor(z, colors=c.label.col, val.range=c.color.map)
  color.point.cl <- mapCoordinatesToColor(z, colors=c.point.col, val.range=c.color.map)

  z <- subset(df, type=="cr", sel=z)
  #cex.cr <- mapCoordinatesToValue(z, c.cex.map)
  cex.label.cr <- mapCoordinatesToValue(z, c.label.cex)
  cex.point.cr <- mapCoordinatesToValue(z, c.point.cex)
  #color.cr <- mapCoordinatesToColor(z, color=c.color, val.range=c.color.map)
  color.label.cr <- mapCoordinatesToColor(z, colors=c.label.col, val.range=c.color.map)
  color.point.cr <- mapCoordinatesToColor(z, colors=c.point.col, val.range=c.color.map)

  #df$cex <- unlist(rbind(cex.e, cex.cl, cex.cr))
  #df$color <- c(color.e, color.cl, color.cr)
  df$label.col <-  c(color.label.e, color.label.cl, color.label.cr)
  df$point.col <-  c(color.point.e, color.point.cl, color.point.cr)
  df$label.cex <-  unlist(c(cex.label.e, cex.label.cl, cex.label.cr))
  df$point.cex <-  unlist(c(cex.point.e, cex.point.cl, cex.point.cr))
  df$devangle <- apply(df, 1, function(x) { 
      a <- as.numeric( c(x["x"], x["y"], x["z"]) )
      n <- c(0,0,1)                              # normal vector for x-y plane
      degreesBetweenVectorAndPlane(a=a, n=n)
  })
  
  # calculate absolute deviation angle from shown plane. If it is bigger than given values
  # the constructs will not be shown on the side and/or the construct points will
  # not be printed. If values >=90 all strokes and points are shown.  

  cs <- subset(df, type %in% c("cl", "cr"))
  draw <- abs(cs$devangle) <= c.labels.devangle   # which angles are smaller or equal than the maximal allowed ones?
  cs$showlabel <- cs$showlabel & draw             # show only labels that are requested and within allowed angle range
  draw <- abs(cs$devangle) <= c.points.devangle   # which angles are smaller or equal than the maximal allowed ones?
  cs$showpoint <- cs$showpoint & draw             # show only labels that are requested and within allowed angle range
  df[df$type %in% c("cl", "cr"), ] <- cs          # show only labels that are requested and within allowed angle range

  # elements #
  # select which element labels to show
  # numerical values for element selection are converted to logical
  seq.e <- seq_len(getNoOfElements(x))
  if (! (identical(e.labels.show, T) | identical(e.labels.show, F) | all(is.numeric(e.labels.show))) )
    stop("'e.labels.show' must either be a logical value or a numeric vector")
  if (all(is.numeric(e.labels.show)))
    e.labels.show <- seq.e %in% seq.e[e.labels.show]  
  df[df$type == "e", "showlabel"] <- e.labels.show   # replace showlabel column for elements

  # select which element points to show
  # numerical values for element selection are converted to logical
  if (! (identical(e.points.show, T) | identical(e.points.show, F) | all(is.numeric(e.points.show))) )
    stop("'e.points.show' must either be a logical value or a numeric vector")
  if (all(is.numeric(e.points.show)))
    e.points.show <- seq.e %in% seq.e[e.points.show]  
  df[df$type == "e", "showpoint"] <- e.points.show     # replace showpoint column for elements

  # constructs #  TODO: mechanism fill fail for single / double mode grids
  # select which construct labels to show (independently from devangle)
  # numerical values for construct selection are converted to logical
  seq.c <- seq_len(getNoOfConstructs(x))    # TODO for single mode grids
  if (! (identical(c.labels.show, T) | identical(c.labels.show, F) | all(is.numeric(c.labels.show))) )
   stop("'c.labels.show' must either be a logical value or a numeric vector")
  if (all(is.numeric(c.labels.show))){
    doubleadd <- c.labels.show + sign(c.labels.show[1]) * getNoOfConstructs(x)  # if double mode
    c.labels.show <- seq.c %in% seq.c[c(c.labels.show, doubleadd)]       
  }
  show.tmp <- df[df$type %in% c("cl", "cr"), "showlabel"] 
  df[df$type %in% c("cl", "cr"), "showlabel"] <- c.labels.show & show.tmp       # replace showlabel column for elements

  # select which construct points to show (independently from devangle)
  # numerical values for construct selection are converted to logical
  if (! (identical(c.points.show, T) | identical(c.points.show, F) | all(is.numeric(c.points.show))) )
    stop("'c.points.show' must either be a logical value or a numeric vector")
  if (all(is.numeric(c.points.show)))
    c.points.show <- seq.c %in% seq.c[c.points.show]  
  points.tmp <- df[df$type %in% c("cl", "cr"), "showpoint"]
  df[df$type %in% c("cl", "cr"), "showpoint"] <- c.points.show &  points.tmp    # replace showpoint column for elements
  
  #list(rg=x, df=df)
  x@calcs$biplot$e.unity <- Eu
  x@calcs$biplot$c.unity <- Cu
  x@plotdata <- df
  x
}


#' biplotDraw is the workhorse doing the drawing of a 2D biplot. 
#'
#' When the number of elements and constructs differs to a large extent, the 
#' absolute values of the coordinates for either constructs or elements 
#' will be much smaller or greater. This is an inherent property of the biplot.
#' In the case it is not necessary to be able to read off the original 
#' data entries from the plot, the axes for elements and constructs
#' can be scaled seperately. The proportional projection values will 
#' stay unaffetced. the absolue will change though. For grid interpretation 
#' the absolze values are usually oh no importance. Thus, there is an optional
#' argument \code{normalize} which is \code{FALSE} as a default which
#' rescales the axes so the longest construct and element vector will be 
#' set to the length of \code{1}.
#' 
#' @param x                    \code{repgrid} object.
#' @param inner.positioning    Logical. Whether to calculate positions to minimize overplotting of 
#'                             elements and construct labels (default is\code{TRUE}). Note that
#'                             the positioning may slow down the plotting.
#' @param outer.positioning    Logical. Whether to calculate positions to minimize overplotting of 
#'                             of construct labels on the outer borders (default is\code{TRUE}). Note that
#'                             the positioning may slow down the plotting.
#' @param c.labels.inside      Logical. Whether to print construct labels next to the points.
#'                             Can be useful during inspection of the plot (default \code{FALSE}).
#' @param flipaxes             Logical vector of length two. Whether x and y axes are reversed 
#'                             (default is \code{c(F,F)}).
#' @param strokes.x            Length of outer strokes in x direction in NDC.  
#' @param strokes.y            Length of outer strokes in y direction in NDC.
#' @param offsetting           Do offsetting? (TODO)
#' @param offset.labels        Offsetting parameter for labels (TODO).
#' @param offset.e             offsetting parameter for elements (TODO).
#' @param axis.ext             Axis extension factor (default is \code{.1}). A bigger value will 
#'                             zoom out the plot.
#' @param mai                  Margins available for plotting the labels in inch 
#'                             (default is \code{c(.2, 1.5, .2, 1.5)}).
#' @param rect.margins         Vector of length two (default is \code{c(.07, .07)}). Two values
#'                             specifying the additional horizontal and vertical margin around each 
#'                             label.      
#' @param srt                  Angle to rotate construct label text. Only used in case \code{offsetting=FALSE}.
#' @param cex.pos              Cex parameter used during positioning of labels if prompted. Does
#'                             usually not have to be changed by user.
#' @param xpd                  Logical (default is \code{TRUE}). Wether to extend text labels 
#'                             over figure region. Usually not needed by the user.      
#' @param c.lines              Logical. Whether construct lines from the center of the biplot
#'                             to the sourrounding box are drawn (default is \code{FALSE}).
#' @param col.c.lines          The color of the construct lines from the center to the borders 
#'                             of the plot (default is \code{gray(.9)}).
#' @param zoom                 Scaling factor for all vectors. Can be used to zoom
#'                             the plot in and out (default \code{1}).
#' @param ...                  Not evaluated.
#' @return                     Invisible return of dataframe used during construction of plot 
#'                             (useful for developers).
#'
#' @author  Mark Heckmann
#' @export
#' @keywords internal
#'
biplotDraw <- function(x, 
                      inner.positioning=TRUE,
                      outer.positioning=TRUE,
                      c.labels.inside=F,
                      flipaxes=c(F,F), 
                      strokes.x=.1, strokes.y=.1, 
                      offsetting=TRUE, offset.labels=.0, offset.e= 1, 
                      axis.ext=.1, mai=c(.2, 1.5, .2, 1.5),
                      rect.margins=c(.01, .01),
                      srt=45,
                      cex.pos=.7,
                      xpd=TRUE,
                      c.lines=TRUE,  ### new
                      col.c.lines=grey(.9),
                      zoom=1,
                       ...)
{
  y <- showpoint <- showlabel <- type <- NULL       # to prevent 'R CMD check' from noting a missing binding 
                                                    # as the variables are provided in object x as default

  x <- x@plotdata          # df = data frame containing the information for printing
  
	max.all <- max(abs(x$x), abs(x$y))
	axis.ext <- 1 + axis.ext
  max.ext <- max.all * axis.ext
	
	x$x <- x$x * zoom     # zoom data
  x$y <- x$y * zoom     # zoom data
  
  # if (! draw.c) 
  #   x$labels[x$type %in% c("cl", "cr")] <- " "
	
	labels.constructs <- x$labels[x$type %in% c("cl", "cr")]
	labels.all <- x$labels
	
	if (flipaxes[1])
	  x$x <- x$x * -1
	if (flipaxes[2])
	  x$y <- x$y * -1
	
	# build plot
	old.par <- par(no.readonly = TRUE)      # save parameters
  #on.exit(par(old.par))                  # reset old par when done
	
	par(mai=mai)
	plot.new()
	plot.window(xlim = c(-max.ext, max.ext), 
	            ylim = c(-max.ext, max.ext), 
	            xaxs="i", yaxs="i", asp=1)
	
	# add center lines and outer rectangle
	segments(-max.ext, 0, max.ext, 0, col="lightgrey")
	segments(0, -max.ext, 0, max.ext, col="lightgrey")
	rect(-max.ext, -max.ext, max.ext, max.ext)
	
	# make standard concentration ellipse # TODO, scaling of ellipse
	#sing <- diag(esa$sing)[dim] / sum(diag(esa$sing))
	#ellipse(sing[1], sing[2], col="lightgrey")

    
  # initial coords for labels for strokes
  str.3 <- calcCoordsBorders(x["x"], x["y"], 
                             xmax=max.ext * (1 + strokes.x + offset.labels), # + rect.margins[1]/2), 
                             ymax=max.ext * (1 + strokes.y + offset.labels))# + rect.margins[2]/2))
  colnames(str.3) <- c("str.3.x", "str.3.y")                           
  x <- cbind(x, str.3)
  
  #segments(0,0,x$str.3.x, x$str.3.y)   # debug
  # calc coordinates for sourrounding rectangles (elements and constructs)
  lb <- calcRectanglesCoordsForLabels(x[, c("str.3.x", "str.3.y")], x$label, 
                                      cex=x$label.cex, x.ext=rect.margins[1], y.ext=rect.margins[2])
  colnames(lb) <- c("str.3.x0", "str.3.y0", "str.3.x1", "str.3.y1")
  x <- cbind(x, lb)
  #segments(x$str.3.x0, x$str.3.y0, x$str.3.x1, x$str.3.y1)   # debug
  
    
  # offset labels in y direction if too close together
  # for labels on the left and on the right seperately
  x$angle <- atan2(x$y, x$x)      # caveat: first y, then y argument!
  x <- x[order(x$angle), ]        # sort by angles

  # assign quandrants for offsetting
  x$quadrant[x$angle >= 0    & x$angle < pi/2] <- "ur"    # 
  x$quadrant[x$angle >= pi/2 & x$angle <= pi] <- "ul"     # 
  x$quadrant[x$angle < 0     & x$angle >= -pi/2] <- "lr"  # 
  x$quadrant[x$angle < -pi/2 & x$angle >= -pi] <- "ll"    # 
  
  
  # calc necessary offset (only correct in case there is overlap!)
  necessaryOffset <- function(a, b, direction=1, margin=.05){
    if (direction >= 0){        # offset upwards
      offset <- a[4] - b[2]     # is always positive >= 0
      if (offset < 0)           # if smaller than zero there should be no overlap anyway
        offset <- 0
    } else {                    # offset downwards
      offset <- a[2] - b[4]     # should always be <= 0
      if (offset > 0)           # if bigger than zero there is no overlap
        offset <- 0
    }
    as.numeric(offset + margin * sign(direction))
  }
  
  # offset quadrants
  #lr <- subset(x.sorted, type %in% c("cr", "cl") & quadrant=="lr")
  #ol <- lr[, c("str.3.x0", "str.3.y0", "str.3.x1", "str.3.y1")]

  #order.lines <- 1:nrow(ol)
  #order.lines <- rev(order.lines)
  
  # lim <- c(min(ol), max(ol))
  # plot(0,  type="n", xlim=lim, ylim=lim)
  # rect(ol[,1], ol[,2], ol[,3], ol[,4])
  # text(ol[,1], ol[,2], order.lines)
     
  offsetQuadrant <- function(x, quadrant="ur", direction=1, 
                             reverse=T, margin=0.02){
    index <- x$type %in% c("cr", "cl") & x$quadrant==quadrant   # get constructs of quandrant
    
    ol <- x[index, ]
    vars <- c("str.3.x0", "str.3.y0", "str.3.x1", "str.3.y1")
    
    order.lines <- 1:nrow(ol)
    if (reverse)
      order.lines <- rev(order.lines)
   
    for (i in order.lines){
      for (i.n in order.lines){
        if(i != i.n){
          overlap <- doRectanglesOverlap(ol[i, vars], ol[i.n, vars])
          if (overlap){    # if overlap is present the rectangles is moved to avoid overlap
            offset <- necessaryOffset(ol[i, vars], ol[i.n, vars], dir=direction, margin=margin)
            ol[i.n, c("str.3.y", "str.3.y0","str.3.y1")] <- 
                ol[i.n, c("str.3.y", "str.3.y0", "str.3.y1")] + offset
          }
        }
      }
    }
    
    x[index, c("str.3.y", vars)] <- ol[, c("str.3.y", vars)]
    x
  }
  
  # code is slow!
  if (outer.positioning){
    x <- offsetQuadrant(x, quadrant="ur", direction=1, reverse=F)    #dir.ur <-  1; reverse <- F 
    x <- offsetQuadrant(x, quadrant="ul", direction=1, reverse=T)    #dir.ul <-  1; reverse <- T
    x <- offsetQuadrant(x, quadrant="ll", direction=-1, reverse=F)   #dir.ll <- -1; reverse <- F 
    x <- offsetQuadrant(x, quadrant="lr", direction=-1, reverse=T)   #dir.lr <- -1; reverse <- T
  }
  # 
  # for (i in order.lines){
  #   #cat("---\n")
  #   for (i.n in order.lines){
  #     if(i != i.n){
  #       overlap <- doRectanglesOverlap(ol[i, ], ol[i.n, ])
  #       #cat("(", i, i.n, ")", "\t\t"); print(overlap)
  #       if (overlap){    # if overlap is present the rectangles is moved to avoid overlap  
  #         offset <- necessaryOffset(ol[i, ], ol[i.n, ], dir=-1, margin=0.02)
  #         #print(offset)
  #         ol[i.n, c("str.3.y0","str.3.y1")] <- 
  #             ol[i.n, c("str.3.y0","str.3.y1")] + offset
  #       }
  #     }
  #   }
  # }

  # lim <- c(min(ol), max(ol))
  #   rect(ol[,1], ol[,2], ol[,3], ol[,4], border="blue", lty=2)
  #   text(ol[,3], ol[,2], order.lines, col="blue" )
  # 
  # 
  #   plot(0:5)
  #   rect(ol[4,1],ol[4,2],ol[4,3],ol[4,4])
  #   rect(ol[3,1],ol[3,2],ol[3,3],ol[3,4])
  #   doRectanglesOverlap(ol[4,], ol[3,])
  # do others overlap? If yes move them
  
  # make outer strokes for all labels (elements and constructs) 
  # select which to draw later
  # coordinates for stroke starts
  str.1 <- calcCoordsBorders(x["x"], x["y"], xmax=max.ext, ymax=max.ext)
  colnames(str.1) <- c("str.1.x", "str.1.y")
  
  # coordinates for stroke ends
  str.2 <- calcCoordsBorders(x["x"], x["y"], xmax=max.ext * (1 + strokes.x), 
                             ymax=max.ext * (1 + strokes.y))
  colnames(str.2) <- c("str.2.x", "str.2.y")
  
  x <- cbind(x, str.1, str.2)
  
  # redo coordinates for stroke ends according to edges of rectangles that have been offsetted
  a <- list()
  for (i in seq_len(nrow(x))){
    a[[i]] <- calcCoordsBorders(x[i, "x"], x[i, "y"], xmax=max.ext * (1 + strokes.x), 
                               ymax=abs(x[i, "str.3.y"]))
  }
  str.4 <- do.call(rbind, a)
  colnames(str.4) <- c("str.4.x", "str.4.y")
  x <- cbind(x, str.4)

  if (!c.labels.inside){   # when constructs labels are  prompted to be outside the plot(default) 
    # rotate labels srt degress on top and bottom for quick printing
    y.max.ext <- max.ext * (1 + strokes.y + offset.labels)  # position of outer strokes to determine side of labels
    x$rotate <- 0                         # default is no rotation of labels in text 
    if (!outer.positioning)               # only for positioning = FALSE to get neater label directions 
      x$rotate[x$str.3.y == y.max.ext |   # replace by standadrd rotation angle
               x$str.3.y == -y.max.ext] <- srt
    
    # only make labels, rectangles and strokes that are prompted
    cl <- subset(x, type %in% c("cl", "cr") & showlabel==T)    # select only labels that should be shown
    segments(cl$str.1.x, cl$str.1.y, cl$str.2.x, cl$str.2.y, xpd=T)
    segments(cl$str.2.x, cl$str.2.y, cl$str.4.x, cl$str.4.y, xpd=T, lty=3)
    rect(cl$str.3.x0, cl$str.3.y0, 
         cl$str.3.x1, cl$str.3.y1, col=grey(1), border=grey(1), xpd=T)

    # print constructs labels (if there are any) and not only inner labels are prompted
    if (nrow(cl) > 0){
      for (i in 1:nrow(cl)){
        if (cl$str.3.x[i] < 0) 
          adj <- c(1, .5) else
          adj <- c(0, .5)
        if (!outer.positioning){      # overwrite adj in case of no positioning
          if (cl$str.3.y[i] == y.max.ext)  
            adj <- c(0, .5)
          if (cl$str.3.y[i] == -y.max.ext)  
            adj <- c(1, .5)
        }
        text(cl$str.3.x[i], cl$str.3.y[i], labels=cl$label[i], col=cl$label.col[i], 
             cex=cl$label.cex[i], adj=adj, xpd=T, srt=cl$rotate[i])   
      }
    }
  }
  
  
  ### plotting of elements and contructs inside plot ###
  
  #make construct lines if prompted 
  if (c.lines){ 
    cli <- subset(x, type %in% c("cl", "cr") & showlabel==T)    # select only labels that should be shown
    segments(0, 0, cli$str.1.x, cli$str.1.y, col=col.c.lines)   # lines form biplot center to outsides
  }

  # make construct symbols
  cs <- subset(x, type %in% c("cl", "cr") & showpoint==T & abs(x)<max.ext & abs(y)<max.ext)
	points(cs[c("x", "y")], col=cs$point.col, pch=4, cex=cs$point.cex, xpd=xpd)  

  # make element symbols
  es <- subset(x, type=="e" & showpoint==T & abs(x)<max.ext & abs(y)<max.ext)
	points(es[c("x", "y")], col=es$point.col, pch=15, cex=es$point.cex, xpd=xpd)

	# positioning of element and constructs labels inside the plot
	if (inner.positioning) {    # positioning by landmark algorithm from maptools package
		
    # dirty hack as I do not understand the problem why some values become NA: 
    #  replace NAs in showlabel and showpoint. Endas grid #10 produces this
		x$showlabel[is.na(x$showlabel)] <- TRUE 
		x$showpoint[is.na(x$showpoint)] <- TRUE 
		
    sh <- subset(x, showlabel==T & showpoint==T)# & 
		lpos <- pointLabel(sh[c("x", "y")], labels=sh$label, doPlot=FALSE, cex=cex.pos)     # package maptools
	  x$x.pos <- NA
	  x$y.pos <- NA
	  sh$x.pos <- lpos$x
	  sh$y.pos <- lpos$y	  
	  x[x$showlabel==T & x$showpoint==T, ] <- sh
	} else {              # simple offsetting in y direction
	  x$x.pos <- x$x
	  x$y.pos <- NA
	  offset.y.pos <- strheight("aaaa", cex=.7)  # string height for dummy string
	  x[x$type=="e", ]$y.pos <- x[x$type=="e", ]$y + offset.y.pos * offset.e  # offset element labels by normal stringheight times x
	  x[x$type %in% c("cl", "cr"), ]$y.pos <- x[x$type %in% c("cl", "cr"), ]$y - .05
	}
	
	# text labels for elements
	#es <- subset(x, type=="e" & showlabel==T & showpoint==T)   # old version
	es <- subset(x, type=="e" & showlabel==T & showpoint==T & 
	             abs(x)<max.ext & abs(y)<max.ext)             # avoid plotting outside plot region
  if (dim(es)[1] > 0)
    text(es[, c("x.pos", "y.pos")], 
         labels=es$label, col=es$label.col, pch=15, cex=es$label.cex, xpd=xpd)
  
  # text labels for constructs inside plot
  if (c.labels.inside){
    cs <- subset(x, type %in% c("cl", "cr") & showlabel==T 
                 & abs(x)<max.ext & abs(y)<max.ext)
  	if (dim(cs)[1] > 0)
  	text(cs[, c("x.pos", "y.pos")], 
          labels=cs$label, col=cs$label.col, pch=4, cex=cs$label.cex, xpd=xpd)
  }
  invisible(x)    # returns plotdata frame
}
# x <- calcBiplotCoords(raeithel, g=1, h=1)
# x <- prepareBiplotData(x)
# biplotDraw(x))    # add amount explained variance to the axes




#' Adds the percentage of the sum-of-squares explained by each axis to the plot.
#' 
#' @param x               \code{repgrid} object containing the biplot coords, i.e. after 
#'                        having called \code{\link{calcBiplotCoords}} and 
#'                        \code{\link{prepareBiplotData}}.
#' @param dim             The dimensions to be printed.
#' @param var.show        Show explained sum-of-squares in biplot? (default \code{TRUE}). 
#' @param var.cex         The cex value for the percentages shown in the plot.
#' @param var.col         The color value of the percentages shown in the plot.
#' @param axis.ext        Axis extension factor (default is \code{.1}). A bigger value will 
#'                        zoom out the plot.
#' @param  center		      Numeric. The type of centering to be performed. 
#'                        0= no centering, 1= row mean centering (construct), 
#'                        2= column mean centering (elements), 3= double-centering (construct and element means),
#'                        4= midpoint centering of rows (constructs).
#'                        The default is \code{1} (row centering).
#' @param normalize       A numeric value indicating along what direction (rows, columns)
#'                        to normalize by standard deviations. \code{0 = none, 1= rows, 2 = columns}
#'                        (default is \code{0}).
#' @param g               Power of the singular value matrix assigned to the left singular 
#'                        vectors, i.e. the constructs.
#' @param h               Power of the singular value matrix assigned to the right singular 
#'                        vectors, i.e. the elements.
#' @param col.active      Columns (elements) that are no supplementary points, i.e. they are used
#'                        in the SVD to find principal components. default is to use all elements.
#' @param col.passive     Columns (elements) that are supplementary points, i.e. they are NOT used
#'                        in the SVD but projecte into the component space afterwards. They do not 
#'                        determine the solution. Default is \code{NA}, i.e. no elements are set 
#'                        supplementary.
#' @param ...             Not evaluated.
#'
#' @author   Mark Heckmann
#' @keywords internal
#' @export
#'
addVarianceExplainedToBiplot2d <- function(x, dim=c(1,2,3), var.cex=.7, 
                                          var.show=TRUE, var.col=grey(.1), 
                                          axis.ext = .1, 
                                          center=1, normalize=0,
                                          g=0, h=1-g, 
                                          col.active=NA, 
                                          col.passive=NA, 
                                          ...){  
  # do only if prompted
  if (var.show){
    # determine way to calculate SSQ proportions. 
    # Different if passive columns are used.
    if(is.na(col.active[1]) & is.na(col.passive[1]))
      standard.calc.ssq <- TRUE else 
      standard.calc.ssq <- FALSE
  
    if (standard.calc.ssq){ 
      # one valid way of calculating the prop SSQ not taking into account passive columns
      sv <- x@calcs$biplot$D        # get singular values from SVD
      sv.exp <- sv^2/sum(sv^2)      # proportion of ssq explained per principal component
      var <- paste("Dim ", dim[1:2], ": ", round(sv.exp[dim[1:2]] * 100, 1), "%", sep="")
    } else {
      # calculating explained variance when passive columns are used
      ssq.out <- ssq(x, along=2, cum=F, g=g, h=h, 
                     center=center, normalize=normalize, 
                     col.active=col.active, 
                     col.passive=col.passive, print=F, ...)
      ssq.prop.dim <- ssq.out[dim(ssq.out)[1], dim]
      var <- paste("Dim ", dim[1:2], ": ", 
                   round(ssq.prop.dim[dim[1:2]], 1), "%", sep="")
    }             
    data <- x@plotdata          # data frame from data prepare function return
    max.all <- max(abs(data$x), abs(data$y))
    axis.ext <- 1 + axis.ext
    max.ext <- max.all * axis.ext
  
    ext <- strheight(var[1], cex=var.cex)
  
    text(max.ext - ext/2, 0, var[1] , cex=var.cex, adj=c(.5,0), col=var.col, srt=90)
    text(0, -max.ext + ext/2, var[2] , cex=var.cex, adj=c(.5,0), col=var.col)
  }
}

# x <- randomGrid(20, 40)
# x <- boeker
# x <- raeithel
# xb <- prepareBiplotData(x, c.labels.show=F, c.points.dev=90, e.points=1:3, e.labels=T)
# biplotDraw(xb, xpd=F, inner=F, outer=F)
# addVarianceExplainedToBiplot(xb$rg, xb$df)
# 
# xb <- prepareBiplotData(x, c.points.dev=5, c.labels.dev=5)
# biplotDraw(xb, xpd=F, inner=F, outer=T, mai=rep(0,4), c.labels.inside=T)
# biplotDraw(xb, xpd=F, inner=F, outer=T)
# 
# dev.new()
# xb <- prepareBiplotData(x, dim=c(1,2), map=4)
# biplotDraw(xb, dev=15)
# addVarianceExplainedToBiplot(x, xb, dim=c(1,2,4))
# 
# 
# dev.new()
# xb <- prepareBiplotData(x, dim=c(2,3), map=1)
# biplotDraw(xb, dev=15)
# addVarianceExplainedToBiplot(x, xb, dim=c(2,3,1))
# 
# dev.new()
# xb <- prepareBiplotData(x, dim=c(3,4), map=1)
# biplotDraw(xb, dev=15)
# addVarianceExplainedToBiplot(x, xb, dim=c(3,4,1))

# x <- boeker
# x <- prepareBiplotData(x, e.col="black", c.col="black", cex.e=.7, cex.c=.7)#, color.e=.8, color.c=.8)
# biplotDraw(x)
# addVarianceExplainedToBiplot(boeker, x)
# 
# x <- boeker
# x <- prepareBiplotData(x, e.col="black", c.col="black", cex.e=.7, cex.c=.7, 
#                         color.e=.3, color.c=.5)
# biplotDraw(x)
# addVarianceExplainedToBiplot(boeker, x)
# 
# x <- boeker
# x <- prepareBiplotData(x, e.col="black", c.col="black", cex.e=c(.3,1))
# x <- prepareBiplotData(x, e.col="black", c.col="black",  cex.e=c(.3,1), cex.c=c(.3,1))
# x <- prepareBiplotData(x, cex.e=c(.5,1.3), cex.c=c(.5,1.3))
# x <- prepareBiplotData(x, cex.e=c(.5,1), cex.c=c(.5,1), color.c.map=c(0, 0))
# biplotDraw2(x)
# 
# x <- boeker
# x <- raeithel
# 
# layout(matrix(1:4, by=T, ncol=2))
# 
# xb <- prepareBiplotData(x, dim=c(1,2), map=3)
# biplotDraw(xb)
# addVarianceExplainedToBiplot(x, xb, dim=1:3)
# 
# xb <- prepareBiplotData(x, dim=c(2,3), map=1)
# biplotDraw(xb)
# addVarianceExplainedToBiplot(x, xb, dim=c(2,3,1))
# 
# xb <- prepareBiplotData(x, dim=c(3,4), map=1)
# biplotDraw(xb)
# addVarianceExplainedToBiplot(x, xb, dim=c(3,4,1))
# 
# xb <- prepareBiplotData(x, dim=c(1,4), map=2)
# biplotDraw(xb)
# addVarianceExplainedToBiplot(x, xb, dim=c(1,4,2))




#' Draw a two-dimensional biplot. 
#'
#' The biplot is the central
#' way to create a joint plot of elements and constructs.
#' Depending on te parameters chosen it contains information
#' on the distances between elements and constructs. Also the 
#' relative values the elements have on a construct can be read off
#' by projetion the element onto the construct vector. 
#' A lot of parameters can be changed rendering
#' different types of biplots (ESA, Slater's) and different 
#' looks (colors, text size).
#' See the example section below to get started.
#'
#' For the construction of a biplot the grid matrix is first
#' centered and normalized according to the prompted options.\cr 
#' Next, the matrix is decomposed by singular value decomposition (SVD)
#' into \deqn{X = UDV^T}{X = UDV^T}
#' The biplot is made up of two matrices 
#' \deqn{X = GH^T}{X = GH^T}
#' These matrices are construed on the basis of the SVD results.
#' \deqn{\hat{X} = UD^gD^hV^T}{X = UD^gD^hV^T}
#' Note that the grid matrix values are only recovered and 
#' the projection property is only given if \eqn{g + h = 1}{g + h = 1}
#' 
#'
#' @param x                   \code{repgrid} object.
#' @param dim                 Dimensions (i.e. principal components) to be used for biplot 
#'                            (default is \code{c(1,2)}).
#' @param map.dim             Third dimension (depth) used to map aesthetic attributes to
#'                            (default is \code{3}).
#' @param  center		          Numeric. The type of centering to be performed. 
#'                            0= no centering, 1= row mean centering (construct), 
#'                            2= column mean centering (elements), 3= double-centering (construct and element means),
#'                            4= midpoint centering of rows (constructs).
#'                            The default is \code{1} (row centering).
#' @param normalize           A numeric value indicating along what direction (rows, columns)
#'                            to normalize by standard deviations. \code{0 = none, 1= rows, 2 = columns}
#'                            (default is \code{0}).
#' @param g                   Power of the singular value matrix assigned to the left singular 
#'                            vectors, i.e. the constructs.
#' @param h                   Power of the singular value matrix assigned to the right singular 
#'                            vectors, i.e. the elements.
#' @param col.active          Columns (elements) that are no supplementary points, i.e. they are used
#'                            in the SVD to find principal components. default is to use all elements.
#' @param col.passive         Columns (elements) that are supplementary points, i.e. they are NOT used
#'                            in the SVD but projecte into the component space afterwards. They do not 
#'                            determine the solution. Default is \code{NA}, i.e. no elements are set 
#'                            supplementary.
#' @param e.point.col         Color of the element symbols. The default is \code{"black"}.
#'                            Two values can be entered that will create a color ramp. The values of 
#'                            \code{map.dim} are mapped onto the ramp.
#'                            If only one color color value is supplied (e.g. \code{"black"}) 
#'                            no mapping occurs and all elements will have the same color 
#'                            irrespective of their value on the \code{map.dim} dimension.
#' @param e.point.cex         Size of the element symbols. The default is \code{.9}.
#'                            Two values can be entered that will create a size ramp. The values of 
#'                            \code{map.dim} are mapped onto the ramp.
#'                            If only one color size value is supplied (e.g. \code{.8}) 
#'                            no mapping occurs and all elements will have the same size 
#'                            irrespective of their value on the \code{map.dim} dimension.
#' @param e.label.col         Color of the element label. The default is \code{"black"}.
#'                            Two values can be entered that will create a color ramp. The values of 
#'                            \code{map.dim} are mapped onto the ramp.
#'                            If only one color color value is supplied (e.g. \code{"black"}) 
#'                            no mapping occurs and all labels will have the same color 
#'                            irrespective of their value on the \code{map.dim} dimension.
#' @param e.label.cex         Size of the element labels. The default is \code{.7}.
#'                            Two values can be entered that will create a size ramp. The values of 
#'                            \code{map.dim} are mapped onto the ramp.
#'                            If only one color size value is supplied (e.g. \code{.7}) 
#'                            no mapping occurs and all labels will have the same size 
#'                            irrespective of their value on the \code{map.dim} dimension.
#' @param e.color.map         Value range to determine what range of the color ramp defined in 
#'                            \code{e.color} will be used for mapping the colors. 
#'                            Default is \code{c(.4, ,1)}. Usually not important for the user. 
#' @param c.point.col         Color of the construct symbols. The default is \code{"black"}.
#'                            Two values can be entered that will create a color ramp. The values of 
#'                            \code{map.dim} are mapped onto the ramp.
#'                            If only one color color value is supplied (e.g. \code{"black"}) 
#'                            no mapping occurs and all construct will have the same color 
#'                            irrespective of their value on the \code{map.dim} dimension.
#' @param c.point.cex         Size of the construct symbols. The default is \code{.8}.
#'                            Two values can be entered that will create a size ramp. The values of 
#'                            \code{map.dim} are mapped onto the ramp.
#'                            If only one color size value is supplied (e.g. \code{.8}) 
#'                            no mapping occurs and all construct will have the same size 
#'                            irrespective of their value on the \code{map.dim} dimension.
#' @param c.label.col         Color of the construct label. The default is \code{"black"}.
#'                            Two values can be entered that will create a color ramp. The values of 
#'                            \code{map.dim} are mapped onto the ramp.
#'                            If only one color color value is supplied (e.g. \code{"black"}) 
#'                            no mapping occurs and all labels will have the same color 
#'                            irrespective of their value on the \code{map.dim} dimension.
#' @param c.label.cex         Size of the construct labels. The default is \code{.7}.
#'                            Two values can be entered that will create a size ramp. The values of 
#'                            \code{map.dim} are mapped onto the ramp.
#'                            If only one color size value is supplied (e.g. \code{.7}) 
#'                            no mapping occurs and all labels will have the same size 
#'                            irrespective of their value on the \code{map.dim} dimension.
#' @param c.color.map         Value range to determine what range of the color ramp defined in 
#'                            \code{c.color} will be used for mapping. Default is \code{c(.4, ,1)}.
#'                            Usually not important for the user.
#' @param c.points.devangle   The deviation angle from the x-y plane in degrees. These can only be calculated
#'                            if a third dimension \code{map.dim} is specified. Only the constructs 
#'                            that do not depart more than the specified degrees from the 
#'                            x-y plane will be printed. This facilitates the visual 
#'                            interpretation, as only vectors represented near the current plane 
#'                            are shown. Set the value to \code{91} (default) 
#'                            to show all vectors.
#' @param c.labels.devangle   The deviation angle from the x-y plane in degrees. These can only be calculated
#'                            if a third dimension \code{map.dim} is specified. Only the labels of constructs 
#'                            that do not depart more than the specified degrees from the 
#'                            x-y plane will be printed. Set the value to \code{91} (default) 
#'                            to show all construct labels.
#' @param c.points.show       Whether the constructs are printed (default is \code{TRUE}).
#'                            \code{FALSE} will surpress the printing of the constructs.
#'                            To only print certain constructs a numeric vector can be 
#'                            provided (e.g. \code{c(1:10)}).
#' @param c.labels.show       Whether the construct labels are printed (default is \code{TRUE}).
#'                            \code{FALSE} will surpress the printing of the labels.
#'                            To only print certain construct labels a numeric vector can be 
#'                            provided (e.g. \code{c(1:10)}).
#' @param e.points.show       Whether the elements are printed (default is \code{TRUE}).
#'                            \code{FALSE} will surpress the printing of the elements.
#'                            To only print certain elements a numeric vector can be 
#'                            provided (e.g. \code{c(1:10)}).
#' @param e.labels.show       Whether the element labels are printed (default is \code{TRUE}).
#'                            \code{FALSE} will surpress the printing of the labels.
#'                            To only print certain element labels a numeric vector can be 
#'                            provided (e.g. \code{c(1:10)}).
#' @param inner.positioning   Logical. Whether to calculate positions to minimize overplotting of 
#'                            elements and construct labels (default is\code{TRUE}). Note that
#'                            the positioning may slow down the plotting.
#' @param outer.positioning   Logical. Whether to calculate positions to minimize overplotting of 
#'                            of construct labels on the outer borders (default is\code{TRUE}). Note that
#'                            the positioning may slow down the plotting.
#' @param c.labels.inside     Logical. Whether to print construct labels next to the points.
#'                            Can be useful during inspection of the plot (default \code{FALSE}).
#' @param c.lines             Logical. Whether construct lines from the center of the biplot
#'                            to the sourrounding box are drawn (default is \code{FALSE}).
#' @param col.c.lines         The color of the construct lines from the center to the borders 
#'                            of the plot (default is \code{gray(.9)}).
#' @param flipaxes            Logical vector of length two. Whether x and y axes are reversed 
#'                            (default is \code{c(F,F)}).
#' @param strokes.x           Length of outer strokes in x direction in NDC.  
#' @param strokes.y           Length of outer strokes in y direction in NDC.
#' @param offsetting          Do offsetting? (TODO)
#' @param offset.labels       Offsetting parameter for labels (TODO).
#' @param offset.e            offsetting parameter for elements (TODO).
#' @param axis.ext            Axis extension factor (default is \code{.1}). A bigger value will 
#'                            zoom out the plot.
#' @param mai                 Margins available for plotting the labels in inch 
#'                            (default is \code{c(.2, 1.5, .2, 1.5)}).
#' @param rect.margins        Vector of length two (default is \code{c(.07, .07)}). Two values
#'                            specifying the additional horizontal and vertical margin around each 
#'                            label.      
#' @param srt                 Angle to rotate construct label text. Only used in case \code{offsetting=FALSE}.
#' @param cex.pos             Cex parameter used during positioning of labels if prompted. Does
#'                            usually not have to be changed by user.
#' @param xpd                 Logical (default is \code{TRUE}). Wether to extend text labels 
#'                            over figure region. Usually not needed by the user.
#' @param unity               Scale elements and constructs coordinates to unit scale in 2D (maximum of 1)
#'                            so they are printed more neatly (default \code{TRUE}).
#' @param unity3d             Scale elements and constructs coordinates to unit scale in 3D (maximum of 1)
#'                            so they are printed more neatly (default \code{TRUE}).
#' @param scale.e             Scaling factor for element vectors. Will cause element points to move a bit more
#'                            to the center. (but only if \code{unity} or \code{unity3d} is \code{TRUE}).
#'                            This argument is for visual appeal only.
#' @param zoom                Scaling factor for all vectors. Can be used to zoom
#'                            the plot in and out (default \code{1}).
#' @param var.show            Show explained sum-of-squares in biplot? (default \code{TRUE}). 
#' @param var.cex             The cex value for the percentages shown in the plot.
#' @param var.col             The color value of the percentages shown in the plot.
#' @param ...                 parameters passed on to  come.
#'
#' @author                    Mark Heckmann
#' @export
#'
#' @seealso   Unsophisticated biplot: \code{\link{biplotSimple}}; \cr
#'            2D biplots:
#'            \code{\link{biplot2d}},
#'            \code{\link{biplotEsa2d}},
#'            \code{\link{biplotSlater2d}};\cr
#'            Pseudo 3D biplots:
#'            \code{\link{biplotPseudo3d}},  
#'            \code{\link{biplotEsaPseudo3d}},
#'            \code{\link{biplotSlaterPseudo3d}};\cr
#'            Interactive 3D biplots:
#'            \code{\link{biplot3d}},
#'            \code{\link{biplotEsa3d}},
#'            \code{\link{biplotSlater3d}};\cr
#'            Function to set view in 3D:
#'            \code{\link{home}}.  
#'
#' @examples \dontrun{
#'
#'    biplot2d(boeker)                # biplot of boeker data
#'    biplot2d(boeker, c.lines=T)     # add construct lines
#'    biplot2d(boeker, center=2)      # with column centering
#'    biplot2d(boeker, center=4)      # midpoint centering
#'    biplot2d(boeker, normalize=1)   # normalization of constructs
#'
#'    biplot2d(boeker, dim=2:3)       # plot 2nd and 3rd dimension
#'    biplot2d(boeker, dim=c(1,4))    # plot 1st and 4th dimension
#'
#'    biplot2d(boeker, g=1, h=1)            # assign singular values to con. & elem.
#'    biplot2d(boeker, g=1, h=1, center=1)  # row centering (Slater)
#'    biplot2d(boeker, g=1, h=1, center=4)  # midpoint centering (ESA)
#'
#'    biplot2d(boeker, e.color="red", c.color="blue")   # change colors
#'    biplot2d(boeker, c.color=c("white", "darkred"))   # mapped onto color range
#'    
#'    biplot2d(boeker, unity=T)                 # scale con. & elem. to equal length
#'    biplot2d(boeker, unity=T, scale.e=.5)     # scaling factor for element vectors
#'
#'    biplot2d(boeker, e.labels.show=F)         # do not show element labels
#'    biplot2d(boeker, e.labels.show=c(1,2,4))  # show labels for elements 1, 2 and 4
#'    biplot2d(boeker, e.points.show=c(1,2,4))  # only show elements 1, 2 and 4
#'    biplot2d(boeker, c.labels.show=c(1:4))    # show constructs labels 1 to 4   
#'    biplot2d(boeker, c.labels.show=c(1:4))    # show constructs labels except 1 to 4
#'
#'    biplot2d(boeker, e.cex.map=1)   # change size of texts for elements 
#'    biplot2d(boeker, c.cex.map=1)   # change size of texts for constructs 
#'
#'    biplot2d(boeker, g=1, h=1, c.labels.inside=T)  # constructs inside the plot
#'    biplot2d(boeker, g=1, h=1, c.labels.inside=T,  # different margins and elem. color 
#'             mai=c(0,0,0,0), e.color="red") 
#'  
#'    biplot2d(boeker, strokes.x=.3, strokes.y=.05)  # change length of strokes
#'
#'    biplot2d(boeker, flipaxes=c(T, F))      # flip x axis
#'    biplot2d(boeker, flipaxes=c(T, T))      # flip x and y axis
#'
#'    biplot2d(boeker, outer.positioning=F)   # no positioning of con.-labels
#'
#'    biplot2d(boeker, c.labels.devangle=20)  # only con. within 20 degree angle
#' } 
#'
biplot2d <- function(x, dim=c(1,2), map.dim=3, 
                    center=1,
                    normalize=0, 
                    g=0, 
                    h=1-g, 
                    col.active=NA, 
                    col.passive=NA,
                    #e.color="black", 
                    #c.color="black",   
                    e.point.col="black",
                    e.point.cex=.9,
                    e.label.col="black",
                    e.label.cex=.7,   
                    e.color.map=c(.4, 1),               
                    c.point.col="black",
                    c.point.cex=.8,
                    c.label.col="black",
                    c.label.cex=.7,
                    c.color.map=c(.4, 1),
                    #e.cex.map=.7,
                    #c.cex.map=.7,
                    c.points.devangle=91,
                    c.labels.devangle=91,
                    c.points.show=TRUE,
                    c.labels.show=TRUE,
                    e.points.show=TRUE,
                    e.labels.show=TRUE,
                    inner.positioning=TRUE,
                    outer.positioning=TRUE,
                    c.labels.inside=FALSE,
                    c.lines=TRUE, 
                    col.c.lines=grey(.9),
                    flipaxes=c(FALSE,FALSE), 
                    strokes.x=.1, strokes.y=.1, 
                    offsetting=TRUE, offset.labels=.0, offset.e= 1, 
                    axis.ext=.1, mai=c(.2, 1.5, .2, 1.5),
                    rect.margins=c(.01, .01),
                    srt=45,
                    cex.pos=.7,
                    xpd=TRUE, 
                    unity=FALSE, 
                    unity3d=FALSE,
                    scale.e=.9,
                    zoom=1,
                    var.show=TRUE, 
                    var.cex=.7, 
                    var.col=grey(.1),  
                    ...)
{
  x <- calcBiplotCoords(x, center=center, normalize=normalize, 
                        g=g, h=h, 
                        col.active=col.active, col.passive=col.passive, ...)
  x <- prepareBiplotData(x, dim=dim, map.dim=map.dim, 
                          e.label.cex=e.label.cex,  c.label.cex=c.label.cex,
                          e.label.col=e.label.col,  c.label.col=c.label.col,
                          e.point.cex=e.point.cex,  c.point.cex=c.point.cex,
                          e.point.col=e.point.col,  c.point.col=c.point.col,
                          #e.color=e.color, c.color=c.color,
                          #e.cex.map=e.cex.map, c.cex.map=c.cex.map, 
                          e.color.map=e.color.map, c.color.map=c.color.map, 
                          c.points.devangle=c.points.devangle,
                          c.labels.devangle=c.labels.devangle, c.points.show=c.points.show,
                          c.labels.show=c.labels.show, 
                          e.points.show=e.points.show,
                          e.labels.show=e.labels.show, 
                          unity=unity, unity3d=unity3d, scale.e=scale.e, ...)
  biplotDraw(x, inner.positioning=inner.positioning, outer.positioning=outer.positioning,
             c.labels.inside=c.labels.inside, 
             c.lines=c.lines, col.c.lines=col.c.lines, flipaxes=flipaxes, 
             strokes.x=strokes.x, strokes.y=strokes.y, 
             offsetting=offsetting, offset.labels=offset.labels, offset.e=offset.e, 
             axis.ext=axis.ext, mai=mai, rect.margins=rect.margins,
             srt=srt, cex.pos=cex.pos, xpd=xpd, zoom=zoom)
  addVarianceExplainedToBiplot2d(x, dim=dim, center=center, normalize=normalize, 
                                 g=g, h=h, col.active=col.active, 
                                 col.passive=col.passive, var.show=var.show,
                                 var.cex=var.cex, var.col=var.col, ...)
  invisible(x)                                           
}


#' See \code{\link{biplotPseudo3d}} for its use.

#' Draws a biplot of the grid in 2D with depth impression (pseudo 3D).
#'
#' This version is basically a 2D biplot. 
#' It only modifies color and size of the symbols in order to create a 3D impression
#' of the data points. 
#' This function will call the standard \code{\link{biplot2d}} function with some 
#' modified arguments. For the whole set of arguments that can be used
#' see \code{\link{biplot2d}}. Here only the arguments special to 
#' \code{biplotPseudo3d} are outlined.
#'
#' @param x             \code{repgrid} object.
#' @param dim           Dimensions (i.e. principal components) to be used for biplot 
#'                      (default is \code{c(1,2)}).
#' @param map.dim       Third dimension (depth) used to map aesthetic attributes to
#'                      (default is \code{3}).
#' @param e.point.col   Color(s) of the element symbols. Two values can be entered that will
#'                      create a color ramp. The values of \code{map.dim} are mapped onto the ramp.
#'                      The default is \code{c("white", "black")}. If only one color color value
#'                      is supplied (e.g. \code{"black"}) no mapping occurs and all elements
#'                      will have the same color irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param e.point.cex   Size of the element symbols. Two values can be entered that will
#'                      represents the lower and upper size of a range of cex the values of \code{map.dim} 
#'                      are mapped onto. The default is \code{c(.6, 1.2)}. If only one cex value
#'                      is supplied (e.g. \code{.7}) no mapping occurs and all elements
#'                      will have the same size irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param e.label.col   Color(s) of the element labels. Two values can be entered that will
#'                      create a color ramp. The values of \code{map.dim} are mapped onto the ramp.
#'                      The default is \code{c("white", "black")}. If only one color color value
#'                      is supplied (e.g. \code{"black"}) no mapping occurs and all element labels
#'                      will have the same color irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param e.label.cex   Size of the element labels. Two values can be entered that will
#'                      represents the lower and upper size of a range of cex the values of \code{map.dim} 
#'                      are mapped onto. The default is \code{c(.6, .8)}. If only one cex value
#'                      is supplied (e.g. \code{.7}) no mapping occurs and all element labels
#'                      will have the same size irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param e.color.map   Value range to determine what range of the color ramp defined in 
#'                      \code{e.color} will be used for mapping the colors. 
#'                      Default is \code{c(.4, ,1)}. Usually not important for the user. 
#' @param c.point.col   Color(s) of the construct symbols. Two values can be entered that will
#'                      create a color ramp. The values of \code{map.dim} are mapped onto the ramp.
#'                      The default is \code{c("white", "darkred")}. If only one color color value
#'                      is supplied (e.g. \code{"black"}) no mapping occurs and all elements
#'                      will have the same color irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param c.point.cex   Size of the construct symbols. Two values can be entered that will
#'                      represents the lower and upper size of a range of cex the values of \code{map.dim} 
#'                      are mapped onto. The default is \code{c(.6, 1.2)}. If only one cex value
#'                      is supplied (e.g. \code{.7}) no mapping occurs and all elements
#'                      will have the same size irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param c.label.col   Color(s) of the construct labels. Two values can be entered that will
#'                      create a color ramp. The values of \code{map.dim} are mapped onto the ramp.
#'                      The default is \code{c("white", "black")}. If only one color color value
#'                      is supplied (e.g. \code{"black"}) no mapping occurs and all construct labels
#'                      will have the same color irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param c.label.cex   Size of the construct labels. Two values can be entered that will
#'                      represents the lower and upper size of a range of cex the values of \code{map.dim} 
#'                      are mapped onto. The default is \code{c(.6, .9)}. If only one cex value
#'                      is supplied (e.g. \code{.7}) no mapping occurs and all construct labels
#'                      will have the same size irrespective of their value on the \code{map.dim}
#'                      dimension.
#' @param c.color.map   Value range to determine what range of the color ramp defined in 
#'                      \code{c.color} will be used for mapping. Default is \code{c(.4, ,1)}.
#'                      Usually not important for the user.
#' @param ...           Additional parameters passed to \code{\link{biplot2d}}.
#'  
#' @author              Mark Heckmann
#' @export
#'
#' @seealso   Unsophisticated biplot: \code{\link{biplotSimple}}; \cr
#'            2D biplots:
#'            \code{\link{biplot2d}},
#'            \code{\link{biplotEsa2d}},
#'            \code{\link{biplotSlater2d}};\cr
#'            Pseudo 3D biplots:
#'            \code{\link{biplotPseudo3d}},  
#'            \code{\link{biplotEsaPseudo3d}},
#'            \code{\link{biplotSlaterPseudo3d}};\cr
#'            Interactive 3D biplots:
#'            \code{\link{biplot3d}},
#'            \code{\link{biplotEsa3d}},
#'            \code{\link{biplotSlater3d}};\cr
#'            Function to set view in 3D:
#'            \code{\link{home}}.
#'
#' @examples \dontrun{
#'    # biplot with 3D impression
#'    biplotPseudo3d(boeker)    
#'    # Slater's biplot with 3D impression                  
#'    biplotPseudo3d(boeker, g=1, h=1, center=1)  
#'
#'    # show 2nd and 3rd dim. and map 4th 
#'    biplotPseudo3d(boeker, dim=2:3, map.dim=4)  
#'
#'    # change elem. colors
#'    biplotPseudo3d(boeker, e.color=c("white", "darkgreen"))
#'    # change con. colors 
#'    biplotPseudo3d(boeker, c.color=c("white", "darkgreen")) 
#'    # change color mapping range
#'    biplotPseudo3d(boeker, c.colors.map=c(0, 1))            
#'
#'    # set uniform con. text size
#'    biplotPseudo3d(boeker, c.cex=1)     
#'    # change text size mapping range        
#'    biplotPseudo3d(boeker, c.cex=c(.4, 1.2))    
#' }
#'
biplotPseudo3d <- function( x, dim=1:2, map.dim=3, 
                            e.point.col=c("white", "black"),
                            e.point.cex=c(.6, 1.2),
                            e.label.col=c("white", "black"),
                            e.label.cex=c(.6, .8),   
                            e.color.map=c(.4, 1),               
                            c.point.col=c("white", "darkred"),
                            c.point.cex=c(.6, 1.2),
                            c.label.col=c("white", "darkred"),
                            c.label.cex=c(.6, .8),
                            c.color.map=c(.4, 1),
                            ...)
{            
  biplot2d(x=x, dim=dim, map.dim=map.dim, 
            e.point.col=e.point.col,
            e.point.cex=e.point.cex,
            e.label.col=e.label.col,
            e.label.cex=e.label.cex, 
            e.color.map=e.color.map,                                   
            c.point.col=c.point.col,
            c.point.cex=c.point.cex,
            c.label.col=c.label.col,
            c.label.cex=c.label.cex,
            c.color.map=c.color.map, 
            ...)
}

    
#' Draws Slater's INGRID biplot in 2D. 
#'
#' The default is to use row centering 
#' and no normalization. Note that Slater's biplot is just a 
#' special case of a biplot
#' that can be produced using the \code{\link{biplot2d}} function with the arguments
#' \code{center=1, g=1, h=1}. The arguments that can be used in this function
#' are the same as in \code{\link{biplot2d}}. 
#' Here, only the arguments that are set for Slater's biplot are described.
#' To see all the parameters that can be changed see \code{\link{biplot2d}}.
#'
#' @param x           \code{repgrid} object.
#' @param  center		  Numeric. The type of centering to be performed. 
#'                    0= no centering, 1= row mean centering (construct), 
#'                    2= column mean centering (elements), 3= double-centering (construct and element means),
#'                    4= midpoint centering of rows (constructs).
#'                    Slater's biplot uses \code{1} (row centering).
#' @param g           Power of the singular value matrix assigned to the left singular 
#'                    vectors, i.e. the constructs.
#' @param h           Power of the singular value matrix assigned to the right singular 
#'                    vectors, i.e. the elements.
#' @param ...         Additional parameters for be passed to \code{\link{biplot2d}}.
#'
#' @author            Mark Heckmann
#' @export
#'
#' @seealso   Unsophisticated biplot: \code{\link{biplotSimple}}; \cr
#'            2D biplots:
#'            \code{\link{biplot2d}},
#'            \code{\link{biplotEsa2d}},
#'            \code{\link{biplotSlater2d}};\cr
#'            Pseudo 3D biplots:
#'            \code{\link{biplotPseudo3d}},  
#'            \code{\link{biplotEsaPseudo3d}},
#'            \code{\link{biplotSlaterPseudo3d}};\cr
#'            Interactive 3D biplots:
#'            \code{\link{biplot3d}},
#'            \code{\link{biplotEsa3d}},
#'            \code{\link{biplotSlater3d}};\cr
#'            Function to set view in 3D:
#'            \code{\link{home}}.
#'
#' @examples \dontrun{
#'    # See examples in \code{\link{biplot2d}} as the same arguments
#'    # can used for this function.
#' }
#'
biplotSlater2d <- function(x, center=1, g=1, h=1, ...){
  biplot2d(x=x, center=center, g=g, h=h, ...)
}


#' Draws Slater's biplot in 2D with depth impression (pseudo 3D).
#'
#' The default is to use row centering 
#' and no normalization. Note that Slater's biplot is just a special 
#' case of a biplot that can be produced using the \code{\link{biplotPseudo3d}} 
#' function with the arguments \code{center=1, g=1, h=1}.
#' Here, only the arguments that are modified for Slater's biplot are described.
#' To see all the parameters that can be changed see \code{\link{biplot2d}}
#' and \code{\link{biplotPseudo3d}}.
#'
#' @param x           \code{repgrid} object.
#' @param  center		  Numeric. The type of centering to be performed. 
#'                    0= no centering, 1= row mean centering (construct), 
#'                    2= column mean centering (elements), 3= double-centering (construct and element means),
#'                    4= midpoint centering of rows (constructs).
#'                    Slater's biplot uses \code{1} (row centering).
#' @param g           Power of the singular value matrix assigned to the left singular 
#'                    vectors, i.e. the constructs.
#' @param h           Power of the singular value matrix assigned to the right singular 
#'                    vectors, i.e. the elements.
#' @param ...         Additional parameters for be passed to \code{\link{biplotPseudo3d}}.
#'
#' @author            Mark Heckmann
#' @export
#'
#' @seealso   Unsophisticated biplot: \code{\link{biplotSimple}}; \cr
#'            2D biplots:
#'            \code{\link{biplot2d}},
#'            \code{\link{biplotEsa2d}},
#'            \code{\link{biplotSlater2d}};\cr
#'            Pseudo 3D biplots:
#'            \code{\link{biplotPseudo3d}},  
#'            \code{\link{biplotEsaPseudo3d}},
#'            \code{\link{biplotSlaterPseudo3d}};\cr
#'            Interactive 3D biplots:
#'            \code{\link{biplot3d}},
#'            \code{\link{biplotEsa3d}},
#'            \code{\link{biplotSlater3d}};\cr
#'            Function to set view in 3D:
#'            \code{\link{home}}.
#'
#' @examples \dontrun{
#'    # See examples in \code{\link{biplotPseudo3d}} as the same arguments
#'    # can used for this function.
#' }
#'
biplotSlaterPseudo3d <- function(x, center=1, g=1, h=1, ...){
  biplotPseudo3d(x=x, center=center, g=g, h=h, ...)
}



#' Plot an eigenstructure analysis (ESA) biplot in 2D.
#' 
#' The ESA is a special type of biplot suggested by Raeithel (e.g. 1998).
#' It uses midpoint centering as a default. Note that the eigenstructure analysis
#' is just a special case of a biplot that can also be produced using the 
#' \code{\link{biplot2d}} function with the arguments 
#' \code{center=4, g=1, h=1}.
#' Here, only the arguments that are modified for the ESA biplot are described.
#' To see all the parameters that can be changed see \code{\link{biplot2d}}.
#'
#' @param x           \code{repgrid} object.
#' @param  center		  Numeric. The type of centering to be performed. 
#'                    0= no centering, 1= row mean centering (construct), 
#'                    2= column mean centering (elements), 3= double-centering (construct and element means),
#'                    4= midpoint centering of rows (constructs).
#'                    Eigenstructure analyis uses midpoint centering (\code{4}).
#' @param g           Power of the singular value matrix assigned to the left singular 
#'                    vectors, i.e. the constructs. Eigenstructure analyis uses  
#'                    \code{g=1}.
#' @param h           Power of the singular value matrix assigned to the right singular 
#'                    vectors, i.e. the elements. Eigenstructure analyis uses  
#'                    \code{h=1}.
#' @param ...         Additional parameters for be passed to \code{\link{biplot2d}}.
#'
#' @references   Raeithel, A. (1998). Kooperative Modellproduktion von Professionellen 
#'                und Klienten. Erlaeutert am Beispiel des Repertory Grid.
#'                In A. Raeithel (1998). Selbstorganisation, Kooperation, 
#'                Zeichenprozess. Arbeiten zu einer kulturwissenschaftlichen, 
#'                anwendungsbezogenen Psychologie (p. 209-254). Opladen: 
#'                Westdeutscher Verlag.
#'
#' @author        Mark Heckmann
#' @export
#'
#' @seealso   Unsophisticated biplot: \code{\link{biplotSimple}}; \cr
#'            2D biplots:
#'            \code{\link{biplot2d}},
#'            \code{\link{biplotEsa2d}},
#'            \code{\link{biplotSlater2d}};\cr
#'            Pseudo 3D biplots:
#'            \code{\link{biplotPseudo3d}},  
#'            \code{\link{biplotEsaPseudo3d}},
#'            \code{\link{biplotSlaterPseudo3d}};\cr
#'            Interactive 3D biplots:
#'            \code{\link{biplot3d}},
#'            \code{\link{biplotEsa3d}},
#'            \code{\link{biplotSlater3d}};\cr
#'            Function to set view in 3D:
#'            \code{\link{home}}.
#'
#' @examples \dontrun{
#'    # See examples in \code{\link{biplot2d}} as the same arguments
#'    # can used for this function.
#' }
#'
biplotEsa2d <- function(x, center=4, g=1, h=1, ...){
  biplot2d(x=x, center=center, g=g, h=h, ...)
}


#' Plot an eigenstructure analysis (ESA) in 2D grid with 3D 
#' impression (pseudo 3D). 
#'
#' The ESA is 
#' a special type of biplot suggested by Raeithel (e.g. 1998).
#' It uses midpoint centering as a default. Note that the eigenstructure analysis
#' is just a special case of a biplot that can also be produced using the 
#' \code{\link{biplot2d}} function with the arguments 
#' \code{center=4, g=1, h=1}.
#' Here, only the arguments that are modified for the ESA biplot are described.
#' To see all the parameters that can be changed see \code{\link{biplot2d}}
#' and \code{\link{biplotPseudo3d}}.
#'
#' @param x           \code{repgrid} object.
#' @param  center		  Numeric. The type of centering to be performed. 
#'                    0= no centering, 1= row mean centering (construct), 
#'                    2= column mean centering (elements), 3= double-centering 
#'                    (construct and element means),
#'                    4= midpoint centering of rows (constructs).
#'                    Eigenstructure analyis uses midpoint centering (\code{4}).
#' @param g           Power of the singular value matrix assigned to the left singular 
#'                    vectors, i.e. the constructs. Eigenstructure analyis uses  
#'                    \code{g=1}.
#' @param h           Power of the singular value matrix assigned to the right singular 
#'                    vectors, i.e. the elements. Eigenstructure analyis uses  
#'                    \code{h=1}.
#' @param ...         Additional parameters for be passed to \code{\link{biplotPseudo3d}}.
#'
#' @author            Mark Heckmann
#' @export
#'
#' @seealso   Unsophisticated biplot: \code{\link{biplotSimple}}; \cr
#'            2D biplots:
#'            \code{\link{biplot2d}},
#'            \code{\link{biplotEsa2d}},
#'            \code{\link{biplotSlater2d}};\cr
#'            Pseudo 3D biplots:
#'            \code{\link{biplotPseudo3d}},  
#'            \code{\link{biplotEsaPseudo3d}},
#'            \code{\link{biplotSlaterPseudo3d}};\cr
#'            Interactive 3D biplots:
#'            \code{\link{biplot3d}},
#'            \code{\link{biplotEsa3d}},
#'            \code{\link{biplotSlater3d}};\cr
#'            Function to set view in 3D:
#'            \code{\link{home}}.
#'
#' @examples \dontrun{
#'    # See examples in \code{\link{biplotPseudo3d}} as the same arguments
#'    # can used for this function.
#' }
#'
biplotEsaPseudo3d <- function(x, center=4, g=1, h=1, ...){
  biplotPseudo3d(x=x, center=center, g=g, h=h, ...)
}




###############################################################################

# x <- boeker
# x <- calcBiplotCoords(x, g=1, h=1)
# x <- prepareBiplotData(x, unity=F)
# biplot2d(x)
# 
# biplot2d(x)


###############################################################################





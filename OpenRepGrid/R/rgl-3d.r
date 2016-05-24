#' Draw standard axes in the origin in an rgl plot.
#'
#' @param max.dim     maximum length of axis.
#' @param lwd         line width.
#' @param a.cex       cex for axis labels.
#' @param a.col       axis color.
#' @param a.radius      radius of spheres at the end points of the axes.
#' @param labels      logical. wether to draw axis labels.
#' @param spheres     logical. wether to draw axis spheres at the end points.  
#' @param ...         not evaluated.
#'
#' @author  Mark Heckmann
#' @export
#' @keywords internal
#'
rglDrawStandardAxes <- function(max.dim=1, lwd=1, a.cex=1.1, a.col="black",
                                a.radius=.05, labels=TRUE, spheres=FALSE, ...){
  lines3d(c(0, max.dim), c(0,0), c(0,0), lwd=lwd, col=a.col)
  lines3d(c(0,0), c(0, max.dim), c(0,0), lwd=lwd, col=a.col)
  lines3d(c(0,0), c(0,0), c(0, max.dim), lwd=lwd, col=a.col)
  if (labels){
    text3d(max.dim, 0, 0, "X", cex=a.cex, adj=c(1,1), col=a.col)
    text3d(0, max.dim, 0, "Y", cex=a.cex, adj=c(1,1), col=a.col)
    text3d(0, 0, max.dim, "Z", cex=a.cex, adj=c(1,1), col=a.col)
  }
  if (spheres){
    spheres3d(max.dim, 0, 0, radius=a.radius, col=a.col)
    spheres3d(0, max.dim, 0, radius=a.radius, col=a.col)
    spheres3d(0, 0, max.dim, radius=a.radius, col=a.col)    
  }
}
# rgl.open()
# rgl.points(rnorm(1000), rnorm(1000), rnorm(1000), color=heat.colors(1000))
# rglDrawStandardAxes(3)


#' Draw standard ellipses in the origin in an rgl plot.
#'
#' @param max.dim   soon
#' @param lwd       soon
#' @param col       soon
#'
#' @author  Mark Heckmann
#' @export
#' @keywords internal
#'
rglDrawStandardEllipses <- function(max.dim=1, lwd=1, col="black"){
  x <- seq(0, 2*pi, len=361)
  x <- data.frame(sin(x), cos(x)) * max.dim
  lines3d(x[,1], x[,2], 0, col=col, lwd=lwd)
  lines3d(x[,1], 0, x[,2], col=col, lwd=lwd)
  lines3d(0, x[,1], x[,2], col=col, lwd=lwd)
}



rglDrawElementPoints <- function(coords, dim=1:3, e.radius=.1, e.sphere.col="black", ...){
  coords <- coords[ ,dim]
  spheres3d(coords[,1], coords[,2], coords[,3], 
            radius=e.radius, color=e.sphere.col, aspect=F)
}


rglDrawElementLabels <- function(coords, labels=FALSE, dim=1:3, e.radius=.1, e.cex=.6, e.text.col="black", ...){
  coords <- coords[ ,dim]
  if (!identical(labels, FALSE)){  
    coords.text <- coords - e.radius/2     # offset text for elements
    texts3d(x= coords.text[,1], 
            y= coords.text[,2], 
            z= coords.text[,3], 
            texts=labels, adj=c(1,1), cex=e.cex, col=e.text.col, aspect =F )
  } 
}


#' draw constructs in rgl
#'
#' @param coords          coordinates for construct points.
#' @param dim             dimensions of coordinates to use.
#' @param c.radius        radius of construct spheres.
#' @param c.sphere.col    color of construct spheres.
#' @param ...             not evaluated.
#'
#' @author  Mark Heckmann
#' @export
#' @keywords internal
#'
rglDrawConstructPoints <- function(coords, dim=1:3, c.radius=.02, c.sphere.col=grey(.4), 
                                   ...){
  coords <- coords[ ,dim]
  coords[is.na(coords)] <- 0    # replace NAs by zero, so Na can be entered as dim for 2d projection
  spheres3d(coords[, dim], radius=c.radius, color=c.sphere.col)
}

#' draw constructs in rgl
#'
#' @param coords      coordinates for constructs labels.
#' @param labels      labels for constructs.
#' @param dim         dimensions of coordinates to use.
#' @param c.cex       cex for construct text.
#' @param c.text.col  color for construct text.
#' @param ...         not evaluated.
#'
#' @author  Mark Heckmann
#' @export
#' @keywords internal
#'
rglDrawConstructLabels <- function(coords, labels=FALSE, dim=1:3,   
                              c.cex=.6, c.text.col=grey(.4), ...){
  coords <- coords[ ,dim]
  coords[is.na(coords)] <- 0        # replace NAs by zero, so Na can be entered as dim for 2d projection
  if (!identical(labels, FALSE)){
    texts3d(coords, texts=labels, adj=c(.5,.5), 
            cex=c.cex, col=c.text.col, aspect=F)  
  }
}


#' biplot3dBase2 is the workhorse to draw a grid in rgl (3D device).
#'
#' @param x               \code{repgrid} object.
#' @param dim             Dimensions to display.
#' @param labels.e        Logical. whether element labels are displayed.
#' @param labels.c        Logical. whether construct labels are displayed.
#' @param lines.c         Numeric. The way lines are drawn through the construct vectors.
#'                        \code{0 =} no lines, \code{1 =} lines from constructs to outer frame,
#'                        \code{2 =} lines from the center to outer frame.
#' @param lef             Construct lines extension factor.
#' @param alpha.sphere    Numeric. alpha blending of the sourrounding sphere (default\code{".05"}).
#' @param col.sphere      Color of sourrouding sphere (default\code{"black"}).
#' @param ext.sphere      Extension factor for sphere
#' @param col.frame       Color of the sourrounding frame.
#' @param zoom            Not yet used. Scaling factor for all vectors. Can be used to zoom
#'                        the plot in and out (default \code{1}). 
#' @param ...             Parameters to be passed on.  
#'
#' @author  Mark Heckmann
#' @keywords internal
#' @export
#'
biplot3dBase2 <- function(x, dim=1:3, labels.e=TRUE, labels.c=TRUE, lines.c=1, 
                     lef=1.1, frame=1, col.frame=grey(.6), 
                     col.sphere ="black", alpha.sphere=.05, zoom=1,
#                           c.points.show=TRUE,
#                           c.labels.show=TRUE,
#                          e.points.show=TRUE,
#                          e.labels.show=TRUE,      
                     ...)
{
  x <- calcBiplotCoords(x, ...)
  x <- prepareBiplotData(x, ...)
  
  showpoint <- showlabel <- type <- NULL          # to prevent 'R CMD check' from noting a missing binding 
                                                  # as the variables are provided in object x as default
  open3d()                                        # open rgl device
  par3d(params=list(
        windowRect=c(100,100,600,600)))           # enlarge and position 3d device
  view3d(theta = 0, phi = 0, zoom=.6)             # change 3d view angle
  rgl.bg(color="white")                           # set background color
     
  # select spheres to draw and labels to show
  # select which elements to show
  if (identical(labels.e, TRUE)) 
    labels.e <- getElementNames(x)  
  if (identical(labels.c, TRUE)){
    labels.l <- getConstructNames(x)$l
    labels.r <- getConstructNames(x)$r
  } else {
    labels.r <- FALSE
    labels.l <- FALSE
  }

  X <- x@calcs$biplot$X             # pre-transformed (centered etc.) grid matrix
  Eu <- x@calcs$biplot$e.unity      # element coordinates scaled/unified
  Cu <- x@calcs$biplot$c.unity      # construct coordinates scaled/unified
  
  pdat <- x@plotdata                # plotdata prepared by prepareBiplotData()
  
  mval <- max(abs(rbind(Eu[, dim], Cu[ ,dim])),   # get maximum value of construct and element coordinates
              na.rm=TRUE) 
  
  #Eu <- Eu * zoom
  #Cu <- Cu * zoom
  
  # prolongation of construct vector to outer side
  Cu.norm <- apply(Cu[, dim]^2, 1, sum, na.rm=TRUE)^.5 
  Cup <- Cu[, dim] / Cu.norm * (lef * mval)
    
  # plot element spheres
  es.p <- subset(pdat, type=="e" & showpoint==T)
  rglDrawElementPoints(es.p[c("x", "y", "z")], e.radius=mval/50, ...)
  # labels for elements 
  es.l <- subset(pdat, type=="e" & showlabel==T & showpoint==T)
  rglDrawElementLabels(es.l[c("x", "y", "z")], labels=es.l$label, e.radius=mval/50, ...)
  
  standardizeCoords <- function(x, dim=1:3){
    x.norm <- apply(x[, dim]^2, 1, sum, na.rm=TRUE)^.5 
    xsc <- x[, dim] / x.norm * (lef * mval)
    xsc
  }
    
  # make construct spheres
  cs.p <- subset(pdat, type %in% c("cl", "cr") & showpoint==T)
  cs.p.xyz <- cs.p[c("x", "y", "z")] 
  # labels for constructs
  cs.l <- subset(pdat, type %in% c("cl", "cr") & showlabel==T)
  cl.l.xyz <- cs.l[c("x", "y", "z")] 
  cl.l.xyz.outer <- standardizeCoords(cs.l[c("x", "y", "z")])
  
  if (lines.c == 0){           # no construct lines labels at cons pos
    rglDrawConstructLabels(cl.l.xyz, labels=cs.l$label, ...)    
    rglDrawStandardAxes(mval, spheres=F)
    #rglDrawConstructLabels(Cu[, dim], labels=labels.r, ...)
    #rglDrawConstructLabels(-Cu[, dim], labels=labels.l, ...)
  } else if (lines.c == 1){     # construct lines from cons pos to outside
    segments3d(interleave(cl.l.xyz, cl.l.xyz.outer), col="grey")
    rglDrawConstructLabels(cl.l.xyz.outer, labels=cs.l$label, ...)        
    rglDrawStandardAxes(lef * mval, a.col="black")
    #segments3d(interleave(-Cu[, dim], -Cup), col="grey")       # Cu and Cup from older implementation without use if x@plotdata
    #rglDrawConstructLabels(Cup, labels=labels.r, ...)  
    #rglDrawConstructLabels(-Cup, labels=labels.l, ...)
  } else if (lines.c == 2){     # construct lines from center to outside
    nm <- matrix(0, ncol=3, nrow=nrow(cl.l.xyz.outer))
    segments3d(interleave(nm, as.matrix(cl.l.xyz.outer)), col="grey")
    rglDrawConstructLabels(cl.l.xyz.outer, labels=cs.l$label, ...)        
    rglDrawStandardAxes(lef * mval, a.col="black")
  } else {
    stop("'lines.c' can only take numeric values from 0 to 2")
  }
  rglDrawConstructPoints(cs.p.xyz, c.radius=mval/200, ...)
  #rglDrawConstructPoints(-Cu[, dim], c.radius=mval/200, ...)
  #rglDrawStandardEllipses(max.dim)
  
  # trick to make user coordinate system's origin the center of rotation
  mval <- max(abs(par3d()$bbox))                   # get max value in x,y,z
  ps <- interleave(mval*diag(3), -mval*(diag(3)))     
  spheres3d(ps, radius=0)                          # draw invisible spheres at the extremes
  
  # select type of frame ariound the whole plot
  # 0=none, 1= simple box, 2= box with grid, 3=sphere.
  if (frame == 1){
    # make box around device
    ss <- matrix(c(  mval,  mval,  mval,    # top
                    -mval,  mval,  mval,
                    -mval,  mval,  mval,
                    -mval, -mval,  mval,
                    -mval, -mval,  mval,
                     mval, -mval,  mval,
                     mval, -mval,  mval,
                     mval,  mval,  mval,
                     mval,  mval, -mval,    # bottom                
                    -mval,  mval, -mval,
                    -mval,  mval, -mval,
                    -mval, -mval, -mval,
                    -mval, -mval, -mval,
                     mval, -mval, -mval,
                     mval, -mval, -mval,
                     mval,  mval, -mval,
                     mval,  mval,  mval,    # sides 
                     mval,  mval, -mval,
                    -mval,  mval,  mval,
                    -mval,  mval, -mval,
                    -mval, -mval,  mval,
                    -mval, -mval, -mval,
                     mval, -mval,  mval,
                     mval, -mval, -mval), ncol=3, byrow=T)
    segments3d(ss, col=col.frame)
  } else if (frame == 2){
    grid3d(c("x+","x-", "y+", "y-", "z+", "z-"))    
  } else if (frame == 3){
    # sphere for easier 3D impression if prompted 
    spheres3d(0, 0, 0, radius=mval, color=col.sphere, 
              alpha=alpha.sphere, aspect=F, front="lines", back="lines")
  }
}



#' Draw grid in rgl (3D device). 
#'
#' The 3D biplot opens an interactive 
#' 3D device that can be rotated and zoomed using the mouse. 
#' A 3D device facilitates the exploration of grid data as 
#' significant proportions of the sum-of-squares are often 
#' represented beyond the first two dimensions. Also, in a lot of 
#' cases it may be of interest to explore the grid space from 
#' a certain angle, e.g. to gain an optimal view onto the set 
#' of elements under investigation (e.g. Raeithel, 1998). 
#'
#' @param x             \code{repgrid} object.
#' @param dim           Dimensions to display.
#' @param labels.e      Logical. whether element labels are displayed.
#' @param labels.c      Logical. whether construct labels are displayed.
#' @param lines.c       Numeric. The way lines are drawn through the construct vectors.
#'                      \code{0 =} no lines, \code{1 =} lines from constructs to outer frame,
#'                      \code{2 =} lines from the center to outer frame.
#' @param lef           Construct lines extension factor
#' 
#' @param  center		    Numeric. The type of centering to be performed. 
#'                      0= no centering, 1= row mean centering (construct), 
#'                      2= column mean centering (elements), 3= double-centering (construct and element means),
#'                      4= midpoint centering of rows (constructs).
#'                      Default is \code{1} (row centering).
#'
#' @param  normalize    A numeric value indicating along what direction (rows, columns)
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
#'
#' @param c.sphere.col  Color of construct spheres.
#' @param c.cex         Size of construct text.
#' @param c.text.col    Color for construct text.
#'
#' @param e.sphere.col  Color of elements.
#' @param e.cex         Size of element labels.
#' @param e.text.col    Color of element labels.
#'
#' @param alpha.sphere  Numeric. alpha blending of the sourrounding sphere (default\code{".05"}).
#' @param col.sphere    Color of sourrouding sphere (default\code{"black"}).
#'
#' @param unity         Scale elements and constructs coordinates to unit scale (maximum of 1)
#'                      so they are printed more neatly (default \code{TRUE}).
#' @param unity3d       To come.
#' @param scale.e       Scaling factor for element vectors. Will cause element points to move a bit more
#'                      to the center (but only if \code{unity} or \code{unity3d} is \code{TRUE}).
#'                      This argument is for visual appeal only.
#' @param zoom          Not yet used. Scaling factor for all vectors. Can be used to zoom
#'                      the plot in and out (default \code{1}).
#' @param ...           Parameters to be passed on.
#'
#' @author  Mark Heckmann
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
#' @references      Raeithel, A. (1998). Kooperative Modellproduktion von 
#'                  Professionellen und Klienten - erlauetert am Beispiel des 
#'                  Repertory Grid. \emph{Selbstorganisation, Kooperation, Zeichenprozess: 
#'                  Arbeiten zu einer kulturwissenschaftlichen, anwendungsbezogenen 
#'                  Psychologie} (pp. 209-254). Opladen: Westdeutscher Verlag.
#'
#' @examples \dontrun{
#'
#'    biplot3d(boeker)
#'    biplot3d(boeker, unity3d=T)
#'
#'    biplot3d(boeker, e.sphere.col="red",
#'             c.text.col="blue")
#'    biplot3d(boeker, e.cex=1)
#'    biplot3d(boeker, col.sphere="red")
#'
#'    biplot3d(boeker, g=1, h=1)    # INGRID biplot
#'    biplot3d(boeker, g=1, h=1,    # ESA biplot
#'             center=4)
#' }
#'
biplot3d <- function(x, dim=1:3, labels.e=TRUE, labels.c=TRUE, lines.c=TRUE, 
                     lef=1.3, center=1, normalize=0, g=0, h=1, col.active=NA, 
                     col.passive=NA, 
                     c.sphere.col =grey(.4), c.cex=.6, c.text.col=grey(.4),
                     e.sphere.col =grey(0), e.cex=.6, e.text.col=grey(0),
                     alpha.sphere=.05, col.sphere="black", 
                     unity=FALSE, 
                     unity3d=FALSE,
                     scale.e=.9, zoom=1, ...)
{
    biplot3dBase2(x=x, dim=dim, labels.e=labels.e, labels.c=labels.c, lines.c=lines.c, 
                lef=lef, center=center, normalize=normalize, g=g, h=h, 
                col.active=col.active, col.passive=col.passive,  
                c.sphere.col =c.sphere.col, c.cex=c.cex, c.text.col=c.text.col,
                e.sphere.col =e.sphere.col, e.cex=e.cex, e.text.col=e.text.col, 
                alpha.sphere=alpha.sphere, col.sphere=col.sphere, 
                unity=unity, unity3d=unity3d, scale.e=scale.e, zoom=zoom, ...)
}
       
  
#' Draw the Slater's INGRID biplot in rgl (3D device).
#'
#' The 3D biplot opens an interactive 
#' 3D device that can be rotated and zoomed using the mouse. 
#' A 3D device facilitates the exploration of grid data as 
#' significant proportions of the sum-of-squares are often 
#' represented beyond the first two dimensions. Also, in a lot of 
#' cases it may be of interest to explore the grid space from 
#' a certain angle, e.g. to gain an optimal view onto the set 
#' of elements under investigation (e.g. Raeithel, 1998).
#' Note that Slater's biplot is just a special case of a biplot
#' that can be produced using the \code{\link{biplot3d}} 
#' function with the arguments \code{center=1, g=1, h=1}.
#'
#' @param x             \code{repgrid} object.
#' @param  center		    Numeric. The type of centering to be performed. 
#'                      0= no centering, 1= row mean centering (construct), 
#'                      2= column mean centering (elements), 3= double-centering (construct and element means),
#'                      4= midpoint centering of rows (constructs).
#'                      Default is \code{1} (row i.e. construct centering).
#' @param g             Power of the singular value matrix assigned to the left singular 
#'                      vectors, i.e. the constructs.
#' @param h             Power of the singular value matrix assigned to the right singular 
#'                      vectors, i.e. the elements.
#' @param ...           Additional arguments to be passed to biplot3d.
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
#'    biplotSlater3d(boeker)
#'    biplotSlater3d(boeker, unity3d=T)
#'
#'    biplotSlater3d(boeker, e.sphere.col="red",
#'                   c.text.col="blue")
#'    biplotSlater3d(boeker, e.cex=1)
#'    biplotSlater3d(boeker, col.sphere="red")
#'
#' }
#'
biplotSlater3d <- function(x, center=1, g=1, h=1, ...){
  biplot3d(x=x, center=center, g=g, h=h, ...)
}                      


#' Draw the eigenstructure analysis (ESA) biplot in rgl (3D device). 
#'
#' The 3D biplot opens an interactive 
#' 3D device that can be rotated and zoomed using the mouse. 
#' A 3D device facilitates the exploration of grid data as 
#' significant proportions of the sum-of-squares are often 
#' represented beyond the first two dimensions. Also, in a lot of 
#' cases it may be of interest to explore the grid space from 
#' a certain angle, e.g. to gain an optimal view onto the set 
#' of elements under investigation (e.g. Raeithel, 1998).
#' Note that the eigenstructure analysisis just a special case 
#' of a biplot that can also be produced using the 
#' \code{\link{biplot3d}} function with the arguments 
#' \code{center=4, g=1, h=1}.
#'
#' @param x             \code{repgrid} object.
#' @param  center		    Numeric. The type of centering to be performed. 
#'                      0= no centering, 1= row mean centering (construct), 
#'                      2= column mean centering (elements), 3= double-centering (construct and element means),
#'                      4= midpoint centering of rows (constructs).
#'                      Default is \code{4} (scale midpoint centering).
#' @param g             Power of the singular value matrix assigned to the left singular 
#'                      vectors, i.e. the constructs.
#' @param h             Power of the singular value matrix assigned to the right singular 
#'                      vectors, i.e. the elements.
#' @param ...           Additional arguments to be passed to \code{\link{biplot3d}}.
#'
#' @author  Mark Heckmann
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
#'    biplotEsa3d(boeker)
#'    biplotEsa3d(boeker, unity3d=T)
#'
#'    biplotEsa3d(boeker, e.sphere.col="red",
#'                c.text.col="blue")
#'    biplotEsa3d(boeker, e.cex=1)
#'    biplotEsa3d(boeker, col.sphere="red")
#'
#' }
#'
biplotEsa3d <- function(x, center=1, g=1, h=1, ...){
  biplot3d(x=x, center=center, g=g, h=h, ...)
}                      
 

#' Rotate the interactive 3D device to default views.
#'
#' Rotate the interactive 3D device to a default viewpoint or
#' to a position defined by \code{theta} and \code{phi} in Euler angles.
#' Three default viewpoints are implemented rendering a view 
#' so that two axes span a plane and the third axis is 
#' poiting out of the screen.
#' 
#' @param view    Numeric. Specifying one of three default views.
#'                1 = XY, 2=XZ and 3=YZ-plane.
#' @param theta   Numeric. Euler angle. Overrides view setting.
#' @param phi     Numeric. Euler angle. Overrides view setting.
#'
#' return \code{NULL}.
#'
#' @author  Mark Heckmann
#' @export
#'
#' @seealso   Interactive 3D biplots:
#'            \code{\link{biplot3d}},      
#'            \code{\link{biplotSlater3d}},
#'            \code{\link{biplotEsa3d}}.
#'
#' @examples \dontrun{
#'
#'    biplot3d(boeker)
#'    home(2)
#'    home(3)
#'    home(1)
#'    home(theta=45, phi=45)
#'
#' }
#'
home <- function(view=1, theta=NULL, phi=NULL){
  if (!view %in% 1:3)
    stop("'view' must take a numeric value between 1 and 3")
  p3d <- par3d()
  if (is.null(theta) & is.null(phi)){
    if (view == 1){
      theta <- 0; phi <- 0
    } else if (view == 2){
      theta <- 0; phi <- 90
    } else if (view == 3){
      theta <- 90; phi <- 0
    }  
  }
  view3d(theta = theta, phi = phi, zoom=p3d$zoom)           # change 3d view angle 
}




###############################################################################
###                              EXAMPLES                                   ###
###############################################################################

# biplot3d(raeithel, labels.c=F)
# 
# x <- raeithel
# x <- calcBiplotCoords(x, g=0, h=1, midp=T, col.active=c(2,4,10))
# x <- prepareBiplotData(x, unity=T)
# biplot3d(x)

# 
# M <- par3d("userMatrix")              # get current position matrix
# dir <- "/Users/markheckmann/Documents/Magic Briefcase/DA openRepgrid/openrepgrid/basic/output/animation"
# M1 <- rotate3d(M, pi/2, 1, 0, 0)
# M2 <- rotate3d(M1, pi/2, 0, 0, 1)
# movie3d(par3dinterp( userMatrix=list(M, M1, M2, M1, M), method="linear"), 
#         duration=4, fps=20, convert=F, clean=F, dir=dir)


# open3d()
# lines3d(c(0, 1), c(0,0), c(0,0))
# lines3d(c(0,0), c(0, 1), c(0,0))
# lines3d(c(0,0), c(0,0), c(0, 1))




# mouseTrackballOrigin <- function(button = 1, dev = rgl.cur(), origin=c(0,0,0) ) {
#    width <- height <- rotBase <- NULL
#    userMatrix <- list()
#    cur <- rgl.cur()
#    offset <- NULL
#    scale <- NULL
# 
#    screenToVector <- function(x, y) {
#      radius <- max(width, height)/2
#      centre <- c(width, height)/2
#      pt <- (c(x, y) - centre)/radius
#      len <- vlen(pt)
# 
#      if (len > 1.e-6) pt <- pt/len
# 
#      maxlen <- sqrt(2)
#      angle <- (maxlen - len)/maxlen*pi/2
#      z <- sin(angle)
#      len <- sqrt(1 - z^2)
#      pt <- pt * len
#      return (c(pt, z))
#    }
# 
#    trackballBegin <- function(x, y) {
#        vp <- par3d("viewport")
#        width <<- vp[3]
#        height <<- vp[4]
#        cur <<- rgl.cur()
#        bbox <- par3d("bbox")
#        center <- c(sum(bbox[1:2])/2, sum(bbox[3:4])/2, sum(bbox[5:6])/2)
#        scale <<- par3d("scale")
#        offset <<- (center - origin)*scale
#        for (i in dev) {
#            if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
#            else userMatrix[[i]] <<- par3d("userMatrix")
#        }
#        rgl.set(cur, TRUE)
#        rotBase <<- screenToVector(x, height - y)
#    }
# 
#    trackballUpdate <- function(x,y) {
#        rotCurrent <- screenToVector(x, height - y)
#        angle <- angle(rotBase, rotCurrent)
#        axis <- xprod(rotBase, rotCurrent)
#        mouseMatrix <- rotationMatrix(angle, axis[1], axis[2], axis[3])
#        for (i in dev) {
#            if (inherits(try(rgl.set(i, TRUE)), "try-error")) dev <<- dev[dev != i]
#            else par3d(userMatrix = t(translationMatrix(-offset[1], -offset[2], -offset[3])) %*% mouseMatrix  %*% t(translationMatrix(offset[1], offset[2], offset[3])) %*%userMatrix[[i]])
#        }
#        rgl.set(cur, TRUE)
#    }
# 
#    for (i in dev) {
#        rgl.set(i, TRUE)
#        rgl.setMouseCallbacks(button, begin = trackballBegin, update = trackballUpdate, end = NULL)
#    }
#    rgl.set(cur, TRUE)
# }

# additioally load functions from demo(). see email from Duncan Murdoch 25.04.2011
# mouseTrackballOrigin()


#########################################################################################
# TODO: rotations of the biplot
# 
# eulerxyz <- function(phi, theta, psi){
#   phi <- phi*180/pi     # conversion from degree to radians
#   theta <- theta*180/pi
#   psi <- psi*180/pi
#   
#   matrix(c(cos(theta)*cos(psi), -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi), sin(phi)*sin(psi)+ cos(phi)*sin(theta)*cos(psi), 
#            cos(theta)*sin(psi),  cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi), -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi),
#           -sin(theta)         ,  sin(phi)*cos(theta)                             , cos(phi)*cos(theta)), ncol=3)
# }
# 
# m <- par3d()$userMatrix[1:3, 1:3]



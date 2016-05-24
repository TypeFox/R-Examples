##########################################################################
#
#  Routines to draw graphs.
#
#  These routines have been extracted from the SNA package (mainly gplot)
#  and reorganized to exploit more conveniently its possibilies.
#
#  Author : G.Grasseau (Statistics and Genome laboratory)
#
#  Entry points:
#  ------------
#
#    Gplot.graphics:     set default graphic parameters.
#                        Main steps:
#                        + store initial matrix as boolean one,
#                        + build the structure (list) which handle all
#                        graphics parameters.
#
#    Gplot.network:      specify the network representation (up to now only
#                        undirect graphs are taken into account).
#                        Main steps:
#                        + compute vertex position,
#                        + graph box limits and scaling factor,
#                        + define which vertices to display.
#
#    Gplot.vertex:       draw vertices.
#
#    Gplot.vertex.label: draw vertex labels.
#
#    Gplot.edge:         draw edges
#                        Main steps:
#                        + identify edges to be drawn.
#                        + remove loops
#                        + draw edges/arrows.
#
#  Utilities:
#  ---------
#   + draw.edge:        (used by Gplot.edge) draw edges/arrows
#   + isolate.vertices: (used by Gplot.network) tag isolate vertices (TRUE)
#
#  To improve:
#  ----------
#   + sides of some boxed labels are not displayed
#   + 2 variables (displayisolates and use.isolates) should be introduced
#     to avoid taking into account of isolate vertices in the vertex coordinate
#     computation (4 cases have to be studied).
#
#  To implement:
#  ------------
#   + drawing the loops
#
##########################################################################

Gplot.graphics<-function( mat, thresh=0, xlim=NULL, ylim=NULL, scale=1.0,
                          margin=0.2, main="", sub=""
                        ) {
# -------- Arguments -----------------------------------------------------
#
# mat        (matrix): initial matrix
# thresh     (scalar): matrix values are set to zero if they are < threshold
#                      (Pb: should be "abs( mat ) < 0")
# xlim, ylim (vector): x minimun and x maximum (same for y)
# scale  (scalar)    : scaling factor for the whole graphics.
# margin (scalar)    : margin value to add to all the graphics box sides.
# main  (char): title
# sub   (char): subtitle
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters (default)
#
# ------------------------------------------------------------------------

   n<-dim(mat)[1]

   # Replace NAs with 0s
   mat[is.na(mat)]<-0

   # Save a copy of mat
   mat.raw<-mat

   # Binary matrix
   mat<-matrix(as.numeric(mat>thresh),n,n)

   l.graphics <- list( xlim=xlim, ylim=ylim, scale=scale, margin=margin,
                       main=main, sub=sub, baserad=0, resolution=0 )

   l.network <- list( mode=NULL, loop.draw=FALSE, vertex.pos.mode=NULL,
                      coord= NULL, displayisolates=TRUE, use=NULL,
                      vertex.reso=0, edge.reso=0, arrow.reso=0)

   l.label   <- list( label=c(1:dim(mat)[1]), cex=1, col=1, pos=0,
                      useboxes=TRUE, box.margin=0.5, box.col=1,
                      box.bg="white", box.lty=NULL, box.lwd=par("lwd") )

   l.vertex <- list( cex=1, sides=8, col=2, border=1, lty=1 ,
                     label = l.label )

   l.edge  <- list( col=1, lty=1, lwd=0,
                    arrow.draw=TRUE, arrow.cex=1, loop.cex=1 )

   graph   <- list( mat=mat, thresh=thresh, mat.raw = mat.raw,
                    graphics=l.graphics, network=l.network,
                    vertex=l.vertex, edge=l.edge )

   return( graph )
}


Gplot.network <-function( graph,
                         # Bazard des arguments
                          mode=NULL, loop.draw=NULL,
                          vertex.pos.mode="default",
                          coord= NULL,
                          random.pos="circle",
                          displayisolates=NULL,
                          class=NULL,
                         ...)
{
# -------- Arguments -----------------------------------------------------
#
# graph (list): structure handling all graphics parameters
# mode  (char): kind of network to draw (not used).
# loop.draw       (bool): draw network loops (not implemented).
# vertex.pos.mode (char): vertex positionning mode (not used).
#                         "default" : simulated annealing
#                         "circle" : nodes are positionned on a circle
# coord         (vector): vertex position.
#                         If NULL these coordinates are computed
# random.pos              : "random" randomize nodes on a circle.
# displayisolates (bool): display isolate vertices.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------

  # Arguments mode and vertex.pos.mode are not implemented
  if( ! is.null(loop.draw) )
    graph$network$loop.draw = loop.draw
  if( ! is.null(coord) &  is.numeric(coord) )
    # Store given coordinates
    graph$network$coord = coord
  if( ! is.null(displayisolates) )
    graph$network$displayisolates = displayisolates

  Env <- get("Gplot.graph", envir=Gplot.env)

  if( is.null(graph$network$coord ) | is.character(coord) ) {
    if (is.null(graph$network$coord))
      vertex.pos.mode <- "default"
    if( is.character(coord) )
      vertex.pos.mode <- coord

    if( vertex.pos.mode == "default" ) {

      # Provide default settings
      n           <- dim(graph$mat)[1]
      niter       <- 500
      max.delta   <- n
      area        <- n^2
      cool.exp    <- 3
      repulse.rad <- area*n

      # Set initial positions (randomly) on the circle
      tempa<-(0:(n-1))/n
      if( random.pos == "random" ) {
        tempa <- sample( tempa )
      }
      x<-n/(2*pi)*sin(2*pi*tempa)
      y<-n/(2*pi)*cos(2*pi*tempa)

      layout<-.C( "vertex_coord_C",
                 as.integer(graph$mat),
                 as.double(n), as.integer(niter),
                 as.double(max.delta), as.double(area),
                 as.double(cool.exp), as.double(repulse.rad),
                 x=as.double(x), y=as.double(y)#,
                # PACKAGE="simone"
                 )

      graph$network$coord <- cbind(layout$x,layout$y)
    } else if ( vertex.pos.mode == "circle" ) {

        n           <- dim(graph$mat)[1]

        if ( is.null( class ) ) {
        # Set initial positions on the circle
        tempa<-(0:(n-1))/n
        x <- n/(2*pi)*sin(2*pi*tempa)
        y <- n/(2*pi)*cos(2*pi*tempa)
        graph$network$coord <- cbind( x, y)
      } else {
         f.class <- factor( class )
         class.dims <- rep(0, nlevels( f.class ))
         n.class <- length( class.dims )
         for (i in 1:nlevels( f.class ) ) {
           class.dims[i] <- length( which( f.class == levels(f.class)[i] ) )
         }
         IDInCircle <- rep( 0, n.class )
         IDInCircle[1] <- 0
         for ( i in 2:n.class) {
            IDInCircle[i] = IDInCircle[i-1] + class.dims[i-1]
         }
         xx <- vector()
         yy <- vector()
         cst <- n/(2*pi)
         for ( i in 1:n ) {
           # compute coord
           class.i = class[i]
           card.class = class.dims[ class.i ]
           x <- cst * sin( IDInCircle[ class.i ] * 2 * pi / n )
           y <- cst * cos( IDInCircle[ class.i ] * 2 * pi / n )
           IDInCircle[ class.i ] = IDInCircle[ class.i ] + 1
           xx <- append( xx, x)
           yy <- append( yy, y)
         }
         graph$network$coord <- cbind( xx, yy)
       }

    } else if ( vertex.pos.mode == "circles" ) {

      n           <- dim(graph$mat)[1]

      # Compute class size
      f.class <- factor( class )
      class.dims <- rep(0, nlevels( f.class ))
      for (i in 1:nlevels( f.class ) ) {
        class.dims[i] <- length( which( f.class == levels(f.class)[i] ) )
      }


      # Set initial positions on the circle
      n.class <- length( class.dims )
      ang <- 2*pi /  n.class
      max.radius <- max( class.dims ) / (2*pi)
      # Space for labels between two circles
      # Space estimate : 2 words of 12 characters
      # each character have 10 pt resolution
      # 4 * radius is assumed to be the physical plot size
      max.radius <- max.radius + 4*max.radius / 1000 * 10 * 12
      # graph$network$vertex.reso
      xx <- vector()
      yy <- vector()
      two.pi <- 2*pi
      one.over.two.pi <- 1/two.pi

      x.shift <- rep( 0, n.class )
      y.shift <- rep( 0, n.class )

      for ( i in 1:n.class ) {

        # shift
        x.shift[i] <-  max.radius * cos( (i-1) * ang + pi*0.5)
        y.shift[i] <-  max.radius * sin( (i-1) * ang + pi*0.5)
      }
      IDInClass <- rep( 0, n.class )
      for ( i in 1:n ) {
        # compute coord
        class.i = class[i]
        card.class = class.dims[ class.i ]
        x <- x.shift[ class.i ] + card.class*one.over.two.pi*
                  sin( 2*pi*( IDInClass[ class.i ] )/card.class)
        y <- y.shift[ class.i ] + card.class*one.over.two.pi*
                  cos( 2*pi*( IDInClass[ class.i ] )/card.class )
        IDInClass[ class.i ] = IDInClass[ class.i ] + 1
        xx <- append( xx, x)
        yy <- append( yy, y)
      }
      graph$network$coord <- cbind( xx, yy)

    }
  }

  # Remove isolated vertex (if displayisolates FALSE)
  use <- displayisolates | ( ! isolate.vertices( graph$mat ) )

  graph$network$use = use

  x <- graph$network$coord[,1]
  y <- graph$network$coord[,2]


  # Set limits for plotting region
  xlim = graph$graphics$xlim
  ylim = graph$graphics$ylim
  margin = graph$graphics$margin
  if(is.null(xlim)) {
    xmin <- min(x[use]); xmax <- max(x[use])
    xlim<-c( xmin - (xmax - xmin) * margin, xmax + (xmax - xmin) * margin)  # Save x, y limits
  }
  if(is.null(ylim)) {
    ymin <- min(y[use]); ymax <- max(y[use])
    ylim<-c( ymin - (ymax - ymin) * margin, ymax + (ymax - ymin) * margin )  # Save x, y limits
  }
  xrng<-diff(xlim)
  yrng<-diff(ylim)
  xctr<-(xlim[2]+xlim[1])/2                 # Get center of plotting region
  yctr<-(ylim[2]+ylim[1])/2

  # Force scale to be symmetric
  if(xrng<yrng)
    xlim<-c(xctr-yrng/2,xctr+yrng/2)
  else
    ylim<-c(yctr-xrng/2,yctr+xrng/2)

  graph$graphics$xlim = xlim
  graph$graphics$ylim = ylim

  # Extract "base radius"
  graph$baserad <- min(diff(xlim),diff(ylim))* graph$graphics$scale

  # Resolution
  graph$resolution <- min(diff(xlim),diff(ylim))* graph$graphics$scale / Env$reso

  # Sizes of reference in resolution unit (in pixels)
  graph$network$vertex.reso = 20
  graph$network$edge.reso   = 0.5
  graph$network$arrow.reso  = 40
  graph$network$arrow.angle  = 20*pi/180

  # Configure the graphic box
  plot( 0,0,
        xlim=graph$graphics$xlim,
        ylim=graph$graphics$ylim, type="n", xlab="",ylab="", asp=1, axes=FALSE,
        main= graph$graphics$main, sub=graph$graphics$sub,
        ...
       )
  return( graph )
}


Gplot.vertex <-function( graph,
                         cex, sides, col, border, lty, ...)
{
# -------- Arguments -----------------------------------------------------
#
# graph           (list): structure handling all graphics parameters
# cex    (scalar/vector): scaling factors.
# sides  (scalar/vector): number of node sides.
# col    (scalar/vector): vertex colors.
# border (scalar/vector): color of node (vertex) borders.
# lty    (scalar/vector): type of nodes borders.
#                         (Pb: would be better with the line width
#                         border parameter "lwd".
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------

  if( ! missing(cex) )
   graph$vertex$cex    = cex
  if( ! missing(sides) )
   graph$vertex$sides  = sides
  if( ! missing(col) )
   graph$vertex$col    = col
  if( ! missing(border) )
   graph$vertex$border = border
  if( ! missing(lty) )
   graph$vertex$lty    = lty

   n <- dim(graph$mat)[1]

   # Build vectors describing vertex
   v.cex    <- rep( graph$vertex$cex                       , length=n)
   v.radius <- graph$resolution * graph$network$vertex.reso * v.cex
   v.sides  <- rep( graph$vertex$sides                     , length=n)
   v.col    <- rep( graph$vertex$col                       , length=n)
   v.border <- rep( graph$vertex$border                    , length=n)
   v.lty    <- rep( graph$vertex$lty                       , length=n)

   # remove unused
   use = graph$network$use
   v.radius <-  v.radius[use]
   v.sides  <-  v.sides[use]
   v.col    <-  v.col[use]
   v.border <-  v.border[use]
   v.lty    <-  v.lty[use]

   x <- graph$network$coord[use,1]
   y <- graph$network$coord[use,2]
   n <- length(x)

 # Compute the coordinates
  coord<-vector()

  for(i in 1:n){
    ang <- (1:v.sides[i])/v.sides[i]*2*pi
    dx <- v.radius[i]*cos(ang)
    dy <- v.radius[i]*sin(ang)
    XY = rbind( cbind( x[i]+dx, y[i]+dy ), c(NA,NA) )
    coord<-rbind(coord, XY)
  }
  # Plot the vertices
  polygon(coord, col=v.col, border=v.border, lty=v.lty, ...)

  return( graph )
}

Gplot.vertex.pie <-function( graph,
                             cex=NULL, pie.coef=NULL, col=-1, border=NULL, lty=NULL, display=TRUE, ...)
{
# -------- Arguments -----------------------------------------------------
#
# graph                (list): structure handling all graphics parameters
# cex         (scalar/vector): scaling factors.
# pie.coef    (scalar/vector): Pie chat coefficients.
# col         (scalar/vector): vertex colors.
# border      (scalar/vector): color of node (vertex) borders.
# lty         (scalar/vector): type of nodes borders.
#                              (Pb: would be better with the line width
#                              border parameter "lwd".
# display            (scalar): Save the nodes characteristics and  if TRUE
#                              displays the pies.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------

  if( ! is.null(cex) )
    graph$vertex$cex <- cex

  draw.pie <- FALSE
  n.sector <- 1
  if( ! is.null( pie.coef ) ) {
    if (dim( pie.coef )[1] == dim(graph$mat)[1]) {
      n.sector <- dim( pie.coef )[2]
      draw.pie <- TRUE
    }
  }

  if ( length(col)  < n.sector ) {
    col <- sample( light.palette(n.sector, black=FALSE) )
  }

  if( ! is.null(border) )
   graph$vertex$border <- border
  if( ! is.null(lty) )
    graph$vertex$lty    <- lty

  graph$vertex$sides = -1

  n <- dim(graph$mat)[1]

  # Build vectors describing vertex
  scale.vertex <- graph$resolution * graph$network$vertex.reso
  v.cex    <- rep( graph$vertex$cex                       , length=n)
  v.radius <- scale.vertex * v.cex
  v.sides  <- rep( graph$vertex$sides                     , length=n)
  v.col    <- rep( graph$vertex$col                       , length=n)
  v.border <- rep( graph$vertex$border                    , length=n)
  v.lty    <- rep( graph$vertex$lty                       , length=n)

  # remove unused
  use = graph$network$use
  v.radius <-  v.radius[use]
  v.sides  <-  v.sides[use]
  v.col    <-  v.col[use]
  v.border <-  v.border[use]
  v.lty    <-  v.lty[use]

  # Display
  if ( display ) {
    x <- graph$network$coord[use,1]
    y <- graph$network$coord[use,2]
    n <- length(x)

    # Compute the coordinates
    if ( draw.pie ) {
      if (col[1] == -1) {
        # Generate new colors for sectors
        v.col <- sample( light.palette(n.sector, black=FALSE) )

      } else {
        # Use the colors of vertices
        v.col <- col
      }
      # Percent normalization
      sum.coef <- rowSums( pie.coef)
      cumul.theta <- 0
      for (i in 1:n){
        sector.theta <- pie.coef[i, ] * 2 * pi / sum.coef[i]
        sector.start <- rep(0, n.sector )

        if ( is.nan( sector.theta[1] ) ) {
          # Null coefficients - No Color
          draw.circle( graph, radius=v.radius[i], sides=-1,
                       center.x=x[i], center.y=y[i],
                       col="white",
                       border=v.border, lty=v.lty, ...  )
        } else if ( any( (sector.theta[] + 1.0e-7) > 2 * pi ) )  {
          # Only one class - One Color
          if( col[1] == -1 ) {
            # Sector colors
            draw.circle( graph, radius=v.radius[i], sides=-1,
                        center.x=x[i], center.y=y[i],
                        col=v.col[ which.max( pie.coef[i,] )],
                        border=v.border, lty=v.lty, ...  )
          } else {
            # Colors of vertices case
            draw.circle( graph, radius=v.radius[i], sides=-1,
                        center.x=x[i], center.y=y[i],
                        col=v.col[ which.max( pie.coef[i,] )],
                        border=v.border, lty=v.lty, ...  )
          }
        } else {
          for( j in 2:n.sector) {
            sector.start[j] <- sector.start[j-1] + sector.theta[j-1]
          }

          draw.sector( graph, radius=v.radius[i],
                      theta.start=sector.start,
                      theta=sector.theta,
                      center.x=x[i], center.y=y[i],
                      col=v.col, border=v.border, lty=v.lty, ...  )
        }
      }
    } else {
      draw.circle( graph, radius=v.radius, sides=-1, center.x=x, center.y=y ,
                  col=v.col, border=v.border, lty=v.lty, ...  )
    }
  }
  return( graph )
}



Gplot.vertex.label <-function( graph,
                               label, cex, col, pos,
                               useboxes, box.margin, box.col, box.bg,
                               box.lty, box.lwd,
                               ... )
{
# -------- Arguments -----------------------------------------------------
#
# graph           (list): structure handling all graphics parameters
# label         (vector): label titles (if NULL a default is provided)
# cex    (scalar/vector): scaling factors.
# col    (scalar/vector): label colors.
# pos           (scalar): label positionning mode
#                         + 0 labels are placed away from the graph
#                         + 1 labels are placed below the vertices
#                         + 2 labels are placed on the vertex left.
#                         + 3 labels are placed above the vertices
#                         + 4 labels are placed on the vertex right.
# useboxes        (scalar): frame (box) the labels
# box.margin      (scalar): margin between the label titles and their boxes
#                           (in character size unit).
# box.col  (scalar/vector): box colors.
# box.bg   (scalar/vector): box backgroung color.
# box.lty  (scalar/vector): boxe line type.
# box.lwd  (scalar/vector): boxe line width.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------
  if( ! missing( label ) )
    if( ! is.null( label ) )
      graph$vertex$label$label = label
  if( ! missing( cex ) )
   graph$vertex$label$cex   = cex
  if( ! missing( col ) )
   graph$vertex$label$col   = col
  if( ! missing( pos ) )
   graph$vertex$label$pos   = pos
  if( ! missing( useboxes ) )
   graph$vertex$label$useboxes   = useboxes
  if( ! missing( box.margin ) )
   graph$vertex$label$box.margin = box.margin
  if( ! missing( box.col ) )
   graph$vertex$label$box.col    = box.col
  if( ! missing( box.bg ) )
   graph$vertex$label$box.bg     = box.bg
  if( ! missing( box.lty ) )
   graph$vertex$label$box.lty    = box.lty
  if( ! missing( box.lwd ) )
   graph$vertex$label$box.lwd    = box.lwd

  # Plot vertex labels
  use  <- graph$network$use
  x <- graph$network$coord[use,1]
  y <- graph$network$coord[use,2]

  if((!all(graph$vertex$label$label==""))&(!all(use==FALSE))){

    # Label display mode
    if ( graph$vertex$label$pos == 0 ){

      # Labels are placed away from the graph
      xoff <- x - mean(x)
      yoff <- y - mean(y)
      roff <- sqrt(xoff^2+yoff^2)
      xhat <- xoff/roff
      yhat <- yoff/roff

    } else if (graph$vertex$label$pos<5) {

      # below (0,-1), left (-1,0), top (0,1) , right (1,0)
      xhat <- switch( graph$vertex$label$pos,  0,-1, 0, 1)
      yhat <- switch( graph$vertex$label$pos, -1, 0, 1, 0)

    } else {
      xhat <- 0
      yhat <- 0
    }

    # Get character size
    l.cex <- graph$vertex$label$cex
    char.len <- par()$cxy * l.cex

    # Get label width and height
    label <- graph$vertex$label$label[use]
    lw <- strwidth( label,cex=l.cex) / 2
    lh <- strheight(label,cex=l.cex) / 2
    b.size   = 1 + graph$vertex$label$box.margin

    scale.vertex <-  graph$resolution * graph$network$vertex.reso
    v.radius <- scale.vertex * rep( graph$vertex$cex,
                                            dim(graph$mat)[1] )

    # Draw boxes
    if( graph$vertex$label$useboxes ){
      rect(x - lw*b.size + xhat*(lw*(b.size+0.2) + v.radius),
           y - lh*b.size + yhat*(lh*(b.size+0.2) + v.radius),
           x + lw*b.size + xhat*(lw*(b.size+0.2) + v.radius),
           y + lh*b.size + yhat*(lh*(b.size+0.2) + v.radius),
           col    = graph$vertex$label$box.bg,
           border = graph$vertex$label$box.border,
           lty    = graph$vertex$label$box.lty,
           lwd    = graph$vertex$label$box.lwd)
    }

    # Draw labels
    text(x + xhat * ( lw*(b.size+0.2) + v.radius ),
         y + yhat * ( lh*(b.size+0.2) + v.radius ),
         label, cex=l.cex, col=graph$vertex$label.col, offset=0, ...)

  }
  return ( graph )
}

Gplot.radial.vertex.label <-function( graph,
                               class.dims=NULL,
                               label, cex, col, pos,
                               useboxes, box.margin, box.col, box.bg,
                               box.lty, box.lwd,
                               ... )
{
# -------- Arguments -----------------------------------------------------
#
# graph           (list): structure handling all graphics parameters
# class.dims    (vector): number of nodes by class/circle
# label         (vector): label titles (if NULL a default is provided)
# cex    (scalar/vector): scaling factors.
# col    (scalar/vector): label colors.
# pos           (scalar): label positionning mode
#                         + 0 labels are placed away from the graph
#                         + 1 labels are placed below the vertices
#                         + 2 labels are placed on the vertex left.
#                         + 3 labels are placed above the vertices
#                         + 4 labels are placed on the vertex right.
# useboxes        (scalar): frame (box) the labels
# box.margin      (scalar): margin between the label titles and their boxes
#                           (in character size unit).
# box.col  (scalar/vector): box colors.
# box.bg   (scalar/vector): box backgroung color.
# box.lty  (scalar/vector): boxe line type.
# box.lwd  (scalar/vector): boxe line width.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------
  if( ! missing( label ) )
    if( ! is.null( label ) )
      graph$vertex$label$label = label
  if( ! missing( cex ) )
   graph$vertex$label$cex   = cex
  if( ! missing( col ) )
   graph$vertex$label$col   = col
  if( ! missing( pos ) )
   graph$vertex$label$pos   = pos
  if( ! missing( useboxes ) )
   graph$vertex$label$useboxes   = useboxes
  if( ! missing( box.margin ) )
   graph$vertex$label$box.margin = box.margin
  if( ! missing( box.col ) )
   graph$vertex$label$box.col    = box.col
  if( ! missing( box.bg ) )
   graph$vertex$label$box.bg     = box.bg
  if( ! missing( box.lty ) )
   graph$vertex$label$box.lty    = box.lty
  if( ! missing( box.lwd ) )
   graph$vertex$label$box.lwd    = box.lwd

  # Plot vertex labels
  use  <- graph$network$use
  x <- graph$network$coord[use,1]
  y <- graph$network$coord[use,2]

  if( is.null(class.dims) )
    class.dims <- length(x)

  if((!all(graph$vertex$label$label==""))&(!all(use==FALSE))){

     # For different circles/classes
    index.end <- 0
    for( i.class in 1:length(class.dims) ) {
      index.start <-  index.end + 1
      index.end <- class.dims[i.class]
      index.size <- index.end - index.start + 1

     # TODO : label placing with circles mode
      xctr = sum( x[ index.start : index.end] ) / (index.size)
      yctr = sum( y[ index.start : index.end] ) / (index.size)

    # Label display mode
      if ( graph$vertex$label$pos == 0 ){


      # Labels are placed away from the graph
        xoff <- x[ index.start : index.end] - xctr
        yoff <- y[ index.start : index.end] - yctr
        roff <- sqrt(xoff^2+yoff^2)
        xhat <- xoff/roff
        yhat <- yoff/roff
        ang <- atan( yoff / xoff )
        ang <- ang[] + 0*(xoff[] < 0) * pi
      } else if (graph$vertex$label$pos<5) {

      # below (0,-1), left (-1,0), top (0,1) , right (1,0)
        xhat <- switch( graph$vertex$label$pos,  0,-1, 0, 1)
        yhat <- switch( graph$vertex$label$pos, -1, 0, 1, 0)

      } else {
        xhat <- 0
        yhat <- 0
      }

    # Get character size
      l.cex <- graph$vertex$label$cex
      char.len <- par()$cxy * l.cex

    # Get label width and height
      label <- graph$vertex$label$label[use]
      label <- label[ index.start : index.end]
    # ??? unused
      lw <- strwidth( label,cex=l.cex) / 2
      lh <- strheight(label,cex=l.cex) / 2
      b.size   = 1 + graph$vertex$label$box.margin

      scale.vertex <-  graph$resolution * graph$network$vertex.reso
      v.radius <- scale.vertex * rep( graph$vertex$cex,
                                            dim(graph$mat)[1] )

    # Shift the label if loops are drawn
      l.shift <- ( graph$network$loop.draw ) * 2 *
        graph$network$vertex.reso * graph$resolution

    # Draw boxes
      if( graph$vertex$label$useboxes ){
        rect(x - lw*b.size + xhat*(lw*(b.size+0.2) + v.radius),
             y - lh*b.size + yhat*(lh*(b.size+0.2) + v.radius),
             x + lw*b.size + xhat*(lw*(b.size+0.2) + v.radius),
             y + lh*b.size + yhat*(lh*(b.size+0.2) + v.radius),
             col    = graph$vertex$label$box.bg,
             border = graph$vertex$label$box.border,
             lty    = graph$vertex$label$box.lty,
             lwd    = graph$vertex$label$box.lwd)
      }

      # Draw labels
      b.size <- 0.8
      for ( i in 1:index.size ) {
        par(srt= ang[i]*180/pi)

      # Text alignment
        just <- 4
        if( xoff[i] < 0) just <- 2

        text(x[index.start + i] + xhat[i] * ( l.shift + v.radius[i] ),
             y[index.start + i] + yhat[i] * ( l.shift + v.radius[i] ),
             label[i], cex=l.cex[index.start+i],
             col=graph$vertex$label.col[index.start+i],
             pos=just, offset=0, ...)
      }
    }
    par(srt=0)
  }
  return ( graph )
}

Gplot.edge <-function( graph,
                       col=1, lty=1, lwd=0,
                       arrow.draw=TRUE, arrow.cex=2, curved=FALSE,
                       loop.draw=FALSE, loop.arrow = TRUE, loop.cex=1,
                       ... )
{
# -------- Arguments -----------------------------------------------------
#
# graph              (list): structure handling all graphics parameters
# col  (scalar/vector/matrix): edge colors.
# lty  (scalar/vector/matrix): edge line types.
# lwd  (scalar/vector/matrix): edge line widths.
# arrow.draw      (bool)   : draw arrows.
# arrow.cex (scalar/vector): arrow head scaling factors.
# curved       (bool): draw curved edges.
# loop.draw          (bool): draw loops.
# loop.arrow         (bool): draw loop arrows.
# loop.cex  (scalar/vector): loop scaling factors.
#
# -------- Return value --------------------------------------------------
#
# graph (list): structure handling all graphics parameters
#
# ------------------------------------------------------------------------

  graph$edge$col = col
  graph$edge$lty = lty
  graph$edge$lwd = lwd
  graph$edge$arrow.draw = arrow.draw
  graph$edge$arrow.cex = arrow.cex
  graph$edge$loop.cex  = loop.cex

  n <- dim( graph$mat )[1]

   # Build vectors describing edges
   # Each edge is a polygon
   px0<-vector()
   py0<-vector()
   px1<-vector()
   py1<-vector()

   e.lwd<-vector()  # Create edge attribute vectors
   e.type<-vector()
   e.col<-vector()
   e.hoff<-vector() # Offset radii for heads
   e.toff<-vector() # Offset radii for tails
   e.diag<-vector() # Indicator for self-ties
   e.curv<-vector() # Curvature value

   l.px0 <- vector()
   l.py0 <- vector()
   l.e.lwd  <- vector()
   l.e.type <- vector()
   l.e.col  <- vector()
   l.e.off  <- vector()
   l.e.rad  <- vector()

   # Coerce edge.col/edge.lty to array form
   if(!is.array(graph$edge$col))
     col <- array(graph$edge$col, dim=dim(graph$mat))
   else
     col <- graph$edge$col
   if(!is.array(graph$edge$lty))
     lty<-array(graph$edge$lty, dim=dim(graph$mat))
   else
     lty = graph$edge$lty
   if(!is.array( graph$edge$lwd)){
     if( graph$edge$lwd>0 )
       lwd<-array( graph$edge$lwd * graph$mat.raw, dim=dim(graph$mat))
     else
       lwd<-array(1, dim=dim(graph$mat))
   }

   # Used for curved edges
   dist   <- as.matrix( dist(graph$network$coord  ))
   tl     <- graph$mat.raw*dist    # Rescaled edge lengths
   tl.max <- max(tl)               # Maximum edge length

   scale.vertex <-  graph$resolution *  graph$network$vertex.reso
   scale.edge   <-  graph$resolution *  graph$network$edge.reso
   scale.arrow   <- graph$resolution *  graph$network$arrow.reso

   v.radius <- scale.vertex * rep( graph$vertex$cex, dim(graph$mat)[1] )

   # Select edges between vertices
   # -----------------------------
   x <- graph$network$coord[,1]
   y <- graph$network$coord[,2]

   for(i in (1:n)[graph$network$use]) {
     for(j in (1:n)[graph$network$use]) {
       if( graph$mat[i,j] ){        # Edge exists
         px0 <- c(px0,(x[i]))  # Store endpoint coordinates
         py0 <- c(py0,(y[i]))
         px1 <- c(px1,(x[j]))
         py1 <- c(py1,(y[j]))
         e.toff <-c ( e.toff, v.radius[i] ) # Store endpoint offsets
         e.hoff <-c ( e.hoff, v.radius[j] )
         e.col  <-c ( e.col , col[i,j])     # Store other edge attributes
         e.type <-c ( e.type, lty[i,j])
         e.lwd  <-c ( e.lwd , lwd[i,j])
         e.diag <-c ( e.diag, i==j)         # Set to true if diagonal
         if(curved){   # Curvature on interpoint distances
           e.curv      <- c( e.curv, 2 )
         }
         if( i == j) {
#           l.e.rad<-c( l.e.rad, scale.edge*lwd[i,i] * loop.cex )
           l.e.rad<-c( l.e.rad, v.radius * loop.cex )
         }
       }
     }
   }

   # Store loops
   # ------------

   if( (length(px0)>0) & loop.draw ){
     l.px0 <- px0[e.diag]
     l.py0 <- py0[e.diag]
     l.e.lwd  <- e.lwd[e.diag]
     l.e.type <- e.type[e.diag]
     l.e.col  <- e.col[e.diag]
     l.e.off  <- e.toff[e.diag]
   }

   # Remove loops
   # ------------
   if(length(px0)>0){
     px0 <- px0[!e.diag]
     py0 <- py0[!e.diag]
     px1 <- px1[!e.diag]
     py1 <- py1[!e.diag]
     e.lwd  <- e.lwd[!e.diag]
     e.type <- e.type[!e.diag]
     e.col  <- e.col[!e.diag]
     e.hoff <- e.hoff[!e.diag]
     e.toff <- e.toff[!e.diag]
     e.curv <- e.curv[!e.diag]
   }

   # Draw edges
   if(length(px0)>0)
     if( !curved ) {
       draw.edges(
            as.vector(px0),as.vector(py0),as.vector(px1),as.vector(py1),
            width=e.lwd * scale.edge,
            col=e.col, lty=e.type,
            o.head=e.hoff, o.tail=e.toff,
            arrow=arrow.draw,
            a.len=pmax( scale.edge * e.lwd * arrow.cex,
                        scale.arrow ),
            a.angle = graph$network$arrow.angle,
            ...)
     } else {
       # 0->1 vector angle
       phi <-  atan( (py1 - py0) /  (px1 - px0) )
       phi[ which( is.nan( phi ) ) ] <- pi/2

       #
       norm  <- sqrt( (py1 - py0)^2 + (px1 - px0)^2)
       radius <- norm * e.curv

       alpha  <- asin( 0.5 / e.curv )
       pi.rot <- pi * ((px1 - px0 ) < 0.0 )

       xc <- px0 - radius * cos( (pi/2 + alpha + phi) + pi.rot )
       yc <- py0 - radius * sin( (pi/2 + alpha + phi) + pi.rot )

       theta0 <-   pi*0.5 + phi + alpha - atan2( e.toff, radius ) + pi.rot
       theta <-  - 2 * alpha + atan2( e.toff, radius ) +  atan2(e.hoff, radius )

       # cat( "px0 = ", px0, "\n")
       # cat( "py0 = ", py0, "\n")
       # cat( "px1 = ", px1, "\n")
       # cat( "py1 = ", py1, "\n")
       # cat( "xc = ", xc, "\n")
       # cat( "yc = ", yc, "\n")
       # cat( "radius = ", radius, "\n")
       # cat( "theta0 = ", theta0, "\n")
       # cat( "theta = ", theta, "\n")
       # cat( "alpha = ", alpha, "\n")
       # cat( "a.len = ", phi, "\n")
       # cat( "graph$baserad = ", graph$baserad, "\n")
       # cat( "arrow.cex = ", arrow.cex, "\n")

       draw.arc( graph,
                 radius = radius,
                 theta.start = theta0, theta=theta,
                 xc = xc, yc = yc,
                 arrow=arrow.draw,
                 a.len=pmax( scale.edge*e.lwd*arrow.cex,scale.arrow),
                 a.angle = graph$network$arrow.angle,
                 col=e.col, width=e.lwd*scale.edge, lty=e.type, ...
       )

     }
   if( (length(l.px0)>0) & loop.draw ){

     draw.loops( graph,
                l.px0, l.py0,
                width   = l.e.lwd * scale.edge,
                col      = l.e.col,
                lty      = l.e.type,
                # loops coef. "cex" must be 2.5 times larger than the node radii
                r.beg = l.e.off, r.end =0,
                cex = pmax( l.e.rad, 2.5*l.e.off),
                theta.beg = 0, theta.end = 0, phi = 0,
                arrows       = loop.arrow,
                # Loop arrows are thiner than edge ones ( twice times)
                a.len       = 2 * pmax( l.e.lwd * scale.edge*arrow.cex, scale.arrow) ,
                a.angle     = -1,
                xctr = mean(x[graph$network$use]),
                yctr = mean(y[graph$network$use])
                )
   }
   return( graph )
}

isolate.vertices<-function( mat. ) {
# -------- Arguments -----------------------------------------------------
#
# mat. (matrix): binary (0/1) square matrix
#
# -------- Return value --------------------------------------------------
#
# isolate (vector): isolated (if TRUE) vertex vector.
#
# ------------------------------------------------------------------------

  n <- dim(mat.)[1];
  mat <- mat.
  isolate <- vector()
  if ( n > 1 ){
    # Set to zero NA and diagonal terms
    for(i in 1:n){
      mat[i,i] = 0
      mat[i,1:n] ==  as.numeric( ! is.na(mat[i,1:n]) )
    }
    for(i in 1:n) {
      isolate = c( isolate, all(( mat[i,] == 0 )) & all(( mat[,i] == 0 )) )
    }
  }
  return( isolate )
}

draw.edges<-function( x0, y0, x1, y1,
                      width=0.01, col=1, lty=1,
                      o.head=0, o.tail=0,
                      arrow=TRUE, a.len=0.4, a.angle=0.2,
                      ... )
{
# -------- Arguments -----------------------------------------------------
#
# x0, y0       (vector): start coordinates of edges to draw.
# x1, y1       (vector): end coordinates of edges to draw.
# width (scalar/vector): edge line widths.
# col   (scalar/vector): edge colors.
# lty   (scalar/vector): edge line types.
# o.head (scalar/vector): offset (vertex size shift) at the start points.
# o.tail (scalar/vector): offset (vertex size shift) at the end points.
# arrow           (bool): arrows are drawn.
# a.len   (scalar/vector): arrow head lengths.
# a.angle (scalar/vector): arrow head angle (in degree).
#
# -------- Return value --------------------------------------------------
#
# No value
#
# ------------------------------------------------------------------------

  if(length(x0)==0)   #Leave if there's nothing to do
    return;

  n<-length(x0)

  # Transform scalars into vectors
  width <- rep(width,length=n)
  col   <- rep(col,length=n)
  lty   <- rep(lty,length=n)

  # Offsets
  o.head  <- rep(o.head,length=n)
  o.tail  <- rep(o.tail,length=n)

  # Arrow parameters
  a.angle <- rep(a.angle,length=n)/360*2*pi
  a.len   <- rep(a.len,length=n)

  # Debug point
  # cat("xy  :",x0, y0, x1,y1, "\n")
  # cat("width :", width, "\n")
  # cat("col :", col, "\n")
  # cat("lty :", lty, "\n")
  # cat("o.tail :",o.tail, "\n")
  # cat("o.head :", o.head , "\n")
  # cat("a.angle :", a.angle, "\n")
  # cat("a.len   :", a.len, "\n")

  # Computes edges/arrows coordinates
  coord<-vector()
  XY <- vector()

  # Needed for curved edges
  sqrt2 <- sqrt(2)/2

  for(i in 1:n) {

    # Edge length
    slen<-sqrt((x0[i]-x1[i])^2+(y0[i]-y1[i])^2)

    # Arrow longuer than edge
    a.len[i] <-  min( 0.9*(slen - (o.head[i] + o.tail[i])), a.len[i])


    # If the node distance is greater than the arrow size
    if (arrow & a.len[i] > 0 ){
      #  With Arrows
      a.sin = sin( a.angle[i] )
      a.cos = cos( a.angle[i] )
      XY<-rbind(
                  c( - width[i]/2      , o.tail[i]),
                  c( - width[i]/2      , slen - 0.5*a.len[i] - o.head[i]),
                  c( - a.len[i] * a.sin, slen - a.len[i]*a.cos - o.head[i]),
                  c(   0               , slen - o.head[i]),
                  c(   a.len[i] * a.sin, slen - a.len[i]*a.cos - o.head[i]),
                  c(   width[i]/2      , slen - 0.5*a.len[i] - o.head[i]),
                  c(   width[i]/2      , o.tail[i] ),
                  c(   NA              , NA)
                )
    } else {

      #  Without Arrows
      XY<-rbind(
                  c( - width[i]/2, o.tail[i]       ),
                  c( - width[i]/2, slen - o.head[i]),
                  c(   width[i]/2, slen - o.head[i]),
                  c(   width[i]/2, o.tail[i]       ),
                  c(   NA,      NA)
                )
    }

    # Rotate
    theta <- atan2(y1[i]-y0[i],x1[i]-x0[i])-pi/2
    rmat  <- rbind(c(cos(theta),sin(theta)),c(-sin(theta),cos(theta)))
    XY    <- XY %*% rmat
    # Translate
    XY[,1] <- XY[,1]+x0[i]
    XY[,2] <- XY[,2]+y0[i]

    coord<-rbind( coord, XY)
  }

  #coord<-coord[-NROW(coord),]

  #Draw polygons
  polygon(coord,col=col,border=col,lty=lty,...)
}

draw.circle <- function ( graph, radius, sides, center.x, center.y ,
                           col, border, lty, ...  ) {
# -------- Arguments -----------------------------------------------------
# graph                  : current graph
# radius (scalar/vector) : circle radii.
# sides (scalar/vector)  : number of sides to draw a regular polygon
#                          if -1 the number of sides is done to approximate a circle
# center.x[.y] (scalar/vector) : coordinates of the circle center .
# col, border, lty   (scalar/vector) : color, color border, line type
  compute.side.nbr <- TRUE
  if( sides[1] != -1 )
    compute.side.nbr <- FALSE

  if ( compute.side.nbr ) {
    xy.size <- min( diff( graph$graphics$xlim ),  diff( graph$graphics$xlim ) )
    sides <- 2 * pi * 500 * radius / xy.size
  }

  coord<-vector()
  n <- length( radius )
  for(i in 1:n){
    ang <- (1:sides[i])/sides[i]*2*pi
    dx  <- radius[i]*cos(ang)
    dy  <- radius[i]*sin(ang)
    XY  <- rbind( cbind( center.x[i]+dx, center.y[i]+dy ), c(NA,NA) )
    coord <- rbind(coord, XY)
  }
  # Plot the vertices
  polygon(coord, col=col, border=border, lty=lty, ...)
}

draw.sector <- function ( graph, radius, theta.start, theta, center.x, center.y ,
                          col, border, lty, ...  ) {
# -------- Arguments -----------------------------------------------------
# graph                  : current graph
# radius (scalar/vector) : circle radii.
# theta.start            : initial angle to start the sector
# theta                  : sector angle
# center.x[.y] (scalar/vector)       : coordinates of the circle center .
# col, border, lty   (scalar/vector) : color, color border, line type

  n <- length( theta.start )
  if( length( radius ) != n ) {
    radius <- rep( radius, n )
    center.x <- rep( center.x[1], n)
    center.y <- rep( center.y[1], n)
  }
  xy.size <- min( diff( graph$graphics$xlim ),  diff( graph$graphics$xlim ) )
  sides <- round( theta * 500 * radius / xy.size )

  coord<-vector()

  for(i in 1:n){
    ang <- (0:sides[i])/sides[i] * theta[i] + theta.start[i]
    dx  <- radius[i]*cos(ang)
    dy  <- radius[i]*sin(ang)
    XY  <- rbind( cbind( center.x[i]+dx, center.y[i]+dy ),
                 c(center.x[i], center.y[i]),
                 c(NA,NA) )
    polygon(XY, col=col[i], border=border, lty=lty, ...)

  }
}

draw.arc <- function ( graph, radius = radius,
                       theta.start = theta.start, theta = theta,
                       xc = xc, yc = yc ,
                       arrow = arrow, a.len = a.len, a.angle = a.angle,
                       col = col, width = width, lty = lty,...  ) {
# -------- Arguments -----------------------------------------------------
# graph                  : current graph
# radius (scalar/vector) : circle radii.
# theta.start            : initial angle to start the circle
# theta                  : arc angle
# xc, yc (scalar/vector) : coordinates of the circle center .
# arrow                  : draw arrow (TRUE/FALSE)
# a.len, a.angle,        : arrow shape
# col, border, lty   (scalar/vector) : color, color border, line type


  # Shape factor for the arrow
  a.shape.factor <- 0.5
  xy.size <- min( diff( graph$graphics$xlim ),  diff( graph$graphics$xlim ) )

  # Number of sides
  a.phi <- 0
  if( arrow ) {
    # remove a.len /2
    a.phi <- sign(theta)*asin( .666 * a.len[]  / radius[] )
    a.sin <- sin( a.angle[] )
  }

  # Remove the arrow
  theta.var <- (theta - a.phi)
  sides <-abs(  round( theta.var * 30 * radius / xy.size ) )

  coord<-vector()
  n <- length( radius )
  for( i in 1:n ){
    ang <- (0:sides[i])/sides[i] * theta.var[i] + theta.start[i]
    cos.ang <- cos(ang)
    sin.ang <- sin(ang)
    dx  <- (radius[i]+0.5*width[i])*cos.ang
    dy  <- (radius[i]+0.5*width[i])*sin.ang
    # cat( "dx[1] = ", dx[1]+xc[i], "\n")
    # cat( "dy[1] = ", dy[1]+yc[i], "\n")
    # cat( "radius[i] = ", radius[i], "\n")

    if ( arrow ) {
      r <- radius[i] + 0.5*width[i] + a.sin*a.len[i]
      th <-  ang[sides[i]+1] - a.shape.factor *a.phi[i]
      dx1  <- r * cos( th )
      dy1  <- r * sin( th )
      dx2 <- radius[i] * cos( theta[i] + theta.start[i] )
      dy2 <- radius[i] * sin( theta[i] + theta.start[i] )
      r <- radius[i] - 0.5*width[i] - a.sin*a.len[i]
      dx3  <- r * cos( th )
      dy3  <- r * sin( th )
    }
    cos.ang <- rev( cos.ang )
    sin.ang <- rev( sin.ang )
    dx.rev  <- (radius[i]-0.5*width[i])*cos.ang
    dy.rev  <- (radius[i]-0.5*width[i])*sin.ang
    # cat( "dx.rev[1] = ", dx.rev[1] + xc[i], "\n")
    # cat( "dy.rev[1] = ", dy.rev[1] + yc[i], "\n")
    # cat( "ang[1]    = ", rev(ang)[1], "\n")

    coord <- rbind( coord, cbind( xc[i]+dx, yc[i]+dy ) )

    if( arrow ) {
      coord <- rbind( coord, cbind( xc[i] + c( dx1, dx2, dx3 ),
                                    yc[i] + c( dy1, dy2, dy3 ) ) )
    }

    coord <- rbind( coord, cbind( xc[i]+dx.rev, yc[i]+dy.rev ),
                    c(NA,NA) )
  }
  # Plot the vertices
  polygon(coord, col=col, border=col, lty=lty, ...)
}


draw.loops <- function ( graph,
                            xc, yc,
                            width    = 0.01, col = 1, lty = 1,
                            r.beg = 0, r.end=0, cex = 1,
                            theta.beg=0, theta.end=0, phi = 0,
                            arrows = FALSE, a.len = 0.1, a.angle = 10,
                            xctr=0, yctr=0, ... ) {

  # Compute "theta" versus "rho" and the coefficient "stroph.coef" of the function
  RhoToTheta <- function( stroph.coef, rho) {
    var <- 0.25/stroph.coef[] *
                 ( rho[] + sqrt( rho[]^2 + 8*stroph.coef[]^2) )
    if ( all( abs(var) > 1) ) {
        cat("Gplot error: bad loop caracteristics \n")
      }
    theta <- acos( var  )
  }

  # Number of loops
  np <- length(xc)

  if (length( a.angle ) == 1) a.angle <- rep( -1, np )

  # Test if arrows are too large
  arrow <- arrows & (a.len[] < 2 * cex)
  stroph.coef <- cex[]

  # Shape factor for the arrow
  a.shape.factor <- 0.6

  # Remove width nul values
  indexes <- which( width[] < 1/400 )
  width[indexes] <- 1/400

  # Remove a.len nul values
  indexes <- which( a.len[] < 1/400 )
  a.len[indexes] <- stroph.coef[indexes] / 4

  a.angle[indexes] <- 20*pi/180

  # Compute the arrow angle
  indexes <- which( a.angle[] == -1 )
  a.angle[indexes] <-
               ( atan( width[indexes] * 0.5 / (a.shape.factor*a.len[indexes]) ))

  # Too small arrow angle
  ang.min <-  20*pi/180
  indexes <- which( a.angle[] < ang.min )
  a.angle[indexes] <- ang.min

  # Too large arrows
  indexes <- which( a.len[] > 0.4 * stroph.coef[] )
  a.len[indexes] <-  0.4 * stroph.coef[indexes]
  a.angle[indexes] <- atan( width[indexes] * 0.5 / (a.shape.factor*a.len[indexes]))

  xy.size <- min( diff( graph$graphics$xlim ),  diff( graph$graphics$xlim ) )

  r.end <- r.beg + arrow[] * a.shape.factor * a.len[]


  # Loop rotation
  theta0.rot <- 0

  index <- which( abs(xc[] - xctr) > .Machine$double.eps)
  theta0.rot[ index ] <- - atan( (yc[index] - yctr) / (xc[index] - xctr) ) +
                           pi * ( (xc[index] - xctr) < 0 )

  # Starting and ending angle of the middle loop
  theta.beg.stroph <- - RhoToTheta( stroph.coef[], r.beg[] )
  theta.end.stroph <- RhoToTheta( stroph.coef[], r.end[] )

  # Starting Angles of the internal(min) and external(max) loops (min)
  x.stroph.min <- r.beg[] * cos( theta.beg.stroph ) - 0.5 * 0.5*width[] * cos(pi- theta.beg.stroph)
  y.stroph.min <- r.beg[] * sin( theta.beg.stroph ) - 0.5 * 0.5*width[] * sin(pi- theta.beg.stroph)
  x.stroph.max <- r.beg[] * cos( theta.beg.stroph ) + 0.5 * 0.5*width[] * cos(pi- theta.beg.stroph)
  y.stroph.max <- r.beg[] * sin( theta.beg.stroph ) + 0.5 * 0.5*width[] * sin(pi- theta.beg.stroph)

  rho.stroph.min <- sqrt(  (x.stroph.min - 0.5*width[])^2 + y.stroph.min^2 )
  rho.stroph.max <- sqrt(  (x.stroph.max + 0.5*width[])^2 + y.stroph.max^2 )
  theta.beg.stroph.min <- - RhoToTheta( stroph.coef[]-0.5*width[], rho.stroph.min   )
  theta.beg.stroph.max <- - RhoToTheta( stroph.coef[]+1.0*width[], rho.stroph.max )

  # End Angles of the internal(min) and external(max) loops (min)
  x.stroph.min <- r.end[] * cos( theta.end.stroph ) + 0.5 * 0.5*width[] * sin(pi- theta.end.stroph)
  y.stroph.min <- r.end[] * sin( theta.end.stroph ) + 0.5 * 0.5*width[] * cos(pi- theta.end.stroph)
  x.stroph.max <- r.end[] * cos( theta.end.stroph ) - 0.5 * 0.5*width[] * sin(pi- theta.end.stroph)
  y.stroph.max <- r.end[] * sin( theta.end.stroph ) - 0.5 * 0.5*width[] * cos(pi- theta.end.stroph)

  rho.stroph.min <- sqrt(  (x.stroph.min - 0.5*width[])^2 + y.stroph.min^2 )
  rho.stroph.max <- sqrt(  (x.stroph.max+0.5*width[])^2 + y.stroph.max^2 )

  theta.end.stroph.min <-   RhoToTheta( stroph.coef[]-0.5*width[], rho.stroph.min   )
  theta.end.stroph.max <-   RhoToTheta( stroph.coef[]+1.0*width[], rho.stroph.max)
  theta.stroph.min <- theta.end.stroph.min  - theta.beg.stroph.min
  theta.stroph.max <- theta.end.stroph.max  - theta.beg.stroph.max

  sides <-abs(  round( theta.stroph.max[] * 400 * stroph.coef[] / xy.size ) )

### brute (julien)

  if (length(sides) < np) {
    sides <- rep(sides,np)
    arrow <- rep(arrow,np)
    width <- rep(width,np)
  }

  coord<-vector()
  a.du <- vector(); a.dv <- vector();
  a.dx <- vector(); a.dy <- vector();

  n <- length( stroph.coef )
  for( i in 1:np ){
    ang <- (0:sides[i])/sides[i] * ( theta.stroph.min[i]  ) +
       ( theta.beg.stroph.min[i])
    cos.ang <- cos(ang)
    sin.ang <- sin(ang)
    cos.rot <- cos(theta0.rot[i])
    sin.rot <- sin(theta0.rot[i])

    # Debug
    # draw.rot.poly( xc[i], yc[i], x.stroph.min[i], y.stroph.min[i],
    #                c(0,0), cos.rot, sin.rot, color="red" )
    # draw.rot.poly( xc[ i], yc[i], x.stroph.max[i], y.stroph.max[i],
    #                c(0,0),  cos.rot, sin.rot, color="blue" )

    du  <- ((stroph.coef[i]-0.5*width[i])*(2*cos.ang - 1.0/cos.ang) ) * cos.ang
    dv  <- ((stroph.coef[i]-0.5*width[i])*(2*cos.ang - 1.0/cos.ang) ) * sin.ang
    # Debug
    # draw.rot.poly( xc[i], yc[i], du[1]+0.5*width[i], dv[1],
    #               c(0,0),  cos.rot, sin.rot )
    dx <-   (du + 0.5*width[i]) *  cos.rot + (dv)* sin.rot
    dy <- - (du + 0.5*width[i]) *  sin.rot + (dv)* cos.rot

    if ( arrow[i] ) {

      # Arrow slope
      x.end <- r.end[i] * cos(theta.end.stroph[i])
      y.end <- r.end[i] * sin(theta.end.stroph[i])
      x.beg <- r.beg[i] * cos(theta.end.stroph[i])
      y.beg <- r.beg[i] * sin(theta.end.stroph[i])
      beta <- atan( (y.end - y.beg) / (x.end - x.beg) )
      # phasis -0.05
      a.du[1] <- x.beg + a.len[i] * cos( beta - 1.0*a.angle[i] - 0.05)
      a.dv[1] <- y.beg + a.len[i] * sin( beta - 1.0*a.angle[i] -0.05)

      a.du[2] <- x.beg
      a.dv[2] <- y.beg

      a.du[3] <- x.beg + a.len[i] * cos( beta  + 1.0*a.angle[i] -0.05)
      a.dv[3] <- y.beg + a.len[i] * sin( beta  + 1.0*a.angle[i] -0.05)

      a.dx <-   a.du[] *  cos.rot + (a.dv[]) * sin.rot
      a.dy <- - a.du[] *  sin.rot + (a.dv[]) * cos.rot

    }

    cos.rot <- rev( cos.rot )
    sin.rot <- rev( sin.rot )

    ang <- rev( (0:sides[i])/sides[i] * ( theta.stroph.max[i] ) +
               theta.beg.stroph.max[i] )

    cos.ang <- cos(ang)
    sin.ang <- sin(ang)
    tan.ang <- tan(ang)

    du  <- ((stroph.coef[i]+1.0*width[i]) * (2*cos.ang - 1.0/cos.ang)) * cos.ang
    dv  <- ((stroph.coef[i]+1.0*width[i]) * (2*cos.ang - 1.0/cos.ang)) * sin.ang
    dx.rev <-   (du - 0.5*width[i]) *  cos.rot + (dv)* sin.rot
    dy.rev <- - (du - 0.5*width[i]) *  sin.rot + (dv)* cos.rot

    coord <- rbind( coord, cbind( xc[i]+dx, yc[i]+dy ) )

    if( arrow[i] ) {
      coord <- rbind( coord, cbind( xc[i] + a.dx[],
                                    yc[i] + a.dy[] ) )
    }

    coord <- rbind( coord, cbind( xc[i]+dx.rev, yc[i]+dy.rev ),
                    c(NA,NA) )
  }
  # Plot the vertices
  polygon(coord, col=col, border=col, lty=lty, ...)
}

# Debug utility
draw.rot.poly <- function(xc,yc,x, y, shift=c(0,0), cos.rot, sin.rot, color="black"){

      xx <- xc + (x[] + shift[1]) *  cos.rot + (y[])* sin.rot
      yy <- yc - (x[] + shift[1]) *  sin.rot + (y[])* cos.rot
      polygon( cbind( xx, yy), col=color )
      # stamp the first point
      points(xx[1], yy[1],col=color)
    }


light.palette <- function( n, black=TRUE, format="rgb" ){

  R1 <- vector()
  G1 <- vector()
  B1 <- vector()

  if ( ! black ) n = n + 1

#  Old values
#  R1.ref <- c( 0    , 0.08 , 0.15 ,
#          0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    ,
#          0    , 0.376, 0.6  , 0.8  , 0.9  , 0.95 , 1    , 1    , 1    ,
#          1    , 1    , 1    , 1
#        )
#  G1.ref <- c( 0    , 0    , 0    ,
#          0    , 0.294, 0.498, 0.725, 0.933, 1    , 1    , 1    ,
#          1    , 1    , 1    , 1    , 1    , 1    , 1    , 0.9  ,0.796,
#          0.7  , 0.569, 0.365,
#          0.0
#        )
#  B1.ref <- c( 0    , 0.5  , 0.750,
#          1    , 1    , 1    , 1    , 1    , 0.863, 0.659, 0.355,
#          0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    ,
#          0    , 0    , 0    , 0
#        )

  R1 <- c( 0    , 0.08 , 0.15 ,
          0    , 0    , 0    , 0.15    , 0.3    , 0.15    ,
          0.6  , 0.8  , 0.9  , 0.95 , 0.98 , 1    , 1    , 1    , 1    ,
          1    , 1    , 1    , 1
        )

  G1 <- c( 0    , 0    , 0    ,
          0    , 0.294, 0.498, 0.725, 0.933, 1    ,
          1    , 1    , 1    , 1    , 1    , 1    , 0.95 , 0.9  , 0.796,
          0.7  , 0.569, 0.365,
          0.0
        )

  B1 <- c( 0    , 0.5  , 0.750,
          1    , 1    , 1    , 1    , 1    , 0.659,
          0    , 0    , 0    , 0.2    , 0.4    , 0.7    , 0.2    , 0.2    , 0    ,
          0    , 0    , 0    , 0
        )
  if( n > 1 ) {

    X1 <- 1:length( R1 )
    X <- (0:(n-1)) * (length(R1)-1) / (n-1) + 1

    R <- approx(X1, R1, X )
    G <- approx(X1, G1, X )
    B <- approx(X1, B1, X )
  } else {
    return( rgb(0,0,0) )
  }

  if ( ! black ) {
    R$y <- R$y[-1]
    G$y <- G$y[-1]
    B$y <- B$y[-1]
  }

 if (format == "rgb" )
  val <- rgb(R$y, G$y, B$y)
 else
  val <- cbind( R$y, G$y, B$y )

  return(val)
}

color.lightening <- function( RGB, contrast=0.1) {

  n <- dim( RGB ) [1]
  for(i in 1:n) {
    dI <- min( RGB[i, ],       contrast )
    dJ <- min( 1.0 - RGB[i, ], contrast )

    if ( (dI == 0) && (dJ == 0) ) {
      if ( sum( RGB[i, ] ) > 1.5 ){
        RGB[i, ] <-  RGB[i, ] - contrast
      } else {
        RGB[i, ] <-  RGB[i, ] + contrast
      }
    }

    if( dI > dJ ){
      # Get darker
      RGB[i, ] <-  RGB[i, ] - dI
    } else {
      # Lightening
      RGB[i, ] <-  RGB[i, ] + dJ
    }
  }

  index <- which( RGB[] < 0.0 )
  RGB[index] <- 0.0

  index <- which( RGB[] > 1.0 )
  RGB[index] <- 1.0

  val <- rgb( RGB[,1], RGB[,2], RGB[,3] )
}

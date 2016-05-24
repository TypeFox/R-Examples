debug <- FALSE

dmsgf <- function( ... ) if( debug ) print( sprintf( ... ) )
dmsg  <- function( ... ) if( debug ) print( ... )
dmsgc <- function( ... ) if( debug ) cat( sprintf( ... ) )

# make sure that N1's position < N2's position
orderConnections <- function( x ) {

  pos <- x$nodes$x

  x$edges[ , c( "N1", "N2" ) ] <-
    t( apply( x$edges[ , c( "N1", "N2" ) ],
     1, function( x ) x[ order( pos[ as.character( x ) ] )  ] ) )

  return( x )
}

# adds new, invisible nodes to better direct output jumping over certain
# stops
add_mid_points <- function( x, default_style= NULL ) {

  # return a non-existing node name
  newnodename <- function( NN, nlist ) {
    i <- 0 ; nn <- paste0( NN, ".", i )
    while( nn %in% nlist ) { i <- i + 1 ; nn <- paste0( NN, ".", i ) }
    nn
  }

  l       <- names( x )
  all.pos <- sort( unique( x$nodes$x ) )
  rownames( x$nodes ) <- x$nodes$ID

  for( l1 in l[ order( x$nodes$x ) ] ) {

    # position of the first node
    p1 <- x$nodes[ l1, "x" ]


    for( p in all.pos[ all.pos < p1 - 1 ] ) {


      left_neighbors <- x$edges$N1[ x$edges$N2 == l1 & x$edges$Value > 0 ]
      todo <- left_neighbors[ x$nodes[ left_neighbors, "x" ] == p ]

      # nothing to do, move to the next position
      if( length( todo ) == 0 ) next ;

      newnode <- newnodename( l1, l )

      # add new row to nodes
      new.row <- x$nodes[ 1,,drop=F ][ NA, ]
      new.row[ , c( "ID", "x" ) ] <- list( ID= newnode, x= p + 1 )
      dmsg( new.row )
      # print( p + 1 )
      x$nodes <- rbind( x$nodes, new.row )
      rownames( x$nodes )[ nrow( x$nodes ) ] <- newnode
      # print( x$nodes )

      x$styles[[ newnode ]] <- getstyle( list( nodestyle= "invisible" ) )

      # we accumulate and add up all incoming traffic to the new node
      rsize <- 0
      for( t in todo ) {

        cur_e <- which( x$edges$N1 == t & x$edges$N2 == l1 )
        t.size <- x$edges$Value[ cur_e ]
        rsize <- rsize + t.size

        # rewire the connection to the new node
        x$edges$N2[ cur_e ] <- newnode 
      }

      # we need to calculate the node's edge colors. If there are multiple
      # nodes on the left side, we will take over the colors from the right

      if( length( todo ) > 1 ) {
        x$styles <- copyattr( x$styles, l1, newnode, "col" )
        dmsgf( "node %s: col %s", newnode, getattr( x$styles, newnode, "col" ) )
       #x$styles[[ newnode ]]$lcol <- 
       #  x$styles[[ newnode ]]$rcol <- 
       #  x$styles[[ l1 ]]$lcol
      } else {
        dp <- p1 - p + 1
        col1 <- getattr( x$styles, todo, "col" )
        col2 <- getattr( x$styles, l1,   "col" )
        tmp <- colorRampPaletteAlpha( c( col1, col2 ) )( dp )
        dmsgf( "node %s: col %s < %s < %s", newnode, col1, tmp[2], col2 )
        x$styles <- setattr( x$styles, newnode, "col", tmp[2] )
        #x$styles[[ newnode ]]$lcol <- 
          #x$styles[[ newnode ]]$rcol <- tmp[2]
      }

      eid <- newnodename( paste0( newnode, "->", l1 ), x$edges$ID )

      x$edges <- rbind( x$edges, rep( NA, ncol( x$edges ) ) )
      nn <- nrow( x$edges )
      rownames( x$edges )[nn] <- eid
      x$edges$ID[nn] <- eid
      x$edges$N1[nn] <- newnode
      x$edges$N2[nn] <- l1
      x$edges$Value[nn] <- rsize

      # update l
      l <- names( x )


    } # finish going over all position for the node l1
  }

  # print( x$nodes )
  # stop( "dupa" )
  return( x )
}

# calculate the sizes of the nodes, left, right and total
calcsizes2 <- function( x ) {

  e      <- x$edges
  nnames <- names( x )

  # total node sizes on the left and on the right
  lefts  <- sapply( nnames, function( n ) { sum( e$Value[ e$N2 == n ] ) } )
  names( lefts ) <- nnames
  rights <- sapply( nnames, function( n ) { sum( e$Value[ e$N1 == n ] ) } )
  names( rights ) <- nnames
  # sizey     <- sapply( nnames, function( n ) max( rights[n], lefts[n] ) )
  sizey <- apply( cbind( lefts, rights ), 1, max )
  names( sizey ) <- nnames
  pos.list  <- sort( unique( x$nodes$x ) )
  pos.maxes <- sapply( pos.list, function( p ) sum( sizey[ nnames[ x$nodes$x == p ] ] ) )

  return( list( lefts= lefts, rights= rights, sizey= sizey, pos.maxes= pos.maxes ) )
}


calcpos <- function( x, s, gravity= "top", node_margin= 0.1 ) {

  pos.list <- sort( unique( x$nodes$x ) )

  nnodes <- names( x )

  # calculate the distance between nodes
  max.y     <- max( s$pos.maxes )
  interstop <- max.y * node_margin # XXX
  dmsg( interstop )

  # the matrix will hold graphical coordinates of the positions
  pos.m <- matrix( 0, nrow= 4, ncol= length( nnodes ) )
  colnames( pos.m ) <- nnodes
  rownames( pos.m ) <- c( "x", "top", "center", "bottom" )
  pos.m[ "x", ]     <- x$nodes$x
  dmsg( x$nodes )

  for( i in 1:length( pos.list ) ) {
    p <- pos.list[ i ]

    cur.y <- 0
    nn <- nnodes[ x$nodes$x == p ] 
    #printf( "pos: %d", p )
    if( gravity == "top" ) nn <- rev( nn )
    for( n in nn ) {
      #print( n )
      if( ! is.null( x$nodes$y ) && ! is.na( x$nodes[n,]$y ) ) {
        pos.m[ "center", n ] <- x$nodes[n,]$y
        pos.m[ "top", n ]    <- x$nodes[n,]$y + s$sizey[ n ] / 2
        pos.m[ "bottom", n ] <- x$nodes[n,]$y - s$sizey[ n ] / 2
      } else if( gravity == "top" ) {
        #print( "top gravity" )
        pos.m[ "bottom", n ] <- cur.y
        pos.m[ "center", n ] <- cur.y - s$sizey[ n ] / 2
        pos.m[ "top"   , n ] <- cur.y - s$sizey[ n ]
        cur.y <- cur.y - s$sizey[ n ] - interstop
      } else if( gravity %in% c( "bottom", "center" ) ) {
        pos.m[ "bottom", n ] <- cur.y
        pos.m[ "center", n ] <- cur.y + s$sizey[ n ] / 2
        pos.m[ "top"   , n ] <- cur.y + s$sizey[ n ]
        cur.y <- cur.y + s$sizey[ n ] + interstop
      }
      #print( cur.y )
    }

  }


  if( gravity == "center" ) {
    ylim <- range( pos.m[ c( "bottom", "top" ), ] )
    dy <- ylim[2] - ylim[1]

    for( p in pos.list ) {
      nn <- nnodes[ x$nodes$x == p ]
      ylim2 <- range( pos.m[ c( "bottom", "top" ), nn ] )
      dy2   <- ylim2[2] - ylim2[1]

      pos.m[ c( "bottom", "center", "top" ), nn ] <- pos.m[ c( "bottom", "center", "top" ), nn ] - dy2/2
    }
  }


  #print( pos.m )

  return( pos.m )
}

draw.edges <- function( x, pos.m, s, col= "#ffcc33", lty= lty, nsteps= 50, boxw= 0.2 ) {
  # for each node, we need to to store the position of the current slot, on
  # the right and on the left

  w <- boxw / 2 
  pos.m <- rbind( pos.m, lpos= 0 )
  pos.m <- rbind( pos.m, rpos= 0 )

  pos.m[ "lpos", ] <- pos.m[ "center", ] - s$lefts  / 2
  pos.m[ "rpos", ] <- pos.m[ "center", ] - s$rights / 2

  for( i in 1:nrow( x$edges ) ) {
    n1 <- x$edges$N1[i]
    n2 <- x$edges$N2[i]
    id <- x$edges$ID[i]

    dx1 <- w
    if( getattr( x$styles, n1, "nodestyle" )  %in% c( "invisible", "point") ) dx1 <- 0
    #if( x$styles[[ n1 ]][[ "nodestyle" ]] %in% c( "invisible", "point" ) ) dx1 <- 0
    col1 <- getattr( x$styles, n1, "col" )
    #col1 <- x$styles[[n1]][["rcol"]] 

    dx2 <- w
    if( getattr( x$styles, n2, "nodestyle" )  %in% c( "invisible", "point") ) dx2 <- 0
    #if( x$styles[[ n2 ]][[ "nodestyle" ]] %in% c( "invisible", "point" ) ) dx2 <- 0
    #col2 <- x$styles[[n2]][["lcol"]] 
    col2 <- getattr( x$styles, n2, "col" )

    ss <- x$edges$Value[i]

    #print( x$styles[[n1]] )

    #if( x$styles[[ n1 ]][[ "edgestyle" ]] == 'straight' ) {
    if( getattr( x$styles, id, "edgestyle" ) == "straight" ) {
      form <- 'line' 
    } else {
      form <- 'sin' 
    }

    # determine the type of edge coloring coloring to use
    grad <- c( col1, col2 )
    col  <- NULL
    #printf( "edgecol: %s", getattr( x$styles, id, "edgecol" ) )
    if( getattr( x$styles, id, "edgecol" ) == "col" ) {
      grad <- NULL
      col  <- getattr( x$styles, id, "col" )
    }

    #print( grad )

    curveseg( pos.m[ "x", n1 ] + dx1, pos.m[ "x", n2 ] - dx2, 
               pos.m[ "rpos", n1 ],    pos.m[ "lpos", n2 ], 
               width= ss, grad= grad, col= col,
               lty= lty, nsteps= nsteps, form= form  )

    pos.m[ "rpos", n1 ] <- pos.m[ "rpos", n1 ] + ss
    pos.m[ "lpos", n2 ] <- pos.m[ "lpos", n2 ] + ss
    

  }

}

draw.nodes <- function( x, pos.m, s, width= 0.2, 
   lty= 1, col= NULL, srt= NULL, textcol= NULL, boxw= 0.2 ) {

  w <- boxw / 2

  for( n in names( x ) ) {
    if( getattr( x$styles, n, "nodestyle" ) == "invisible" ) next ;

    # if specific values are provided, they override the styles
    if( is.null( .lty <- lty ) ) .lty <- getattr( x$styles, n, "lty" )  
    if( is.null( .col <- col ) ) .col <- getattr( x$styles, n, "col" )  
    if( is.null( .srt <- srt ) ) .srt <- getattr( x$styles, n, "srt" )  
    if( is.null( .textcol <- textcol ) ) .textcol <- getattr( x$styles, n, "textcol" )

    if( is.null( x$nodes$labels ) || is.na( x$nodes[n,]$labels ) ) lab <- n
    else                                              lab <- x$nodes[n,]$labels
    
    if( getattr( x$styles, n, "nodestyle" ) == "point" ) {
      points( pos.m[ "x", n ], pos.m[ "center", n ], pch= 19, col= col )
    } else {

      rect( 
        pos.m[ "x", n ] - w, pos.m[ "bottom", n ], 
        pos.m[ "x", n ] + w, pos.m[ "top", n ], 
        lty= .lty, col= .col )
    }

    text( pos.m[ "x", n ], pos.m[ "center", n ], lab, col= .textcol, srt= .srt )
  }
}


# check whether vertical information is present
ypos_present <- function( x ) {
  if( is.null( x$nodes$y ) 
    || all( is.na( x$nodes$y ) ) ) return( FALSE )
  yrange <- range( x$nodes$y )

  if( yrange[1] == yrange[2] ) return( FALSE )
  TRUE
}

# scale node sizes automatically
autoscale <- function( x ) {

  if( ! ypos_present( x ) ) return( 1 )

  yrange <- range( x$nodes$y )
  yrange <- yrange[2] - yrange[1]

  ns <- max( x$edges$Value )

  if( ns == 0 ) return( 1 ) # not our problem

  yscale <- 0.15 * yrange / ns
  return( yscale )
}

#' @rdname riverplot
#' @method plot riverplot
#' @S3method plot riverplot
#' @export
plot.riverplot <- function( x, ...  ) riverplot( x, ... )

#' Create a Sankey plot
#'
#' Create a Sankey plot
#'
#' This functions create a Sankey plot given a riverplot object
#' (\code{plot} is just a wrapper for the \code{riverplot} function.
#' The object to be drawn is a list specifying the plot; see
#' the \code{\link{makeRiver}} function for exact specifications and
#' the \code{\link{riverplot.example}} to see how it can be created.
#' Whether or not the list used to plot is exactly of class
#' \code{riverplot-class} does not matter as long as it has the correct
#' contents. 
#'
#' Style information which is missing from the riverplot object \code{x} (for example, if the
#' node style is not specified for each node in the object) is taken from the \code{default.style} parameter.
#' See \code{\link{riverplot-styles}()} and
#' \code{\link{updateRiverplotStyle}()} to learn how to create and
#' modify the styles.
#'
#' Whether or not the list used to plot is exactly of class
#' \code{riverplot-class} does not matter as long as it has the correct
#' contents. These functions here are for the convenience of checking that
#'
#' The nodes are drawn from bottom to top in the order they are found in
#' the riverplot object. There is no clever algorithm for placing the nodes
#' minimizing the number of crossing edges yet; you need to manipulate the
#' object directly to achieve the desired effect.
#'
#'@param x An object of class riverplot
#'@param lty Line style to use
#'@param srt Rotation of the node labels
#'@param default_style default graphical style
#'@param gravity how the nodes are placed vertically. No effect if node
#'       vertical positions are specified via \var{node_ypos} member
#'@param node_margin how much vertical space should be kept between the
#'       nodes
#'@param nodewidth width of the node (relative to font size)
#'@param plot_area fraction of vertical space to be used as main plot area
#'@param nsteps number of interpolating steps in drawing the segments
#'@param add_mid_points attempt to get a smoother plot by adding additional
#'       nodes. Set this parameter to \code{FALSE} if you are setting node
#'       vertical position manually. If add_mid_points is equal to NULL (the
#'       default), then the
#'       mid points are added only if \var{node_ypos} is empty.
#'@param yscale scale the edge width values by multiplying with this
#'       factor. If \var{yscale} is equal to "auto", scaling is done
#'       automatically such that the vertical size of the largest node is
#'       approximately equal to 15% of the range of ypos (if present). 
#'       If no \var{node_ypos} is specified in the riverplot object, no scaling is
#'       done.
#'       If \var{yscale} is equal to 1, no scaling is done.
#'@param ... further plotting parameters
#'@return \code{riverplot} return invisibly a matrix containing the
#'       actual positions (in user coordinates) of the nodes drawn on the screen.
#'       Note that it also may contain additional, invisible nodes that have been
#'       created by the algorithm to better fit on the screen.
#'@seealso riverplot-styles updateRiverplotStyle minard
#'@examples 
#' x <- riverplot.example()
#' plot(x)
#' plot(x, srt=90, lty=1)
#'@export
riverplot <- function( x, lty= 0, srt= NULL, 
                             default_style= NULL, gravity= "top", 
                             node_margin= 0.1, 
                             nodewidth= 1.5,
                             plot_area= 0.5,
                             nsteps= 50,
                             add_mid_points= NULL,
                             yscale= "auto" 
                             ) {

  default_style <- getstyle( NULL, default_style )
  plot.new()
  dmsgc( "--------------\nDefault style:\n-----------\n" )
  dmsg( default_style )
  dmsgc( "--------------\n" )

  x2 <- x
  x2$nodes$ID <- as.character( x2$nodes$ID )

  # check sanity of the edge information
  dmsg( "checking edges" )
  x2$edges <- checkedges( x2$edges, names( x2 ) )

  # N1 must be the node on the left, N2 on the right
  x2 <- orderConnections( x2 )

  # add mid points automatically depending whether y-positions of nodes are specified
  if( is.null( add_mid_points ) ) {
    if( ypos_present( x2 ) ) add_mid_points <- FALSE
    else                     add_mid_points <- TRUE
  }

  # add mid points
  dmsg( "adding mid points" )
  if( add_mid_points ) x2 <- add_mid_points( x2 )

  # update styles
  for( n in c( x2$nodes$ID, x2$edges$ID ) ) { 
    x2$styles[[ n ]] <- getstyle( x2$styles[[ n ]], default_style, update.missing= FALSE ) 
  }

  dmsgc( "Updated styles:\n" )
  dmsg( x2$styles )
  dmsgc( "--------------\n" )

  #for( n in names( x2$nodes ) ) { x2$styles[[ n ]] <- getstyle( x2$styles[[ n ]], default.style ) }

  if( yscale == "auto" ) yscale <- autoscale( x2 )
  if( yscale != 1 ) x2$edges$Value <- x2$edges$Value * yscale

  dmsg( "calculating sizes" )
  sizes <- calcsizes2( x2 )
  dmsg( sizes )

  dmsg( "calculating positions" )
  positions <- calcpos( x2, sizes, gravity= gravity, node_margin= node_margin )
  dmsg( "done" )
  xrange <- range( x2$nodes$x )
  xlim <- xrange + (xrange[2]-xrange[1]) * c( -0.1, 0.1 ) 
  ylim <- range( positions[ c( "bottom", "top" ), ] )
  b <- (ylim[2] - ylim[1]) * (1-plot_area)/plot_area / 2
  ylim <- ylim + c( -b, b )

  par( mar= c( 0, 0, 0, 0 ) )
  par( usr= c( xlim, ylim ) )

  l <- names( x2$nodes )[ order( x2$nodes ) ]

  dev.hold()
  on.exit( dev.flush() )

  w <- strwidth( "hjKg" ) * nodewidth / 2

  dmsg( "drawing edges" )
  draw.edges( x2, positions, sizes, lty= lty, nsteps= nsteps, boxw= w )
  dmsg( "drawing nodes" )
  draw.nodes( x2, positions, sizes, srt= srt, boxw= w, lty= lty )

  return( invisible( positions ) )
}

Gplot.env <- new.env( parent=baseenv() )
assign("Gplot.graph",
       list( pie.gr=NULL, default.gr=NULL, reso=2000,
            col= sample( light.palette(10,black=FALSE) ),
            dust.bin = NULL ),
           envir=Gplot.env ) 
# __________________________________________________________
##
## Gplot
##
## Public
##
## author : G. Grasseau
##
## INPUT : 
##	Matrix   : adjacency/valued matrix
##                 or edge list
##      type     : type of display ("default, "pie.nodes")
##
## OUTPUT :
##	gplot node coordinates
##
## __________________________________________________________
##
Gplot <- function ( X, type="default", coord=NULL,... ) {

  Env <- get("Gplot.graph", envir=Gplot.env)
  res <- -1

  if( type == "default" ) {
    same.node.pos <- all ( dim( Env$default.gr$mat) ==  dim(X) )
    if( ( is.null(X) ) | ! same.node.pos ) {
      # Clear memory of stored node positions
      Env$default.gr$network$coord=NULL
    }
    if ( ! is.null(X) & ! is.null(coord) ) {
      res <- Gplot.default(  X,  coord=coord, ...)
    } else if ( ! is.null(X) ) {
      res <- Gplot.default(  X,  coord=Env$default.gr$network$coord, ...)
    }
  } else if(  type == "pie.nodes" ) {
    same.node.pos <- all ( dim( Env$pie.gr$mat) ==  dim(X) )
    if( ( is.null(X) ) | !same.node.pos ) {
      # Clear memory of stored node positions
      Env$pie.gr$network$coord=NULL
    }
    if ( ! is.null(X) ) {
      res <- Gplot.pie.nodes( X, coord=Env$pie.gr$network$coord, ... )
    }
  } else {
    cat("Gplot error: bad type of display <" , type, "> \n")
  }
  
  invisible(res)
}

default.color.map <- function( X ) {
  e.col <- X
  e.col[X<0] <- "red"  # red : negative 
  e.col[X>0] <- "blue" # blue : positive
  e.col
 }

## __________________________________________________________
##
## Gplot.default
##
## Public
##
## author : A. Smith and G. Grasseau
##
## INPUT : 
##    Matrix   : adjacency/valued matrix
##               or edge list
##
##    class     : node classification vector 
##    class.col : colours to use (per class if classes is specified) 
##    coord	: matrix of node coordinates for plot
##    label    : vector of node labels (override dimnames(Matrix))
##    threshold : threshold under which node degrees are ignored.
##                The value 0 means that all nodes are displayed
##
##    max.edges : max number of edges below which the graph is plotted
##
##    main	: graphical parameter main
##    sub       : graphical parameter sub
##
## OUTPUT :
##	gplot node coordinates
## __________________________________________________________
##
Gplot.default <- function ( X, ... ) {

  
  #    Defaults for hidden optional arguments 
  #    --------------------------------------
  
  #  Colours for each class
  directed		<- sub.param("directed"		, NULL	, ...)

  class 		<- sub.param("class"		, NULL	, ...)

  #  Colours for each class
  class.col		<- sub.param("class.col"	, NULL	, ...)
  
  #  Coordinates for nodes
  coord			<- sub.param("coord"		, NULL	, ...)
  
  #  Coordinates for nodes
  color.map	        <- sub.param("color.map"	, default.color.map , ...)

  #  Labels for nodes                                       
  label		<- sub.param("label"	, NULL	, ...)
  label.display=FALSE
  if( ! is.null(label) ) {
    if( label[1] == "default" ) {
      label.display=TRUE
      # Get dimnames as label
      label <- NULL
    } else if ( is.character(label) )
      label.display=TRUE
  }

  # Max. number of plottable edges
  max.edges		<- sub.param("max.edges"	, 10000	, ...)
  
  # Graph title
  main			<- sub.param("main"		, "Gplot", ...)
  
  # Graph subtitle
  sub			<- sub.param("sub"		, NULL	, ...)	

  #    Check arguments
  #    ---------------

  no.edge <-FALSE
  if ( sum( (X > .Machine$double.eps) )==0) {
    no.edge <-TRUE
  }

  # Edge matrix case
  if ( dim(X)[1] == 2) {
    if ( dim(X)[2] == 2 ) {
      cat("Gplot warning : Matrix(2,2), treated as an adjacency matrix \n")
    } else {
      NNodes <- max(X)
      Y <- matrix(0, NNodes, NNodes )
      Y[ cbind( X[1,], X[2,] ) ] <- 1
      if ( ! directed )
          Y[ cbind( X[2,], X[1,] ) ] <- 1
       X <- Y
    }
  }

  #     Restore context 
  #     ---------------
  Env <- get("Gplot.graph", envir=Gplot.env)
  save.class.col <- Env$col
  
  #  Matrix check
  NNodes <- dim(X)[1]
  if (NNodes != dim(X)[2]) {
    cat("Gplot error : X matrix must be square \n")
    return(coord)
  }
  
  # Set "directed" parameter when not specified
  if ( is.null( directed ) ) {
    symetric <- (is.symetric(X, eps=1.0e-3))
    if ( symetric ) {
       X[lower.tri(X)] <- 0
       directed <- FALSE
    } else {
       directed <- TRUE
    }
  }

  # Check valid labels
  if( (length( label ) != NNodes ) && ( length( label ) != 0 ) ) {
     cat("Gplot error : length of  <label> differs from the matrix size \n")
     return(coord)
  }
  
  # Check valid class
  if ( length(class) != NNodes) {
    if ( length(class) == 0)
      class <-1
    class <- rep(class, NNodes) [1: NNodes]
  }
  
  # Check valid class.col
  class <- as.factor(class)
  nb.lev <- length( levels(class) )
  if ( nb.lev > length( class.col )) {
    class.col <- sample( light.palette( nb.lev, black=FALSE ) )
  }
  vertex.col <- class.col[as.numeric(class)]

  coord.save <- coord
  
  #    Process dustbin nodes
  #    ---------------------

  # Identify dustbin nodes
  degree.threshold <- 0
  b.dust <- rep(FALSE, NNodes)
  b.dust <- ( getDegree( X, upper=TRUE ) < degree.threshold )
  Env$dust.bin <- b.dust
  
  dustbin.not.empty <- any( b.dust )

  if( sum(b.dust) == NNodes ) {
    cat("Gplot warning : no node to display \n")
    return(coord)
  }

  
  # Remove dustbin nodes
  if (dustbin.not.empty) {
    X <- X[!b.dust,!b.dust]
    NNodes <- dim(X)[1]
    if (!is.null(coord)) {
      coord <- coord [!b.dust,]
    }
    if (!is.null(class)) {
      class <- class [!b.dust]
    }
    if (!is.null(label)) {
      label <- label [!b.dust]
    }
    vertex.col <- vertex.col[!b.dust]
  }

  #    Set graphic parameters
  #    ----------------------


  # Edge color handling
  e.col <- color.map( X )

  # Edge width
  lwd <- array(0.5, dim( X ) )
  #   lwd <- abs(X) 
  #   
  #   if( length(lwd) > 100 ) {
  #     # Number of edges to high ( width set to 1) )
  #     lwd[] <- ( lwd[] > 0 ) * 1
  #   } else {
  #     # max line weight is 5
  #     lwd <- lwd * 5 / max( lwd )
  #   }

  # Edge type
  e.lty <- "solid"

  # Get node labels
  if ( is.null(label) ) {
   label <- dimnames(X)[[1]]
  }

  
  # Strictly undirected graph, matrix considered symmetric
  if ( ! directed ) {
    X[lower.tri(X)] <- 0
  }
  
  X <- abs(X)

  # Plot if not too many edges
  num.edges <- sum( X != 0 )
  if ( num.edges > max.edges ) {

    cat( "Gplot warning : too many edges to plot (",
         num.edges,">",max.edges,") \n")

  } else {
    gr <- Gplot.graphics( X,
                          thresh		= 0,	    # no edge blocking
                          margin		= 0.1,
                          main		= main,
                          sub		= sub
                        )

    f.class <- factor( class )
    class.dims <- rep(0, nlevels( f.class ))
    for (i in 1:nlevels( f.class ) ) {
       class.dims[i] <- length( which( f.class == levels(f.class)[i] ) )
    }

    gr <- Gplot.network( gr,
                         loop.draw       = directed,
                        # Warning : conflict with coord
                         vertex.pos.mode = "circle",
                         displayisolates = TRUE,
                         coord           = coord,
                         class = class
                        )
    
    if( ! no.edge ) 
      gr <- Gplot.edge( gr,
                       col       = e.col,
                       lty       = e.lty,
                       lwd       = lwd,
                       curved    = (NNodes < 60) ,
                       arrow.draw = directed,
                       arrow.cex  = 10,
                       loop.draw  = directed,
                       loop.cex   = 5
                       )

    gr <- Gplot.vertex( gr, col= vertex.col )

    if ( label.display )
      gr <- Gplot.vertex.label( gr,
                                label = label,
                                cex   = 0.7,
                                useboxes = FALSE
                               )
  }
  #     Save context 
  #     ------------
  Env$col <- class.col
  Env$default.gr <- gr
  assign("Gplot.graph", Env, envir=Gplot.env )
  
  invisible(coord.save)
  
}


## __________________________________________________________
##
## Gplot.add
##
## Public
##
## author : G. Grasseau
##
##   Add a nex edges/arc in the graph plot
##
## INPUT : 
##    Matrix   : adjacency/valued matrix
##               or edge list
##
##    label    : vector of node labels (override dimnames(Matrix))
##    class     : node classification vector 
##    class.col : colours to use (per class if classes is specified)
##    label    : vector of node labels (override dimnames(Matrix))
##    max.edges : max number of edges below which the graph is plotted
##
##    main	: graphical parameter main
##    sub       : graphical parameter sub
##
## OUTPUT :
##	gplot node coordinates
## __________________________________________________________
##
Gplot.add <- function ( X, ... ) {

  #    Defaults for hidden optional arguments 
  #    --------------------------------------
  
  #  Colours for each class
  directed		<- sub.param("directed"		, NULL	, ...)

  class 		<- sub.param("class"		, NULL	, ...)

  #  Colours for each class
  class.col		<- sub.param("class.col"	, NULL	, ...)
  
  #  Labels for nodes                                       
  label		<- sub.param("label"	, NULL	, ...)
  label.display=FALSE
  if( ! is.null(label) ) {
    if( label[1] == "default" )
      label.display=TRUE
    else if ( is.character(label) )
      label.display=TRUE
  }

  #                                   ...)
  # Max. number of plottable edges
  max.edges		<- sub.param("max.edges"	, 10000	, ...)
  
  # Graph title
  main			<- sub.param("main"		, "Gplot", ...)
  
  # Graph subtitle
  sub			<- sub.param("sub"		, NULL	, ...)	

  #    Check arguments
  #    ---------------

  no.edge <-FALSE
  if ( sum( ( abs( X ) > .Machine$double.eps) )==0) {
    no.edge <-TRUE
  }

  # Edge matrix case
  if ( dim(X)[1] == 2) {
    if ( dim(X)[2] == 2 ) {
      cat("Gplot warning : Matrix(2,2), treated as an adjacency matrix \n")
    } else {
      NNodes <- max(X)
      Y <- matrix(0, NNodes, NNodes )
      Y[ cbind( X[1,], X[2,] ) ] <- 1
      if ( ! directed )
          Y[ cbind( X[2,], X[1,] ) ] <- 1
       X <- Y
    }
  }

  #     Restore context 
  #     ---------------
  Env <- get("Gplot.graph", envir=Gplot.env)
  save.class.col <- Env$col
  
  #  Matrix check
  NNodes <- dim(X)[1]
  if (NNodes != dim(X)[2]) {
    cat("Gplot error : X matrix must be square \n")
  }
  
  # Set "directed" parameter when not specified
  if ( is.null( directed ) ) {
    symetric <- (is.symetric(X, eps=1.0e-3))
    if ( symetric ) {
       X[lower.tri(X)] <- 0
       directed <- FALSE
    } else {
       directed <- TRUE
    }
  }

  # Check valid labels
  if( (length( label ) != NNodes ) && ( length( label ) != 0 ) ) {
     cat("Gplot error : length of  <label> differs from the matrix size \n")
  }
  # Get node labels
  if ( is.null(label) ) {
   label <- dimnames(X)[[1]]
  }
  
  # Check valid class
  if ( length(class) != NNodes) {
    if ( length(class) == 0)
      class <-1
    class <- rep(class, NNodes) [1: NNodes]
  }
  
  # Check valid class.col
  class <- as.factor(class)
  nb.lev <- length( levels(class) )
  if ( nb.lev > length( class.col )) {
    class.col <- sample( light.palette( nb.lev, black=FALSE ) )
  }
  vertex.col <- class.col[as.numeric(class)]

  #    Process dustbin nodes
  #    ---------------------

  # Restore dustbin nodes
  b.dust <- Env$dust.bin
  
  dustbin.not.empty <- any( b.dust )

  # Remove dustbin nodes
  if (dustbin.not.empty) {
    X <- X[!b.dust,!b.dust]
    NNodes <- dim(X)[1]
    if (!is.null(class)) {
      class <- class [!b.dust]
    }
    if (!is.null(label)) {
      label <- label [!b.dust]
    }
     vertex.col <- vertex.col[!b.dust]
  }
  
  #    Set graphic parameters
  #    ----------------------


  # Edge color handling
  e.col <- X
  e.col[X<0] <- "red"  # red : negative 
  e.col[X>0] <- "blue" # blue : positive 

  # Edge width
  lwd <- array(0.5, dim( X ) )

#  if( length(lwd) > 100 ) {
#    # Number of edges to high ( width set to 1) )
#    lwd[] <- ( lwd[] > 0 ) * 1
#  } else {
#    # max line weight is 5
#    lwd <- lwd * 5 / max( lwd )
#  }

  # Edge type
  e.lty <- "solid"

  # Add new edges
  
  # Strictly undirected graph, matrix considered symmetric
  if ( ! directed ) {
    X[lower.tri(X)] <- 0
  }
  # Save new + old matrix in raw field
  gr     <- Env$default.gr
  new.gr <- gr 
  new.gr$mat.raw[ which( X !=0 ) ] <- X[ which( X !=0 ) ]
  
  X <- abs(X)

  # Plot if not too many edges
  num.edges <- sum( X != 0 )
  if ( num.edges > max.edges ) {

    cat( "Gplot warning : too many edges to plot (",
         num.edges,">",max.edges,") \n")
  } else {

    # Draw new edges
    gr$mat <- X
    if( ! no.edge ) {
      gr <- Gplot.edge( gr,
                       col       = e.col,
                       lty       = e.lty,
                       lwd       = lwd,
                       curved    = (NNodes < 60) ,
                       arrow.draw = directed,
                       arrow.cex  = 10,
                       loop.draw  = directed,
                       loop.cex   = 5
                       )
    }
    
    gr <- Gplot.vertex( gr, col= vertex.col )

    gr <- Gplot.vertex.label( gr,
                             label = label,
                             cex   = 0.7,
                             useboxes = FALSE
                             )
  }
  #     Save context 
  #     ------------
  Env$col <- class.col
  Env$default.gr <- new.gr
  assign("Gplot.graph", Env, envir=Gplot.env )
  
}


## __________________________________________________________
##
## Gplot.delete
##
##   Remove edges/arc in a graph.
##
## Public
##
## author : A. Smith and G. Grasseau
##
## INPUT : 
##    Matrix   : adjacency/valued matrix
##               or edge list
##
##    del.col  : edge removed are colored in del.col
## 
##    class     : node classification vector 
##    class.col : colours to use (per class if classes is specified) 
##    label    : vector of node labels (override dimnames(Matrix))
##
##    max.edges : max number of edges below which the graph is plotted
##
##    main	: graphical parameter main
##    sub       : graphical parameter sub
##
## OUTPUT :
##	empty
## __________________________________________________________
##
Gplot.delete <- function ( X, del.col="white", ... ) {

  #    Defaults for hidden optional arguments 
  #    --------------------------------------
  
  #  Colours for each class
  directed		<- sub.param("directed"		, NULL	, ...)

  class 		<- sub.param("class"		, NULL	, ...)

  #  Colours for each class
  class.col		<- sub.param("class.col"	, NULL	, ...)
  
  #  Labels for nodes                                       
  label		<- sub.param("label"	, NULL	, ...)
  label.display=FALSE
  if( ! is.null(label) ) {
    if( label[1] == "default" )
      label.display=TRUE
    else if ( is.character(label) )
      label.display=TRUE
  }


  # Max. number of plottable edges
  max.edges		<- sub.param("max.edges"	, 10000	, ...)
  
  # Graph title
  main			<- sub.param("main"		, "Gplot", ...)
  
  # Graph subtitle
  sub			<- sub.param("sub"		, NULL	, ...)	

  #    Check arguments
  #    ---------------

  no.edge <-FALSE
  if ( sum( ( abs( X ) > .Machine$double.eps) )==0) {
    no.edge <-TRUE
  }

  # Edge matrix case
  if ( dim(X)[1] == 2) {
    if ( dim(X)[2] == 2 ) {
      cat("Gplot warning : Matrix(2,2), treated as an adjacency matrix \n")
    } else {
      NNodes <- max(X)
      Y <- matrix(0, NNodes, NNodes )
      Y[ cbind( X[1,], X[2,] ) ] <- 1
      if ( ! directed )
          Y[ cbind( X[2,], X[1,] ) ] <- 1
       X <- Y
    }
  }

  #     Restore context 
  #     ---------------
  Env <- get("Gplot.graph", envir=Gplot.env)
  save.class.col <- Env$col
  
  #  Matrix check
  NNodes <- dim(X)[1]
  if (NNodes != dim(X)[2]) {
    cat("Gplot error : X matrix must be square \n")
  }
  
  # Set "directed" parameter when not specified
  if ( is.null( directed ) ) {
    symetric <- (is.symetric(X, eps=1.0e-3))
    if ( symetric ) {
       X[lower.tri(X)] <- 0
       directed <- FALSE
    } else {
       directed <- TRUE
    }
  }

  # Check valid labels
  if( (length( label ) != NNodes ) && ( length( label ) != 0 ) ) {
     cat("Gplot error : length of  <label> differs from the matrix size \n")
  }
  
  # Check valid class
  if ( length(class) != NNodes) {
    if ( length(class) == 0)
      class <-1
    class <- rep(class, NNodes) [1: NNodes]
  }
  
  # Check valid class.col
  class <- as.factor(class)
  nb.lev <- length( levels(class) )
  if ( nb.lev > length( class.col )) {
    class.col <- sample( light.palette( nb.lev, black=FALSE ) )
  }
  vertex.col <- class.col[as.numeric(class)]

  #    Process dustbin nodes
  #    ---------------------
  
  # Restore dustbin nodes
  b.dust <- Env$dust.bin
  
  dustbin.not.empty <- any( b.dust )

  # Remove dustbin nodes
  if (dustbin.not.empty) {
    X <- X[!b.dust,!b.dust]
    NNodes <- dim(X)[1]
    if (!is.null(class)) {
      class <- class [!b.dust]
    }
    if (!is.null(label)) {
      label <- label [!b.dust]
    }
    vertex.col <- vertex.col[!b.dust]
  }

  #    Set graphic parameters
  #    ----------------------


  # Edge color handling
  e.col <- X
  e.col[X!=0] <- del.col  
  
  # Edge width
  lwd.del <- array(10.0, dim( X ) )
  lwd <- array(0.5, dim( X ) )
  
  # Edge type
  e.lty <- "solid"

  # Get node labels
  if ( is.null(label) ) {
   label <- dimnames(X)[[1]]
  }

  # Strictly undirected graph, matrix considered symmetric
  if ( ! directed ) {
    X[lower.tri(X)] <- 0
  }

  # Save previous graph (matrix) and set the new one
  gr     <-  Env$default.gr
  old.gr <- gr 
  old.X <- X
  gr$mat <- abs(X)

  # Plot if not too many edges
  num.edges <- sum( X != 0 )
  if ( num.edges > max.edges ) {

    cat( "Gplot warning : too many edges to plot (",
         num.edges,">",max.edges,") \n")
  } else {
    if( ! no.edge ) {
      
      # Edges/arcs to delete
      # Bypas to erase all points
      gr$network$arrow.angle  = 40*pi/180
      gr$network$arrow.reso  = 60
      gr <- Gplot.edge( gr,
                       col       = e.col,
                       lty       = e.lty,
                       lwd       = lwd.del,
                       curved    = (NNodes < 60) ,
                       arrow.draw = directed,
                       arrow.cex  = 10, 
                       loop.draw  = directed,
                       loop.cex   = 5
                       )
      gr$network$arrow.angle  = 40*pi/180
      gr$network$arrow.reso  = 40
      gr <- Gplot.edge( gr,
                       col       = e.col,
                       lty       = e.lty,
                       lwd       = lwd,
                       curved    = (NNodes < 60) ,
                       arrow.draw = directed,
                       arrow.cex  = 10, 
                       loop.draw  = directed,
                       loop.cex   = 5
                       )
     old.gr$network$arrow.angle  = 20*pi/180
    
    #
    # Redraw old edges
    #

    # Resotre the previous matrix to redraw
    old.gr$mat.raw[ which( X != 0 ) ] <- 0
    
    # Edge color handling
    
    e.col <- old.gr$mat.raw
    e.col[ old.gr$mat.raw < 0 ] <- "red"  # red : negative 
    e.col[ old.gr$mat.raw > 0 ] <- "blue" # blue : positive
    
    old.gr$mat <- abs( old.gr$mat.raw )

      
    old.gr <- Gplot.edge( old.gr,
                       col       = e.col,
                       lty       = e.lty,
                       lwd       = lwd,
                       curved    = (NNodes < 60) ,
                       arrow.draw = directed,
                       arrow.cex  = 10,
                       loop.draw  = directed,
                       loop.cex   = 5
                       )

      
      old.gr <- Gplot.vertex( old.gr, col= vertex.col )

      old.gr <- Gplot.vertex.label( old.gr,
                             label = label,
                             cex   = 0.7,
                             useboxes = FALSE
                             )
    }

  }
  #     Save context 
  #     ------------
  Env$col <- class.col
  Env$default.gr <- old.gr
  assign("Gplot.graph", Env, envir=Gplot.env )
  
}



Gplot.pie.nodes <- function ( Mat, 
		    ... ) {
 
  #    Defaults for hidden optional arguments 
  #    --------------------------------------
  
  #  Colours for each class
  directed		<- sub.param("directed"		, NULL	, ...)

  #  Node weight - the node drawing is proportionnal to this weight  
  node.weight			<- sub.param("node.weight", 1, ...) 
  
  #  Pie coefficients 
  node.pie.coef			<- sub.param("node.pie.coef"	  , 1	, ...)

  #  Coordinates for nodes
  coord			<- sub.param("coord"		, NULL	, ...)
  
  #  Labels for nodes                                       
  label		<- sub.param("label"	        , NULL	, ...)

  label.display=FALSE
  if( ! is.null(label) ) {
    if( label == "default" )
      label.display=FALSE
    else if ( is.character(label) )
      label.display=FALSE
  }
  
  # Quantile value under which edges are considered.
  quantile.val	<- sub.param("quantile.val"	, 0.0,	  ...)
  
  # Max. number of plottable edges
  max.edges		<- sub.param("max.edges"	, 10000	, ...)
  
  # Graph title
  main			<- sub.param("main"		, "Gplot", ...)
  
  # Graph subtitle
  sub			<- sub.param("sub"		, NULL	, ...)	


  #     Restore context 
  #     ---------------
  Env <- get("Gplot.graph", envir=Gplot.env)
  save.class.col <- Env$col
  
  no.edge <-FALSE
  if ( sum( (Mat >.Machine$double.eps) ) == 0) {
    no.edge <-TRUE
  }

  # Dimension check
  p <- dim(Mat)[1]
  if (p != dim(Mat)[2]) {
    cat("Gplot warning : Mat matrix must be squared \n")
    return(coord)
  }
  NNodes <- p

  # Set default colors for pie nodes ie automatics
  pie.col <- -1
  # Set "directed" parameter when not specified
  if ( is.null( directed ) ) {
    symetric <- (is.symetric(Mat, eps=1.e-6) )
    if ( symetric ) {
       Mat[lower.tri(Mat)] <- 0
       directed <- FALSE
    } else {
       directed <- TRUE
    }
  }
  
  # Check valid node.weight
  if ( length(node.weight) != NNodes) {
    node.weight[] <- rep(node.weight, NNodes) [1: NNodes]
  }
  
  # Check valid node.pie.coef
  if ( is.vector( node.pie.coef ) ){
    if( length( node.pie.coef ) == 1) {
      if( length(save.class.col) < NNodes ) {
        save.class.col <- sample( light.palette( NNodes, black=FALSE ) )

      }
      pie.col <- save.class.col
      node.pie.coef <- diag( 1, NNodes, NNodes)
    }
    else {
      # Not tested
      node.pie.coef <- matrix( rep( node.pie.coef, NNodes),
                            NNodes, length( node.pie.coef ),  byrow=TRUE )
    }
  } else if ( is.matrix( node.pie.coef ) ) {
    dim.coef <- dim( node.pie.coef )
    if ( dim.coef[2] > 30 )  {
      cat("Gplot error : second dimension of node.pie.coef matrix too large \n")
      return(coord)
    }
    if( dim( node.pie.coef)[1] != NNodes ) {
       node.pie.coef <- rep( node.pie.coef, NNodes) [1: NNodes, dim.coef[2]]
    }
    
    # Select palette
    if( length(save.class.col) < dim(node.pie.coef)[2] ) {
        save.class.col <- sample( light.palette( dim(node.pie.coef)[2],
                                                 black=FALSE ) )

      }
    pie.col <- save.class.col

  } else {
    cat("Gplot error : a matrix is required for node.pie.coef \n")
    return(coord)
  }

  # Check valid labels
  if( (length( label ) != NNodes ) && ( length( label ) != 0 ) ) {
     cat("Gplot error : length of  <label> differs from the matrix size \n")
     return(coord)
  }
  
  # Node weight
  weight.max <- max( node.weight[] )
  node.weight[ which( node.weight[] <= 0 ) ] <- 0
  if ( weight.max > 0)
    node.weight[] <-  node.weight[] / weight.max
  max.node.weight <- weight.max 
  min.node.weight <- 0
  
  # Node coefficient
  n.nodes <- dim( node.pie.coef )[1]
  node.pie.coef[ which( node.pie.coef[] <= 0 ) ] <- 0
  l.sum <- rowSums( node.pie.coef )  
  for ( i in 1:n.nodes) {
    if ( l.sum[i] > 0)
      node.pie.coef[i,] <-  node.pie.coef[i,] / l.sum[i]
  }
  
  # Quantile filter
  if ( directed ) {
    threshold <- quantile( Mat[ ], probs=c(quantile.val) )
  }else{
    threshold <- quantile( Mat[ upper.tri( Mat[], diag=TRUE ) ],
                           probs=c(quantile.val) )
  }
  indexes <- which(Mat[] < threshold )
  Mat[indexes] <- 0

  # gplot threshold and line weight handling
  e.col <- "blue"
  # e.col[X<0] <- "red"  # red : negative partial correlation (inhibition)
  # e.col[X>0] <- "blue" # blue : positive partial correlation (activation)
  # e.col[X==0]<- "black"	


  # Line width
  # Max line width is 40
  abs.Mat <- abs( Mat )
  max.Mat <- max( abs.Mat )
  min.Mat <- min( abs.Mat )
  e.lwd <- matrix(0, dim(Mat)[1], dim(Mat)[2] )
  if ( abs(max.Mat - min.Mat) < .Machine$double.eps ) {
    e.lwd[] <- 1 * (Mat[] >  .Machine$double.eps )
  } else {
    norm.Mat <-  40 / ( max.Mat - min.Mat )
    e.lwd[] <- (Mat[] -  min.Mat) * norm.Mat
  }

  e.lty <- "solid"

  # Get node label
  if (is.null(label)) {
   label <- dimnames(Mat)[[1]]
  } 

  # Strictly undirected graph, matrix considered symmetric
  if ( ! directed ) {
    Mat[ lower.tri(Mat) ] <- 0
  }

  # Only plot if not too many edges
  num.edges <- sum( Mat > 0)

  if ( num.edges > max.edges ) {

    cat( "Gplot warning : too many edges to plot (",
         num.edges,">",max.edges,") \n")
    g  <- coord
    coord.save <- g

  } else {
    gr <- Gplot.graphics( Mat,
                          thresh		= 0,	    # no edge blocking
                          margin		= 0.3,
                          main		= main,
                          sub		= sub
                        )
    gr <- Gplot.network( gr,
                         loop.draw       = FALSE,
                         displayisolates = TRUE,
                         coord           = coord
                       )

    sqrt.min.Mat    <- sqrt( min.Mat)
    sqrt.min.node.weight <- sqrt( min.node.weight )
    cst <- 5 * 1 / ( sqrt( max.node.weight ) - sqrt.min.node.weight )

    gr <- Gplot.vertex.pie( gr,
                           cex = (sqrt( node.weight )*10), display=FALSE )

   if( ! no.edge )
     gr <- Gplot.edge( gr,
                      col        = "blue",
                      lty        = e.lty,
                      lwd        = e.lwd,
                      curved     = (NNodes < 60) ,
                      arrow.draw = directed,
                      arrow.cex  = 3.0,
                      loop.draw  = TRUE,
                      loop.cex   = 3.0
                     )

    gr <- Gplot.vertex.pie( gr, col = pie.col,
                           cex = (sqrt( node.weight )*10), pie.coef=node.pie.coef )
    
    gr <- Gplot.vertex.label( gr,
                                label=label,
                                cex=0.7,
                                useboxes=FALSE
                               )
    coord.save <- gr$network$coord

  }
  
  #     Save context 
  #     ------------
  Env$col <- save.class.col
  Env$pie.gr <- gr
  assign("Gplot.graph", Env, envir=Gplot.env )
  
  invisible(coord.save)

}


getDegree <- function ( K, upper=FALSE, diagonal=FALSE ) {

  K <- as.matrix( K )
  if ( upper )
    K[lower.tri(K)] <- 0
  K <- (K != 0) * 1
  diag(K) <- ( diag(K) !=0 ) * diagonal
  d.in  <- rowSums( K )
  d.out <- colSums( K )
  res <- (d.in + d.out) 
  return( res )
}

sub.param <- function (param, default=NULL, ... ) {

  if (missing(param)) { return(default) }

  l <- list(...)
  res <- l[[  which( names(l) %in% param)[1] ]] 

  if (is.null(res)) { res <- default }

  return ( res )

}

is.symetric <- function( X, eps=0 ) {
  
  if ( dim(X)[1] != dim(X)[2]) {
    # Matrix of edges
    n      <- max(X)
    A      <- matrix( 0, n, n )
    A[ t(X) ] <- 1
  } else {
    A <- X
  }
  
  sym <- all ( abs(A- t(A)) <= eps )
  return (sym )
  
}





Gplot.env <- new.env( parent=baseenv() )
assign("Gplot.graph",
       list( pie.gr=NULL, default.gr=NULL, reso=1000,
            col= sample( light.palette(10,black=FALSE) ) ),
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
Gplot <- function ( X, type="default", colors=NULL, ... ) {
  
  Env <- get("Gplot.graph", envir=Gplot.env)

  if( ! is.null( colors ) ) {
    Env$col <- colors
    assign("Gplot.graph", Env, envir=Gplot.env )
  }
  res <- -1

  if( type == "default" ) {
    dim.X      <- dim(X)
    dim.gr.mat <- dim( Env$default.gr$mat)
    if ( length( dim.X ) == length( dim.gr.mat ) ) {
      same.node.pos <- all ( dim( Env$default.gr$mat) ==  dim(X) )
    } else {
      same.node.pos <- FALSE
    }

    if( ( is.null(X) ) & !same.node.pos ) {
      # Clear memory of stored node positions
      Env$default.gr$network$coord=NULL
    } else {
      if( same.node.pos ) {
        # Keep the same node coordinates
        res <- Gplot.default(  X,  coord=Env$default.gr$network$coord, ...)
      } else {
        res <- Gplot.default(  X, ...)
      }
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
##    label    : vector of node labels (override dimnames(PrecMat))
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
  
  directed		<- sub.param("directed"		, NULL	, ...)

  class 		<- sub.param("class"		, NULL	, ...)

  #  Colours for each class
  class.col		<- sub.param("class.col"	, NULL	, ...)
  
  #  Coordinates for nodes
  coord			<- sub.param("coord"		, NULL	, ...)
  
  #  Labels for nodes                                       
  label		<- sub.param("label"	, NULL	, ...)
  if( is.null(label) )
    label.display=FALSE
  else
    label.display <- ( label == "default" )
  

  # Filter on node degree (threshold under which edges are not taking into account)
  degree.threshold	<- sub.param("degree.filter"	, 0,	
                                     ...)
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
  class.col <- Env$col
  
  #  Matrix check
  NNodes <- dim(X)[1]
  if (NNodes != dim(X)[2]) {
    cat("Gplot error : X matrix must be square \n")
    return(coord)
  }
  
  # Set "directed" parameter when not specified
  if ( is.null( directed ) ) {
    symetric <- (is.symmetric(X, eps=1.0e-3))
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
  levels <- sort( unique( class ) )
  nb.lev <- length( levels )
  # Isolated nodes don't belong to a class
  if( levels[1] == -1 )
     nb.lev <- nb.lev - 1
  if ( nb.lev > length( class.col )) {
    class.col <- sample( light.palette( nb.lev, black=FALSE ) )
  }

  vertex.col <- vector(length=length( class ) )
  vertex.col[( class  != -1)] <- class.col[ class[ (class != -1) ] ]
  # Unconnected nodes are "white" 
  vertex.col[ ( class  == -1)] <- "white"
  
  coord.save <- coord
  
  #    Process dustbin nodes
  #    ---------------------

  # Identify dustbin nodes
  # b.dust <- rep(FALSE, NNodes)
  # b.dust <- ( getDegree( X, upper=TRUE ) < degree.threshold )
  # dustbin.not.empty <- any( b.dust )

  if( NNodes == 0) {
    cat("Gplot warning : no node to display \n")
    return(coord)
  }
  
  # Remove dustbin nodes
  # if (dustbin.not.empty) {
  #  X <- X[!b.dust,!b.dust]
  #  NNodes <- dim(X)[1]
  #  if (!is.null(coord)) {
  #    coord <- coord [!b.dust,]
  #  }
  #  if (!is.null(class)) {
  #    class <- class [!b.dust]
  #  }
  #  if (!is.null(label)) {
  #    label <- label [!b.dust]
  #  }
  #  vertex.col <- vertex.col[!b.dust]
  # }

  #    Set graphic parameters
  #    ----------------------


  # Edge color handling
  e.col <- X
  e.col[X<0] <- "red"  # red : negative 
  e.col[X>0] <- "blue" # blue : positive 
  e.col[X==0]<- "black"	

  # Edge width
  lwd <- abs(X)
  
  if( length(lwd) > 100 ) {
    # Number of edges to high ( width set to 1) )
    lwd[] <- ( lwd[] > 0 ) * 1
  } else {
    # max line weight is 5
    lwd <- lwd * 5 / max( lwd )
  }

  # Edge type
  e.lty <- "solid"
  # e.lty[abs(lwd)>1]  <- "solid"
  # e.lty[abs(lwd)<=1] <- "dotted"
  # e.lty[abs(lwd)==0] <- "blank"

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
    gr <- Gplot.network( gr,
                         loop.draw       = directed,
                         displayisolates = TRUE,
                         coord           = coord
                       )
    if( ! no.edge ) 
      gr <- Gplot.edgeMIXER( gr,
                       col       = e.col,
                       lty       = e.lty,
                       lwd       = lwd,
                       curved    = (NNodes < 50) ,
                       arrow.draw = directed,
                       loop.draw  = directed,
                       loop.cex   = 10
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
  if( is.null(label) )
    label.display=FALSE
  else
    label.display <- ( label == "default" )
  
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
    symetric <- (is.symmetric(Mat, eps=1.e-6) )
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
     gr <- Gplot.edgeMIXER( gr,
                      col        = "blue",
                      lty        = e.lty,
                      lwd        = e.lwd,
                      curved     = (NNodes < 50) ,
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

## __________________________________________________________
##
## CalcDegrees
##
## Internal
##
## author : A. SMITH
##
## INPUT : 
##	K : precision matrix
##	MARGIN : change this to get in- and out- degrees
##		note : on a symmetric matrix, 
##		makes no difference
## OUTPUT :
##	vector of node degrees
## __________________________________________________________
##
CalcDegrees <- function ( K, MARGIN=1 ) {

	K <- as.matrix(K)
	diag(K) <- 0
	return( apply( abs(K), MARGIN, sum ) )
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

Symmetrize <- function (X) {

  n <- dim(X)[1]
  for (j in 1:n) X[j:n, j] <- X[j, j:n]

  return(X)
}

sub.param <- function (param, default=NULL, ... ) {

  if (missing(param)) { return(default) }

  l <- list(...)
  res <- l[[  which( names(l) %in% param)[1] ]] 

  if (is.null(res)) { res <- default }

  return ( res )

}

is.symmetric <- function( X, eps=0 ) {
  
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






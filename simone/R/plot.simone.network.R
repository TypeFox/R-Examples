is.simone.network <- function(x){
  if (class(x)=="simone.network") TRUE else FALSE
}

plot.simone.network<-function(x, y=NULL, type="default", last.coord=FALSE, ...){

  # Test arguments x, y
  if ( ! is.simone.network( x ) |
      ( ! is.null( y ) & (!is.simone.network( y ) ) ) ) {
      stop("x or y are simone.network object")
  }

  # Test if x and y are type class
  directed <- x$directed
  
  layout.class <- c( "circle", "circles", "cluster" )
  compare.class <- c( "overlap", "4graphs" )
  
  # Test x
  #  if (!is.simone( x ))
  #    stop("Not a simone object")

  #  Final layout (final.layout.class) :
  #  "graph.class"  graph class :
  #    "circle", "circles", "graph"
  #  "compare.class" :
  #    "single.compare.graph"
  #  four graphs to compare

  final.layout.class <- "graph.class"

  # According to y, select between "graph.class" or "compare.class"
  if( ! is.null( y ) ){
    if( all( dim(x$Theta) == dim(y$Theta)) ) {
      final.layout.class <- "compare.class"
    } else {
      stop("The two graph do not share the same dimension")
    }
  }

  if( type == "default" ){
    if( final.layout.class == "graph.class") {
      final.layout <- "cluster" 
    }
    else if( final.layout.class == "compare.class" ){
      final.layout <- "4graphs" 
    }
    else {
      stop("Not implemented")
    }
  } else {
    if( type %in% layout.class ) 
      l.c.n <-  "graph.class"
    else  if( type %in% compare.class ) 
      l.c.n <-  "compare.class"
    
    if( final.layout.class == l.c.n ) {
      final.layout <- type
      
    } else {
      stop("The arguments are not suitable with graphic type")
    }
    
  }
  
  if (final.layout == "cluster" ) {
    # Reduce margins
    old <- par(mar=c(1,1,1,1))
    coord <- "default"
    # Use old nodecoordinates
    if ( last.coord )  coord <- NULL
    Gplot( x$Theta, class= as.vector(x$clusters),
          label="default", coord=coord,
          directed=x$directed,
          main=x$name )
    par(mar=old$mar)
  } else if (final.layout == "circle" ) {
    # Reduce margins
    old <- par(mar=c(1,1,1,1))
    coord <- "circle"
    # Use old nodecoordinates
    if ( last.coord )  coord <- NULL
    Gplot( x$Theta, class= as.vector(x$clusters),
          directed=x$directed,
          label="default", coord=coord,
          main=x$name, ... )
    par(mar=old$mar)
  } else if( (final.layout == "circles") ){
    old <- par(mar=c(1,1,1,1))
    coord <- "circles"
    # Use old nodecoordinates
    if ( last.coord )  coord <- NULL
    Gplot( x$Theta, class= as.vector(x$clusters),
          label="default", coord=coord,
          directed=x$directed,
          main=x$name, ... )
    par(mar=old$mar)
  } else if( (final.layout == "4graphs") ){

    old <- par(mar=c(1,1,1,1), mfrow=c(2, 2))
    coord <- "default"
    # Use old nodecoordinates
    if ( last.coord )  coord <- NULL
    Gplot( x$Theta, class= as.vector(x$clusters),
            label="default",
            coord=coord,
            directed=x$directed,
            main=x$name, ... )

    Gplot( y$Theta, class= as.vector(y$clusters),
            label="default",
            coord=NULL,
            directed=x$directed,
            main=y$name, ... )

      xbin <- (x$Theta != 0 )
      ybin <- (y$Theta != 0 )
      
      Gplot( xbin & ybin, class= as.vector(x$clusters),
            label="default",
            coord=NULL,
            directed=x$directed,
            main="Intersection", ... )

      Gplot( xbin != ybin, class= as.vector(x$clusters),
            label="default",
            coord=NULL,
            directed=x$directed,
            main="Difference", ... )
      
    par(mar=old$mar, mfrow=old$mfrow)
    
  } else if( (final.layout == "overlap") ){

    old <- par(mar=c(1,1,1,1))

    xcmp <- ( x$Theta != 0 ) * 1
    xcmp[ which(y$Theta != 0 ) ] <- -1
    xcmp[ which((y$Theta != 0) & (x$Theta != 0)) ] <- 2
    main=paste("Overlapping between", x$name, "and", y$name)
    Color.Map <- function( X ) {
      e.col <- X
      e.col[X==1] <- "blue"
      e.col[X==-1] <- "red"  # red : negative 
      e.col[X==2] <- "grey"   # blue : positive
      return( e.col )
    }
    coord <- "default"
    # Use old nodecoordinates
    if ( last.coord )  coord <- NULL

    Gplot( xcmp, class= as.vector(x$clusters),
            label="default",
            coord=coord,
            color.map = Color.Map,
            directed=x$directed,
            main=main, ... )

    par(mar=old$mar)

  } else {
    stop("Not implemented")
  }
}


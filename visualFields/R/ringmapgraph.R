ringmapgraph <- function( ncol = 3, mapval = NULL, txtfont = "mono", pointsize = 7,
                          outerSymbol = "circles", innerSymbol = "circles",
                          outerSize = 1, innerSize = 1,
                          outerInch = 0.2, innerInch = 0.1,
                          outerBorderThickness = 2, innerBorderThickness = 2 ) {
  if( is.null( mapval ) ) {
    mapval             <- NULL
    mapval$cutoffs     <- c( "0.5", "1", "5" )
    mapval$innerCircle <- c( 1, 0, 1 )
    mapval$outerCircle <- c( 1, 1, 0 )
    mapval             <- as.data.frame( mapval )
  }
  total <- nrow( mapval )
  nrow <- ceil( total / ncol )

  # get coordinates to plot
  coords   <- NULL
  coords$x <- ( c( 1:total ) - 1 ) %% ncol + 1
  coords$y <- ( c( 1:total ) ) %% nrow + 1

  coords   <- as.data.frame( coords )
  coords   <- coords[order( coords$x ),]
  coords   <- coords[order( coords$y, decreasing = TRUE ),]

  xmin     <- min( coords$x ) - 1 / ncol
  xmax     <- max( coords$x ) + 1 / ncol
  ymin     <- min( coords$y ) - 1 / nrow
  ymax     <- max( coords$y ) + 1 / nrow
# get rgb and text to plot
  innerBorderColor          <- as.data.frame( matrix( rep( 0.5, total * 3 ), total ) )
  names( innerBorderColor ) <- c( "red", "green", "blue" )
  outerBorderColor          <- innerBorderColor
  txtval                    <- NULL
  for( i in 1:nrow( mapval ) ) {
    if( mapval$innerCircle[i] == 1 ) innerBorderColor[i,] <- 0
    if( mapval$outerCircle[i] == 1 ) outerBorderColor[i,] <- 0
    txtval[i] <- as.character( mapval$cutoffs[i] )
  }
  innerBorderColor <- rgb( innerBorderColor )
  outerBorderColor <- rgb( outerBorderColor )

# opar <- par( no.readonly = TRUE )
  oplt    <- par()$plt
  ops     <- par()$ps
  ofamily <- par()$family
  par( plt = c( 0, 1, 0, 1 ) )
  par( ps = pointsize )
  par( family = txtfont )

# legend
  plot( coords$x, coords$y, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c( xmin, xmax ), ylim = c( ymin, ymax ) )
  outerDimensions <- t( matrix( data = rep( outerSize, nrow( coords ) ),nrow = length( outerSize ), ncol = nrow( coords ) ) )
  evaltxt <- paste( "symbols( coords$x, coords$y, " , outerSymbol, " = outerDimensions, add = TRUE, inches = outerInch, bg = 'white', fg = outerBorderColor, lwd = 1 )", sep = "" )
  eval( parse( text = evaltxt ) )
  innerDimensions <- t( matrix( data = rep( innerSize,nrow( coords ) ),nrow = length( innerSize ), ncol = nrow( coords ) ) )
  evaltxt <- paste( "symbols( coords$x, coords$y, ", innerSymbol," = innerDimensions, add = TRUE, inches = innerInch, bg = 'white', lwd = 1, fg = innerBorderColor )", sep = "" )
  eval( parse( text = evaltxt ) )

  text( coords$x, coords$y, labels = txtval, adj = 0.525 )

  par( plt = oplt )
  par( ps = ops )
  par( family = ofamily )
}

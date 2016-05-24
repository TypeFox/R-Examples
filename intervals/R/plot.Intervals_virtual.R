plot.Intervals_full <- function(
                                x, y = NULL,
                                axes = TRUE,
                                xlab = "", ylab = "",
                                xlim = NULL, ylim = NULL,
                                col = "black", lwd = 1,
                                cex = 1,
                                use_points = TRUE,
                                use_names = TRUE,
                                names_cex = 1,
                                ...
                                )
{
  # Subset first, for efficiency, and so that the maximal y plotted is
  # appropriate for the region shown.

  if ( any( is.na( x ) ) ) x <- x[ !is.na(x), ]
  
  if ( is.null(xlim) )
    xlim <- range( x@.Data )
  else
    x <- x[ x[,2] >= xlim[1] & x[,1] <= xlim[2], ]
  
  if ( is.null(y) )
    y <- .Call( "_plot_overlap", x@.Data, closed(x), is( x, "Intervals_full" ) )

  if ( is.null(ylim) )
    ylim <- c( 0, max( y ) )

  plot(
       0, 0,
       type = "n",
       xlim = xlim, ylim = ylim,
       axes = FALSE,
       xlab = xlab, ylab = ylab,
       ...
       )
  # Careful with non-finite endpoints, which segments() ignores.
  segments(
           pmax( x[,1], par("usr")[1] ), y,
           pmin( x[,2], par("usr")[2] ), y,
           col = col,
           lwd = lwd
           )
  if ( use_points ) {   
    # Careful with points...
    adjust <- ( x[,1] == x[,2] ) & !closed(x)[,1]
    closed(x)[ adjust, 2 ] <- FALSE
    points(
           x@.Data, rep( y, 2 ),
           pch = 21, cex = cex,
           col = col, bg = ifelse( closed(x), col, "white" )           
           )
  }
  if ( use_names && !is.null( rownames(x) ) ) {
    mids <- ( x[,1] + x[,2] ) / 2
    text(
         mids, y,
         rownames( x ),
         pos = 3, offset = .5,
         cex = names_cex,
         xpd = NA
         )         
  }
  if ( axes )
    axis( 1 )
}  

plot.Intervals <- function( x, y = NULL, ... ) {
  plot( as( x, "Intervals_full" ), y, ... )
}

setMethod( "plot", c( "Intervals", "missing" ), function( x, y, ... ) plot.Intervals( x, ... ) )
setMethod( "plot", c( "Intervals", "ANY" ), function( x, y, ... ) plot.Intervals( x, y, ... ) )

setMethod( "plot", c( "Intervals_full", "missing" ), function( x, y, ... ) plot.Intervals_full( x, ... ) )
setMethod( "plot", c( "Intervals_full", "ANY" ),  function( x, y, ... ) plot.Intervals_full( x, y, ... ) )


progols <- function( tdate, index, projyears = 5,
                     xlab = "years from first visit", ylab = "md",
                     txtfont = "mono", pointsize = 12, cex = 1 ) {

  lt      <- as.POSIXlt( tdate )
  mon     <- lt$year * 12 + lt$mon
  yeardif <- ( mon - mon[1] ) / 12
  xreg    <- c( min( yeardif ), max( yeardif ) + projyears )
  xlim    <- NULL
  ylim    <- NULL
  xlim[1] <- xreg[1]
  xlim[2] <- xreg[2] + 1
  ylim <- c( max( index ) - 11, max( index ) + 1 )
# get regression
  mdreg <- lm( index ~ yeardif )
  pval  <- summary( mdreg )$coefficients[2,4]
  pval  <- round( pval * 1000 ) / 1000
  ylab  <- paste( ylab, ", p = ", as.character( pval ), sep = "" )
  yreg  <- mdreg$coefficients[1] + mdreg$coefficients[2] * xreg

  ops     <- par()$ps
  ofamily <- par()$family
  oplt    <- par()$plt
  par( ps     = pointsize )
  par( family = txtfont )
  par( plt    = c( 0.29, 1.0, 0.35, 0.9 ) )
  
  firstTick <- round( max( index ) - 11 )
  tickMarks <- c( firstTick, firstTick + 4, firstTick + 8, firstTick + 12 )
  
  plot( yeardif, index, type = "n", axes = FALSE, ann = FALSE, xlim = xlim, ylim = ylim )
  axis( 1, las = 1, tcl = -.3, lwd = 0.5, lwd.ticks = 0.5 )
  axis( 2, las = 1, tcl = -.3, lwd = 0.5, lwd.ticks = 0.5, at = tickMarks )
  grid( nx = NA, ny = NULL, lty = "solid", "lightgray" )
  points( yeardif, index )
  lines( xreg, yreg )
  box()
  title( xlab = xlab, mgp = c( 2, 1, 0 ) )
  title( ylab = ylab, mgp = c( 2.7, 1, 0 ) )

  par( new    = FALSE )
  par( ps     = ops )
  par( family = ofamily )
  par( plt    = oplt )
}

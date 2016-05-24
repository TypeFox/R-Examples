bebie <- function( tdr, type = "conventional", diff = TRUE, percentiles = TRUE,
                   correction = TRUE, txtfont = "mono", pointsize = 7, cex = 1 ) {

  
  if( type != "conventional" & type != "ghrank" ) stop( "wrong type of TD-rank curve plot" )

  xlab <- "rank"
  linewdt <- 1.5
  if( diff ) {
    ylab <- "dB difference"
  } else {
    ylab <- "dB"
  }

  if( type == "conventional" ) {
    if( nrow( tdr ) > 1 ) {
      stop( "Error! only one visual field to use here" )
    }
    tdrval    <- as.numeric( tdr[,visualFields::vfsettings$locini:ncol( tdr )] )
    rank      <- c( 1:length( tdrval ) )
# get reference
    evaltxt        <- paste( "vfenv$nv$", tdr$tpattern, "_", tdr$talgorithm, "$nvtdrank", sep = "" )
    tdrvalref      <- as.numeric( eval( parse( text = evaltxt ) )$mtdr )
    tdrvalsubtract <- tdrvalref
    gh             <- ghpostd( tdr, correction = TRUE )
    evaltxt        <- paste( "vfsettings$", tdr$tpattern, sep = "" )
    settings       <- eval( parse( text = evaltxt ) )
  } else if ( type == "ghrank" ) {
    rank           <- tdr$rank
    tdrval         <- tdr$tdrn
    evaltxt        <- paste( "vfenv$nv$", tdr$tdr$tpattern, "_", tdr$tdr$talgorithm, "$nvtdrank", sep = "" )
    tdrvalref      <- as.numeric( eval( parse( text = evaltxt ) )$mtdr )
    tdrvalsubtract <- tdr$mtdr
    gh             <- tdr$gh
    evaltxt        <- paste( "vfsettings$", tdr$tdr$tpattern, sep = "" )
    settings       <- eval( parse( text = evaltxt ) )
  }
  if( diff ) ylim <- c( -21, 5 ) else ylim <- c( -25, 10 )
# set limits
  xlim <- c( 1, settings$locnum - length( settings$bs ) )

# get differences
  if( diff ) tdrval <- tdrval - tdrvalsubtract
# get correction
  if( correction ) tdrvalc <- tdrval + gh
# get percentiles
  if( percentiles ) {
    if( diff ) {
      if( type == "conventional" ) {
        evaltxt <- paste( "vfenv$nv$", tdr$tpattern, "_", tdr$talgorithm, "$perctdrankadj7", sep = "" )
      } else{
        evaltxt <- paste( "vfenv$nv$", tdr$tdr$tpattern, "_", tdr$tdr$talgorithm, "$perctdrankadjghr", sep = "" )
      }
    } else {
      if( type == "conventional" ) {
        evaltxt <- paste( "vfenv$nv$", tdr$tpattern, "_", tdr$talgorithm, "$perctdrank", sep = "" )
      } else{
        evaltxt <- paste( "vfenv$nv$", tdr$tdr$tpattern, "_", tdr$tdr$talgorithm, "$perctdrank", sep = "" )
      }
    }
    tdrperc <- eval( parse( text = evaltxt ) )
  }
  ops     <- par()$ps
  ofamily <- par()$family
  par( ps     = pointsize )
  par( family = txtfont )

  if( diff ) {
    plot(c( xlim[1], xlim[2] ), c( 0, 0 ), axes = FALSE, ann = FALSE, xlim = xlim, ylim = ylim, type = "l" )
  } else {
    plot( tdrvalref, axes = FALSE, ann = FALSE, xlim = xlim, ylim = ylim, type = "l", lwd = linewdt )
  }

  axis( 1, las = 1, tcl = -.3, lwd = 0.5, lwd.ticks = 0.5 )
  axis( 2, las = 1, tcl = -.3, lwd = 0.5, lwd.ticks = 0.5 )
  grid( nx = NA, ny = NULL, lty = "solid", "gray" )
  box()
  title( xlab = xlab, mgp = c( 2, 1, 0 ) )
  title( ylab = ylab, mgp = c( 2.3, 1, 0 ) )

  if( percentiles ){
    for( i in 1:( ncol( tdrperc ) - 1 ) ) {
      lines( tdrperc[,i], col = rgb( red = visualFields::vfenv$nv$pmapsettings$red[i], green = visualFields::vfenv$nv$pmapsettings$green[i], blue = visualFields::vfenv$nv$pmapsettings$blue[i] ), lwd = linewdt )
    }
    lines( tdrperc[,ncol( tdrperc )], col = rgb( red = visualFields::vfenv$nv$pmapsettings$red[nrow( visualFields::vfenv$nv$pmapsettings )], green = visualFields::vfenv$nv$pmapsettings$green[nrow( visualFields::vfenv$nv$pmapsettings )], blue = visualFields::vfenv$nv$pmapsettings$blue[nrow( visualFields::vfenv$nv$pmapsettings )] ), lwd = linewdt )
  }

  points( rank, tdrval, xlim = xlim, ylim = ylim, pch = 1, cex = cex )
  if( correction ) points( rank, tdrvalc, xlim = xlim, ylim = ylim, pch = 16, cex = cex )

  par( new    = FALSE )
  par( ps     = ops )
  par( family = ofamily )

}
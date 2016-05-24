vflayoutghr <- function( vf, pwidth = 8.27, pheight = 11.69, margin = 0.25,
                      filename = NULL ) {

  if( nrow( vf ) > 1 ) {
    stop("Error! vf cannot have more than 1 rows")
  }
  
  ffamily        <- "serif"
  sizetxt        <- 12
  sizetxtSmall   <- 8
  ffamilyvf      <- "serif"
  pointsize      <- 7
  outerSymbol    <- "circle"
  outerInch      <- 0.13
  innerSymbol    <- "circle"
  innerInch      <- outerInch / 1.9
  inch2axisunits <- 12.528
  lengthLines    <- 1.35 * 2 * outerInch * inch2axisunits
  thicknessLines <- 1.5

# open window wiht A4 page
  if( is.null( filename ) ) {
    device <- options( "device" )
    if( .Platform$OS.type == "unix" ) {
      if( Sys.info()["sysname"] == "Darwin" ) {
        options( device = "quartz" )
        dev.new( width = pwidth, height = pheight, dpi = 85 )
      } else {
        options( device = "x11" )
        dev.new( width = pwidth, height = pheight )
      }
    } else{
      options( device = "windows" )
      dev.new( width = pwidth, height = pheight, rescale = "fixed" )
    }
    options( device = device )
  } else {
    pdf( width = pwidth, height = pheight, file = filename )
  }

# define the margins
  mwidth  <- pwidth  - 2 * margin
  mheight <- pheight - 2 * margin

# create the layout of the printout
  printout <- createviewport( "printout", left = margin, top = margin, height = mheight, width = mwidth )

######################################################
# first plot all graphs
######################################################
  if( vf$tpattern == "p24d2" ) {
    xminmax <- 29
    yminmax <- 29
  } else if( vf$tpattern == "p30d2" ) {
    xminmax <- 29
    yminmax <- 29    
  } else if( vf$tpattern == "p10d2" ) {
    xminmax <- 10
    yminmax <- 10
  } else if( vf$tpattern == "sgrnfl" ) {
    xminmax <- 29
    yminmax <- 29
  } else {
    xminmax <- 100
    yminmax <- 100
  }
# total-deviation plot
  opar <- par( no.readonly = TRUE )
  par( fig = c( 0.3869, 0.9915, 0.5483, 0.9760 ) )
  vfplot( vf, plotType = "td", txtfont = ffamilyvf, pointsize = pointsize,
          xminmax = xminmax, yminmax = yminmax,
          outerSymbol = outerSymbol, innerSymbol = innerSymbol,
          outerInch = outerInch, innerInch = innerInch,
          lengthLines = lengthLines, thicknessLines = thicknessLines )
# sensitivity plot
  par( new = TRUE )
  par( fig = c( 0.0109, 0.6155, 0.2746, 0.7023 ) )
  vfplot( vf, plotType = "vf", txtfont = ffamilyvf, pointsize = pointsize,
          xminmax = xminmax, yminmax = yminmax,
          outerSymbol = outerSymbol, innerSymbol = innerSymbol,
          outerInch = outerInch, innerInch = innerInch,
          lengthLines = lengthLines, thicknessLines = thicknessLines )
# pattern-deviation plot
  par( new = TRUE )
  par( fig = c( 0.3869, 0.9915, 0.0060, 0.4337 ) )
  vfplot( vf, plotType = "pdghr", txtfont = ffamilyvf, pointsize = pointsize,
          xminmax = xminmax, yminmax = yminmax,
          outerSymbol = outerSymbol, innerSymbol = innerSymbol,
          outerInch = outerInch, innerInch = innerInch,
          lengthLines = lengthLines, thicknessLines = thicknessLines )
# stimulus locations
  par( new = TRUE )
  par( fig = c( 0.769, 0.981, 0.42, 0.57 ) )
  par( mar = c( 0, 0, 0.5, 0.5 ) )
  stimLoc( perimetry = vf$tperimetry, pattern = vf$tpattern, eye = vf$seye,
          txtfont = ffamilyvf, pointsize = pointsize,
          xminmax = xminmax, yminmax = yminmax )
# Bebie difference curve
  par( new = TRUE )
  par( fig = c( 0.007, 0.4039, 0.05, 0.2891 ) )
  par( mar = c( 3.25, 4.2, 0.5, 0.5 ) )
  tdr <- ghranktd( tdval( vf ) )
  bebie( tdr, type = "ghrank", txtfont = ffamily, pointsize = sizetxt, cex = 0.75 )
# color-code map
  par( new = TRUE )
  par( fig = c( 0.03, 0.3869, 0.015, 0.060 ) )
  colormapgraph( ncol = 6, txtfont = ffamilyvf, pointsize = pointsize, outerSymbol = outerSymbol, innerSymbol = innerSymbol,
                 outerInch = outerInch, innerInch = innerInch )
  par( opar )

######################################################
# create the text elements in the printouts
######################################################
# The two above are to delete once the graphs are generated!!!
  mainInfo  <- createviewport( "mainInfo",  left =  0.00, top =  0.00, width = 4.75, height = 0.40, pheight = mheight, pwidth = mwidth )
  tdtext    <- createviewport( "tdtext",    left =  4.54, top =  0.00, width = 1.83, height = 0.22, pheight = mheight, pwidth = mwidth )
  infobox2  <- createviewport( "infobox2",  left =  6.37, top =  0.00, width = 1.40, height = 0.40, pheight = mheight, pwidth = mwidth )
  infobox3  <- createviewport( "infobox3",  left =  6.37, top = 10.89, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  infobox1  <- createviewport( "infobox1",  left =  0.25, top =  0.50, width = 2.95, height = 0.40, pheight = mheight, pwidth = mwidth )
  infotest1 <- createviewport( "infotest1", left =  0.25, top =  1.10, width = 1.20, height = 0.65, pheight = mheight, pwidth = mwidth )
  infotest2 <- createviewport( "infotest2", left =  1.55, top =  1.10, width = 0.60, height = 0.65, pheight = mheight, pwidth = mwidth )
  results1  <- createviewport( "results1",  left =  0.25, top =  1.95, width = 0.50, height = 1.10, pheight = mheight, pwidth = mwidth )
  results2  <- createviewport( "results2",  left =  0.85, top =  1.95, width = 0.70, height = 1.10, pheight = mheight, pwidth = mwidth )
  results3  <- createviewport( "results3",  left =  1.60, top =  1.95, width = 1.00, height = 1.10, pheight = mheight, pwidth = mwidth )
  pdtext    <- createviewport( "pdtext",    left =  4.54, top =  6.34, width = 1.83, height = 0.22, pheight = mheight, pwidth = mwidth )

# create the list and then generate the tree and "push" it
  list <- vpList( mainInfo, infobox1, infotest1, infotest2, results1, results2, results3, infobox2, infobox3, tdtext, pdtext )
  tree <- vpTree( printout, list )

  pushViewport( tree )

#  seekViewport( "printout" )
#  grid.rect( gp = gpar( col = "blue" ) )

######################################################
# perimetry information
######################################################
  seekViewport( "mainInfo" )
  text <- vf$tperimetry
  if( text == "sap" ) {
    text = "Standard Automatic Perimetry."
  } else if( text == "fdp" ) {
    text = "Frequency-doubling Perimetry."
  } else if( text == "csp" ) {
    text = "Contrast Sensitivity Perimetry."
  }
  text <- paste( text, "Single field analysis", sep = " " )
# ID
  text <- paste( text, "Subject ID: ", sep = "\n" )
  text <- paste( text, vf$id, ",", sep = "" )
# age
  text <- paste( text, " age: ", round( vf$sage ), ",", sep = "" )
  # eye
  texteye <- paste( "eye:", vf$seye, sep = " " )
  if( vf$seye == "OD" ) {
    texteye <- paste( texteye, "(right)", sep = " " )
  } else if ( vf$seye == "OS" ) {
    texteye <- paste( texteye, "(left)", sep = " " )
  } else {
    texteye <- paste( texteye, "(which?)", sep = " " )
  }
  text <- paste( text, texteye, sep = " " )
  grid.text( text, x = 0.0, y = 1.0, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt, fontface = "bold" ) )

######################################################
# text: total deviation and pattern deviation
######################################################
  seekViewport( "tdtext" )
  grid.text( "total deviation", x = 0.5, y = 1.0, just = c( "center", "top" ), gp = gpar( fontfamily = ffamily,fontsize = sizetxt ) )
  seekViewport( "pdtext" )
  grid.text( "pattern deviation", x = 0.5, y = 1.0, just = c( "center", "top" ), gp = gpar( fontfamily = ffamily,fontsize = sizetxt ) )

######################################################
# Details about printouts
######################################################
  seekViewport( "infobox2" )
  
  if( vf$tpattern == "p24d2" ) {
    textpattern <- "Central 24-2"
  } else if( vf$tpattern == "p30d2" ) {
    textpattern <- "Central 30-2"
  } else if( vf$tpattern == "p10d2" ) {
    textpattern <- "Central 10-2"
  } else if( vf$tpattern == "sgrnfl" ) {
    textpattern <- "CSP-SG-RNFL-57"
  } else {
    textpattern <- "Unknown"
  }
  # algorithm
  if( vf$talgorithm == "sitas" ) {
    textalgorithm <- "SITA standard"
  } else if( vf$talgorithm == "sitaf" ) {
    textalgorithm <- "SITA fast"
  } else if( vf$talgorithm == "fullt" ) {
    textalgorithm <- "Full threshold"
  } else if( vf$talgorithm == "zest" ) {
    textalgorithm <- "ZEST"
  } else {
    textalgorithm <- "Unknown"
  }
  
  text <- paste( textpattern, textalgorithm, sep = "\n" )
  grid.text( text, x = 1.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

  ######################################################
  # Details about printouts
  ######################################################
  seekViewport( "infobox3" )
  
  text <- paste( "norm vals: ", visualFields::vfenv$nv$nvname, sep = "" )
  text <- paste( text, substr( packageDescription( "visualFields" )$Date, 1, 4 ), sep = "\n" )
  text <- paste( text, "visualFields", packageDescription( "visualFields" )$Version, sep = " " )
  grid.text( text, x = 1.00, y = 0.00, just = c( "right", "bottom" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxtSmall ) )

######################################################
# subject and test information
######################################################
  seekViewport( "infobox1" )

  timetxt <- substr( vf$ttime, 1, 5 )
  if( substr( timetxt, 1, 1 ) == "0" ) substr( timetxt, 1, 1 ) <- ""
  text <- paste( "Date:", format( vf$tdate, "%m/%d/%Y" ), "at", timetxt, sep = " " )
# duration and pause of test
  timetxt         <- substr( vf$sduration, 3, nchar( vf$sduration ) )
  if( substr( timetxt, 1, 1 ) == "0" ) substr( timetxt, 1, 1 ) <- ""
  text <- paste( text, paste( "Duration: ", timetxt, sep = " " ), sep = "\n" )
  timetxt         <- substr( vf$spause, 3, nchar( vf$sduration ) )
  if( timetxt != "59:59" ) {
    if( substr( timetxt, 1, 1 ) == "0" ) substr( timetxt, 1, 1 ) <- ""
    text <- paste( text, paste( ", pause: ", timetxt, sep = "" ), sep = "" )
  }
  text <- paste( text, "", sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
######################################################
# add false positives and negatives, fixation losses
######################################################
  seekViewport( "infotest1" )

  text <- "fixation losses"
  text <- paste( text, "false positives", sep = "\n" )
  text <- paste( text, "false negatives", sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

  seekViewport( "infotest2" )

  sfp <- paste( sprintf( "%.1f", round( 1000 * vf$sfp ) / 10 ), "%", sep = " " )
  sfn <- paste( sprintf( "%.1f", round( 1000 * vf$sfn ) / 10 ), "%", sep = " " )
  sfl <- paste( sprintf( "%.1f", round( 1000 * vf$sfl ) / 10 ), "%", sep = " " )

  text <- sfp
  text <- paste( text, sfn, sep = "\n" )
  text <- paste( text, sfl, sep = "\n" )
  grid.text( text, x = 1.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

######################################################
# visual-field results
######################################################
  vfs  <- vfstats( vf )
  vfi  <- vfindex( vf )
  vfsp <- vfstatspmap( vfs )
  vfip <- vfindexpmap( vfi )
# general-height difference, if the used normative values have one.
  texteval <- paste( "vfenv$nv$", vf$tpattern, "_", vf$talgorithm, "$nvtdrank$mtdr", sep = "" )
  tdr <- NULL
  tdr <- eval( parse( text = texteval ) )
  if( !is.null( tdr ) ) {
    ghd <- ghranktd( tdval( vf ) )$gh
    ghd <- paste( sprintf( "%.1f", round( 10 * ghd ) / 10 ), "dB", sep = " " )
  }

  ms  <- paste( sprintf( "%.1f", round( 10 * vfs$msens ) / 10 ), "dB", sep = " " )
  md  <- paste( sprintf( "%.1f", round( 10 * vfs$mtdev ) / 10 ), "dB", sep = " " )
  psd <- paste( sprintf( "%.1f", round( 10 * vfs$spdev ) / 10 ), "dB", sep = " " )
  vfi <- paste( sprintf( "%.1f", round( 10 * vfi$mvfi  ) / 10 ), " %", sep = " " )

  seekViewport( "results1" )

  text <- "MS"
  text <- paste( text, "MD", sep = "\n" )
  text <- paste( text, "PSD", sep = "\n" )
  text <- paste( text, "VFI", sep = "\n" )
  if( !is.null( tdr ) ) {
    text <- paste( text, "GHD", sep = "\n" )
  }
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

  seekViewport( "results2" )

  text <- ms
  text <- paste( text, md, sep = "\n" )
  text <- paste( text, psd, sep = "\n" )
  text <- paste( text, vfi, sep = "\n" )
  if( !is.null( tdr ) ) {
    text <- paste( text, ghd, sep = "\n" )
  }
  grid.text( text, x = 1.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

  seekViewport( "results3" )

  text <- ""
  textp <- paste( "(p < ", vfsp$mtdev, " %)", sep = "" )
  text <- paste( text, textp, sep = "\n" )
  textp <- paste( "(p < ", vfsp$mtdev, " %)", sep = "" )
  text <- paste( text, textp, sep = "\n" )
  textp <- paste( "(p < ", vfip$mvfi, " %)", sep = "" )
  text <- paste( text, textp, sep = "\n" )
  grid.text( text, x = 0.00, y = 1.00, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

# only if in save mode, then set device to off
  if( !is.null( filename ) ) {
    dev.off()
  }
}
vflayout_legoplot <- function( vf, grp = 3, pwidth = 8.27, pheight = 11.69,
                               margin = 0.25, filename = NULL,
                               ffamily = "serif", sizetxt = 12,
                               sizetxtSmall = 8,
                               ffamilyvf = "serif", pointsize = 6,
                               txtcolorlego = "red", pointsizelego = 16,
                               outerSymbol = "circle", outerInch = 0.12,
                               innerSymbol = "circle", innerInch = outerInch / 1.9,
                               inch2axisunits = 12.528,
                               lengthLines = 0, thicknessLines = 0,
                               outerSymbollego = "square", outerInchlego = 0.36,
                               innerSymbollego = "circle", innerInchlego = 0.16 ) {
##############
# input checks
##############
# check that all rows in vf belong to the same subject, the same test, the same perimetry
# testing and the same algorithm, and the same eye
  if( length( unique( vf$tperimetry ) ) > 1 |
      length( unique( vf$tpattern   ) ) > 1 |
      length( unique( vf$talgorithm ) ) > 1 |
      length( unique( vf$id ) ) > 1         |
      length( unique( vf$seye ) ) > 1 ) {
    stop( "all visual fields should belong to the same subject tested with the same perimeter and algorithm on the same locations" )
  }
  if( nrow( vf ) < 2 * grp ) stop( "the number of visual fields needs to be at least twice the number of visual fields to group for the display" )

# get settings for the pattern of test locations
  texteval <- paste( "vfsettings$", vf$tpattern[1], sep = "" )
  settings <- eval( parse( text = texteval ) )
# get x and y locations
  texteval <- paste( vf$tperimetry[1], "locmap$",  vf$tpattern[1], sep = "" )
  locmap   <- eval( parse( text = texteval ) )
# get TD and PD values
  td  <- tdval( vf )
# remove blind spot
  vf     <- vf[,-( settings$bs + visualFields::vfsettings$locini - 1 )]
  td     <- td[,-( settings$bs + visualFields::vfsettings$locini - 1 )]
  locmap <- locmap[-settings$bs,]
# add color for the text of the legoplots
txtcolorlego <- t( matrix( rep( col2rgb( txtcolorlego ) / 255, nrow( locmap ) ), 3, nrow( locmap ) ) )
txtcolorlego <- as.data.frame( txtcolorlego )
names( txtcolorlego ) <- c( "red", "green", "blue" )
######################
# Analysis
######################
# init
  vfinfo0 <- vf[1,1:( visualFields::vfsettings$locini - 1 )]
  vfinfo1 <- vf[1,1:( visualFields::vfsettings$locini - 1 )]
# get indices for averages
  locvalsidx <- visualFields::vfsettings$locini:( visualFields::vfsettings$locini + settings$locnum - length( settings$bs ) - 1 )
  idx0 <- c( 1:grp )
  idx1 <- c( ( nrow( vf ) - grp + 1 ):nrow( vf ) )
# get averages
  vfinfo0$sage <- mean( vf$sage[idx0] )
  vfinfo0$sbsx <- mean( vf$sbsx[idx0] )
  vfinfo0$sbsy <- mean( vf$sbsy[idx0] )
  vfinfo1$sage <- mean( vf$sage[idx1] )
  vfinfo1$sbsx <- mean( vf$sbsx[idx1] )
  vfinfo1$sbsy <- mean( vf$sbsy[idx1] )
  vf0          <- round( colMeans( vf[idx0, locvalsidx] ) )
  vf1          <- round( colMeans( vf[idx1, locvalsidx] ) )
  td0          <- colMeans( td[idx0, locvalsidx] )
  td1          <- colMeans( td[idx1, locvalsidx] )

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
  if( vfinfo0$tpattern == "p24d2" ) {
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
  opar <- par( no.readonly = TRUE )
# legoplot
  color0 <- vfgrayscale( vf0, vfinfo0$sage, pattern = vfinfo0$tpattern, algorithm = vfinfo0$talgorithm )
  color1 <- vfgrayscale( vf1, vfinfo1$sage, pattern = vfinfo1$tpattern, algorithm = vfinfo1$talgorithm )
  par( fig = c( 0.5000, 0.985, 0.5833, 0.9200 ) )
  vfplotloc( vf1 - vf0, eye = vfinfo0$seye, patternMap = locmap, outerColor = color0, innerColor = color1, axesCol = "white",
             txtfont = ffamilyvf, pointsize = pointsizelego, txtcolor = txtcolorlego,
             xminmax = xminmax, yminmax = yminmax,
             outerSymbol = outerSymbollego, innerSymbol = innerSymbollego,
             outerInch = outerInchlego, innerInch = innerInchlego,
             lengthLines = lengthLines, thicknessLines = thicknessLines )
# sensitivity plot first n visits
  par( new = TRUE )
  par( fig = c( 0.0150, 0.5000, 0.2332,  0.5700 ) )
  color <- vfgrayscale( vf0, vfinfo0$sage, pattern =  vfinfo0$tpattern, algorithm = vfinfo0$talgorithm )
  vf0[which( vf0 < 0 )] <- "<0"
  vfplotloc( vf0, eye = vfinfo0$seye, patternMap = locmap , outerColor = color, bs = c( vfinfo0$sbsx, vfinfo0$sbsy ), 
             txtfont = ffamilyvf, pointsize = pointsize,
             xminmax = xminmax, yminmax = yminmax,
             outerSymbol = outerSymbol, innerSymbol = innerSymbol,
             outerInch = outerInch, innerInch = innerInch,
             lengthLines = lengthLines, thicknessLines = thicknessLines )
  # sensitivity plot last n visits
  par( new = TRUE )
  par( fig = c( 0.5000, 0.985, 0.2332,  0.5700 ) )
  color <- vfgrayscale( vf1, vfinfo0$sage, pattern =  vfinfo0$tpattern, algorithm = vfinfo0$talgorithm )
  vf1[which( vf1 < 0 ) ] <- "<0"
  vfplotloc( vf1, eye = vfinfo1$seye, patternMap = locmap , outerColor = color, bs = c( vfinfo0$sbsx, vfinfo0$sbsy ), 
             txtfont = ffamilyvf, pointsize = pointsize,
             xminmax = xminmax, yminmax = yminmax,
             outerSymbol = outerSymbol, innerSymbol = innerSymbol,
             outerInch = outerInch, innerInch = innerInch,
             lengthLines = lengthLines, thicknessLines = thicknessLines )
  # Bebie first n visit
  par( new = TRUE )
  par( fig = c( 0, 0.4, 0, 0.25 ) )
  par( mar = c( 4, 4.4, 0.5, 0.5 ) )
  bebie(  tdrank( as.data.frame( c( vfinfo0, td0 ) ) ), correction = FALSE, txtfont = ffamily, pointsize = sizetxt, cex = 0.75 )
  # Bebie last n visits
  par( new = TRUE )
  par( fig = c( 0.5, 0.9, 0, 0.25 ) )
  par( mar = c( 4, 4.4, 0.5, 0.5 ) )
  bebie( tdrank( as.data.frame( c( vfinfo1, td1 ) ) ), correction = FALSE, txtfont = ffamily, pointsize = sizetxt, cex = 0.75 )
  par( opar )

######################################################
# create the text elements in the printouts
######################################################
# The two above are to delete once the graphs are generated!!!
  mainInfo   <- createviewport( "mainInfo",   left =  0.00, top =  0.00, width = 4.75, height = 0.40, pheight = mheight, pwidth = mwidth )
  infobox1   <- createviewport( "infobox1",   left =  0.25, top =  0.50, width = 2.95, height = 3.75, pheight = mheight, pwidth = mwidth )
  infobox2   <- createviewport( "infobox2",   left =  6.37, top =  0.00, width = 1.40, height = 0.40, pheight = mheight, pwidth = mwidth )
  infobox3   <- createviewport( "infobox3",   left =  3.00, top = 10.89, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textlego   <- createviewport( "textlego",   left =  5.20, top =  0.50, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textvisit0 <- createviewport( "textvisit0", left =  1.20, top =  4.60, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textvisit1 <- createviewport( "textvisit1", left =  5.20, top =  4.60, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )

# create the list and then generate the tree and "push" it
  list <- vpList( mainInfo, infobox1, infobox2, infobox3, textvisit0, textvisit1, textlego )
  tree <- vpTree( printout, list )

  pushViewport( tree )

#  seekViewport( "printout" )
#  grid.rect( gp = gpar( col = "blue" ) )

######################################################
# perimetry information
######################################################
  seekViewport( "mainInfo" )
  text <- vfinfo0$tperimetry
  if( text == "sap" ) {
    text = "Static Automated Perimetry."
  } else if( text == "fdp" ) {
    text = "Frequency-doubling Perimetry."
  } else if( text == "csp" ) {
    text = "Contrast sensitivity Perimetry."
  }
  text <- paste( text, "Legoplot progression display", sep = " " )
# ID
  text <- paste( text, "Subject ID: ", sep = "\n" )
  text <- paste( text, vfinfo0$id, ",", sep = "" )
# age
  text <- paste( text, " age: from ", round( vfinfo0$sage ), " to ", round( vfinfo1$sage ), ",", sep = "" )
# eye
  texteye <- paste( "eye:", vfinfo0$seye, sep = " " )
  if( vfinfo0$seye == "OD" ) {
    texteye <- paste( texteye, "(right)", sep = " " )
  } else if ( vfinfo0$seye == "OS" ) {
    texteye <- paste( texteye, "(left)", sep = " " )
  } else {
    texteye <- paste( texteye, "(which?)", sep = " " )
  }
  text <- paste( text, texteye, sep = " " )
  grid.text( text, x = 0.0, y = 1.0, just = c( "left", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt, fontface = "bold" ) )

######################################################
# Details about printouts
######################################################
  seekViewport( "infobox1" )
  
  text <- "info about the display, etc"
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  grid.rect( gp = gpar( col = "red" ) )

  seekViewport( "infobox2" )
  
  if( vfinfo0$tpattern == "p24d2" ) {
    textpattern <- "Central 24-2"
  } else if( vfinfo0$tpattern == "p30d2" ) {
    textpattern <- "Central 30-2"
  } else if( vfinfo0$tpattern == "p10d2" ) {
    textpattern <- "Central 10-2"
  } else if( vfinfo0$tpattern == "rnfl" ) {
    textpattern <- "RNFL-57"
  } else {
    textpattern <- "Unknown"
  }
  # algorithm
  if( vfinfo0$talgorithm == "sitas" ) {
    textalgorithm <- "SITA standard"
  } else if( vfinfo0$talgorithm == "sitaf" ) {
    textalgorithm <- "SITA fast"
  } else if( vfinfo0$talgorithm == "fullt" ) {
    textalgorithm <- "Full threshold"
  } else {
    textalgorithm <- "Unknown"
  }
  
  text <- paste( textpattern, textalgorithm, sep = "\n" )
  grid.text( text, x = 1.00, y = 1.00, just = c( "right", "top" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

######################################################
# Text for legoplots
######################################################
  seekViewport( "textlego" )
  
  text <- paste( "legoplot", sep = "" )
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

######################################################
# Details about first and last visits
######################################################
  seekViewport( "textvisit0" )
  
  text <- paste( "first ", as.character( grp ), " exams", sep = "" )
  if( grp == 1 ) text <- "first exam"
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

  seekViewport( "textvisit1" )
  
  text <- paste( "last ", as.character( grp ), " exams", sep = "" )
  if( grp == 1 ) text <- "last exam"
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

######################################################
# Details about printouts
######################################################
  seekViewport( "infobox3" )

  text <- paste( "norm vals: ", visualFields::vfenv$nv$nvname, sep = "" )
  text <- paste( text, substr( packageDescription( "visualFields" )$Date, 1, 4 ), sep = "\n" )
  text <- paste( text, "visualFields", packageDescription( "visualFields" )$Version, sep = " " )
  grid.text( text, x = 0.50, y = 0.00, just = c( "center", "bottom" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxtSmall ) )
  
# only if in save mode, then set device to off
  if( !is.null( filename ) ) {
    dev.off()
  }
}
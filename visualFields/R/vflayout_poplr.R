vflayout_poplr <- function( vf, grp = 3, nperm = 5000,
                            plotType = "vf", truncVal = 1,
                            type = "slr", typecomb = "fisher",
                            pwidth = 8.27, pheight = 11.69,
                            margin = 0.25, filename = NULL,
                            colorMapType = "pval", colorScale = NULL,
                            ringMapType  = NULL,  ringScale  = NULL,
                            imparedVision = 10, rangeNormal = NULL,
                            ffamily = "serif", sizetxt = 12,
                            sizetxtSmall = 8,
                            ffamilyvf = "serif", pointsize = 7,
                            outerSymbol = "circle", outerInch = 0.12,
                            innerSymbol = "circle", innerInch = outerInch / 1.9,
                            inch2axisunits = 12.528,
                            lengthLines = 0, thicknessLines = 0,
                            outerInchpoplr = 0.185,
                            innerInchpoplr = outerInchpoplr / 1.9,
                            lengthLinespoplr = 0, borderThickness = 1.5 ) {
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
    stop( "all visual fields should belong to the same subject and eye tested with the same perimeter and algorithm on the same locations" )
  }
  if( nrow( vf ) < 2 * grp ) stop( "the number of visual fields needs to be at least twice the number of visual fields to group for the analysis" )
  if( nrow( vf ) < 8 )       warning( "permutation analysis may not be very precise for less than 8 visual fields" )
# types of color map and ring map
  if( is.null( colorMapType) ) stop( "colorMapType must be 'slope', 'pval', or 'blind'" )
  if( colorMapType != "pval" & colorMapType != "slope" & colorMapType  != "blind" ) stop( "wrong colorMapType. Must be 'slope', 'pval', or 'blind'" )
  if( !is.null( ringMapType ) && ( ringMapType  != "pval" & ringMapType  != "slope" & ringMapType  != "blind" ) ) stop( "wrong ringMapType. Must be 'slope', 'pval', or 'blind'" )
# truncation must be between zero and one
  if( truncVal <= 0 | truncVal > 1 ) stop("truncation must be between 0 and 1")

# get the conventional color scale
  if( colorMapType == "pval" & is.null( colorScale ) ) {
    colorScale  <- visualFields::vfenv$nv$pmapsettings
  }
  if( colorMapType == "slope" & is.null( colorScale ) ) {
    colorScale         <- NULL
    colorScale$cutoffs <- c( -1.5, -1.0, -0.5, 0.5, 1 )
    colorScale$red     <- c( 0.8914642, 0.9999847, 0.9999847, 0.9742432, 0.0000000 )
    colorScale$green   <- c( 0.0000000, 0.5706177, 0.9041748, 0.9355011, 0.9999847 )   
    colorScale$blue    <- c( 0.1622925, 0.1513214, 0.0000000, 0.9213409, 0.9999847 )
    colorScale         <- as.data.frame( colorScale )
  }
  if( colorMapType == "blind" & is.null( colorScale ) ) {
    colorScale         <- NULL
    colorScale$cutoffs <- c( 5, 10, 15, 20 )
    colorScale$red     <- c( 0.8914642, 0.9999847, 0.9999847, 0.9742432 )
    colorScale$green   <- c( 0.0000000, 0.5706177, 0.9041748, 0.9355011 )   
    colorScale$blue    <- c( 0.1622925, 0.1513214, 0.0000000, 0.9213409 )
    colorScale         <- as.data.frame( colorScale )
  }
  if( !is.null( ringMapType ) && ( ringMapType == "pval" & is.null( ringScale ) ) ) {
    ringScale             <- NULL
    ringScale$cutoffs     <- c( 0.5, 1, 5 )
    ringScale$innerCircle <- c( 1, 0, 1 )
    ringScale$outerCircle <- c( 1, 1, 0 )
    ringScale             <- as.data.frame( ringScale )
  }
# get settings for the pattern of test locations
  texteval <- paste( "vfsettings$", vf$tpattern[1], sep = "" )
  settings <- eval( parse( text = texteval ) )
# get x and y locations
  texteval <- paste( vf$tperimetry[1], "locmap$",  vf$tpattern[1], sep = "" )
  locmap   <- eval( parse( text = texteval ) )
######################
# analysis
######################
# get global indices
  vfindices <- vfstats( vf )
  
# get poplr analysis
  if( plotType == "vf" ) vals  <- vf
  if( plotType == "td" ) vals  <- tdval( vf )
  if( plotType == "pd" ) vals  <- pdval( tdval( vf ) )
  pres <- poplr( vals, nperm = nperm, type = type, truncVal = truncVal, typecomb = typecomb )
# remove blind spot
  vf     <- vf[,-( settings$bs + visualFields::vfsettings$locini - 1 )]
  locmap <- locmap[-settings$bs,]
# init
  vfinfo0 <- vf[1,1:( visualFields::vfsettings$locini - 1 )]
  vfinfo1 <- vf[1,1:( visualFields::vfsettings$locini - 1 )]
  # get indices for averages
  locvalsidx <- visualFields::vfsettings$locini:( visualFields::vfsettings$locini + settings$locnum - length( settings$bs ) - 1 )
  idx0 <- c( 1:grp )
  idx1 <- c( ( nrow( vf ) - grp + 1 ):nrow( vf ) )
############################################################################
# CAK BEGIN
############################################################################
# Initialize to 1 
  nonSeenLocations <- NULL
  nonSeenLocations[1:(ncol( vf )-visualFields::vfsettings$locini+1)] <- 1
  # find locations where stimulus is not seen in all of the last n exams
  # mark these locations as 0 
  for(i in 1: length( idx1 ) )
  {  
    nonSeenLocations[which( vf[idx1[i],visualFields::vfsettings$locini:ncol( vf )] > 0 )] <- 0
  }
  # get all indices of locations which are not seen
  idxNotSeen <- which( nonSeenLocations == 1 )
#############################################################################
# CAK END
#############################################################################
# get averages  
# for blind sopts at (0,0), do not consider them for analysis
#############################################################
# CAK BEGIN
#############################################################
# get all x values of blind spot locations from the first n exams that are not zero
# compute a mean on the set of locations obtained
  idx <- which( vf$sbsx[idx0] != 0 )
  if( length( idx ) > 0 ) vfinfo0$sbsx <- mean( vf$sbsx[idx0[idx]] )
  if( length( idx ) == 0 ) vfinfo0$sbsx <- 0

# get all x values of blind spot locations from the first n exams that are not zero
# compute a mean on the set of locations obtained
  idx <- which( vf$sbsy[idx0] != 0 )
  if( length( idx ) > 0 ) vfinfo0$sbsy <- mean( vf$sbsy[idx0[idx]] )
  if( length( idx ) == 0 ) vfinfo0$sbsy <- 0

# get all x values of blind spot locations from the last n exams that are not zero
# compute a mean on the set of locations obtained  
  idx <- which( vf$sbsx[idx1] != 0 )
  if( length( idx ) > 0 ) vfinfo1$sbsx <- mean( vf$sbsx[idx1[idx]] )
  if( length( idx ) == 0 ) vfinfo1$sbsx <- 0

# get all x values of blind spot locations from the last n exams that are not zero
# compute a mean on the set of locations obtained  
  idx <- which( vf$sbsy[idx1] != 0 )
  if( length( idx ) > 0 ) vfinfo1$sbsy <- mean( vf$sbsy[idx1[idx]] )
  if( length( idx ) == 0 ) vfinfo1$sbsy <- 0
#############################################################
# CAK END
############################################################
  vfinfo0$sage <- mean( vf$sage[idx0] )
  vfinfo1$sage <- mean( vf$sage[idx1] )
  vf0          <- round( colMeans( vf[idx0, locvalsidx] ) )
  vf1          <- round( colMeans( vf[idx1, locvalsidx] ) )
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
#      windows( xpos = 0, ypos = 0, width = pwidth, height = pheight, rescale = "fixed" )
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
# sensitivity plot first n visits
  opar <- par( no.readonly = TRUE )
  par( fig = c( 0.0150, 0.5000, 0.5833, 0.9200 ) )
  color <- vfgrayscale( vf0, vfinfo0$sage, pattern = vfinfo0$tpattern, algorithm = vfinfo0$talgorithm )
  vf0[which( vf0 < 0 )] <- "<0"
  vfplotloc( vf0, eye = vfinfo0$seye, patternMap = locmap , outerColor = color, bs = c( vfinfo0$sbsx, vfinfo0$sbsy ), 
             txtfont = ffamilyvf, pointsize = pointsize,
             xminmax = xminmax, yminmax = yminmax,
             outerSymbol = outerSymbol, innerSymbol = innerSymbol,
             outerInch = outerInch, innerInch = innerInch,
             lengthLines = lengthLines, thicknessLines = thicknessLines )
# sensitivity plot last n visits
  par( new = TRUE )
  par( fig = c( 0.5000, 0.985, 0.5833, 0.9200 ) )
  color <- vfgrayscale( vf1, vfinfo1$sage, pattern =  vfinfo0$tpattern, algorithm = vfinfo0$talgorithm )
  vf1[which( vf1 < 0 )] <- "<0"
  vfplotloc( vf1, eye = vfinfo1$seye, patternMap = locmap , outerColor = color, bs = c( vfinfo1$sbsx, vfinfo1$sbsy ), 
             txtfont = ffamilyvf, pointsize = pointsize,
             xminmax = xminmax, yminmax = yminmax,
             outerSymbol = outerSymbol, innerSymbol = innerSymbol,
             outerInch = outerInch, innerInch = innerInch,
             lengthLines = lengthLines, thicknessLines = thicknessLines )
# PoPLR plot
  par( new = TRUE )
  par( fig = c( 0.1461, 0.8539, 0.1400, 0.6315 ) )
  vfplot_poplr( pres$sl, pres$pval, pres$vfdata,
                txtfont = ffamilyvf, pointsize = pointsize,
                xminmax = xminmax, yminmax = yminmax,
                outerSymbol = outerSymbol, innerSymbol = innerSymbol,
                outerInch = outerInchpoplr, innerInch = innerInchpoplr,
                lengthLines = lengthLinespoplr, thicknessLines = thicknessLines,  
                colorMapType = colorMapType, colorScale = colorScale,
                ringMapType  = ringMapType, ringScale = ringScale,
                borderThickness = borderThickness,
                idxNotSeen = idxNotSeen,
                imparedVision = imparedVision,
                rangeNormal = rangeNormal )
# plot permutation histogram
  par( new = TRUE )
  par( fig = c( 0.02, 0.40, 0.015, 0.20 ) )
  hist_poplr( pres$scomb_obs, pres$pcomb_obs, pres$scomb, txtfont = ffamily, pointsize = sizetxt )
  par( opar )
# plot md on age
  par( new = TRUE )
  par( fig = c( 0.65, 0.97, 0.015, 0.20 ) )
  # regression analysis
  progols( vfindices$tdate, vfindices$mtdev, txtfont = ffamily, pointsize = sizetxt, cex = 0.75 )
# color-code map
  if( !is.null( ringMapType ) ) {
    par( new = TRUE )
    par( fig = c( 0.82, 0.97, 0.27, 0.34 ) )
    ringmapgraph( ncol = nrow( ringScale ), mapval = ringScale, txtfont = ffamilyvf, pointsize = pointsize,
                  outerSymbol = outerSymbol, innerSymbol = innerSymbol,
                  outerInch = outerInchpoplr, innerInch = innerInchpoplr,
                  outerBorderThickness = borderThickness, innerBorderThickness = borderThickness )
  }

  if( colorMapType == "slope" ) colorScale$cutoffs <- 10 * colorScale$cutoffs
  par( new = TRUE )
  par( fig = c( 0.82, 0.97, 0.21, 0.28 ) )
  colormapgraph( ncol = 3, mapval = colorScale, notSeenAsBlack = TRUE, txtfont = ffamilyvf, pointsize = pointsize,
                 outerSymbol = outerSymbol, innerSymbol = innerSymbol,
                 outerInch = outerInchpoplr, innerInch = innerInchpoplr )
  par( opar )
######################################################
# create the text elements in the printouts
######################################################
# The two above are to delete once the graphs are generated!!!
  mainInfo      <- createviewport( "mainInfo",      left =  0.00, top =  0.00, width = 4.75, height = 0.40, pheight = mheight, pwidth = mwidth )
  infobox2      <- createviewport( "infobox2",      left =  6.37, top =  0.00, width = 1.40, height = 0.40, pheight = mheight, pwidth = mwidth )
  infobox3      <- createviewport( "infobox3",      left =  3.15, top = 10.89, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textvisit0    <- createviewport( "textvisit0",    left =  1.20, top =  0.50, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textvisit1    <- createviewport( "textvisit1",    left =  5.20, top =  0.50, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  textpoplar    <- createviewport( "textpoplar",    left =  3.20, top =  3.97, width = 1.40, height = 0.30, pheight = mheight, pwidth = mwidth )
  texthistogram <- createviewport( "texthistogram", left =  0.50, top =  9.00, width = 1.00, height = 0.30, pheight = mheight, pwidth = mwidth )
  textmdprogols <- createviewport( "textmdprogols", left =  6.00, top =  9.00, width = 1.00, height = 0.30, pheight = mheight, pwidth = mwidth )
  
# create the list and then generate the tree and "push" it
  list <- vpList( mainInfo, infobox2, infobox3, textvisit0, textvisit1, textpoplar, texthistogram, textmdprogols )
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
  text <- paste( text, "PoPLR progression analysis", sep = " " )
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
# Text for PoPLR plot
######################################################
  seekViewport( "textpoplar" )
  
  text <- paste( "PoPLR analysis", sep = "" )
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )
  
######################################################
# Text for permutation histogram
######################################################
  seekViewport( "texthistogram" )
  
  text <- paste( "permutation histogram", sep = "" )
  grid.text( text, x = 0.50, y = 0.50, just = c( "center", "center" ), gp = gpar( fontfamily = ffamily, fontsize = sizetxt ) )

######################################################
# Text for permutation histogram
######################################################
  seekViewport( "textmdprogols" )
  
  text <- paste( "mean deviation", sep = "" )
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
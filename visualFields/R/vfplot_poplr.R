vfplot_poplr <- function( sl, pval, vfinfo, newWindow = FALSE, txtfont = "mono",
                          pointsize = 7, width = 6,
                          xminmax = 29, yminmax = 29,
                          outerSymbol = "circles", innerSymbol = "circles",
                          outerSize = 1, innerSize = 1,
                          outerInch = 0.24, innerInch = 0.12,
                          lengthLines = 0, thicknessLines = 0,
                          colorMapType = "pval", colorScale = NULL,
                          ringMapType  = NULL, ringScale  = NULL,
                          imparedVision = 10, borderThickness = 2,
                          idxNotSeen = NULL, rangeNormal = NULL,
                          conormal = NULL ) {

# check that symbol input arguments are consistent
# sizes must be arrays, not matrices
  if( is.matrix(outerSize ) ) {
    stop( "Error! outerSize cannot be a matrix" )
  }
  if( is.matrix( innerSize ) ) {
    stop( "Error! innerSize cannot be a matrix" )
  }
# circles and squares are specified by one number
  if( outerSymbol == c( "circles" ) || outerSymbol == c( "squares" ) ) {
    if( length( outerSize ) != 1 ) {
      stop( "Error! length of outerSize should be 1 for circles or squares" )
    }
  }
  if( innerSymbol == c( "circles" ) || innerSymbol == c( "squares" ) ) {
    if( length( innerSize ) != 1 ) {
      stop( "Error! length of outerSize should be 1 for circles or squares" )
    }
  }
# rectangles are specified by two numbers
  if( outerSymbol == c( "rectangles" ) ) {
    if( length( outerSize ) != 2 ) {
      stop( "Error! length of outerSize should be 2 for rectangles" )
    }
  }
# rectangles are specified by two numbers
  if( innerSymbol == c( "rectangles" ) ) {
    if( length( innerSize ) != 2 ) {
      stop( "Error! length of outerSize should be 2 for rectangles" )
    }
  }
# stars are specifed by more than two numbers
  if( outerSymbol == c( "stars" ) ) {
    if( length( outerSize ) < 3 ) {
      stop( "Error! outerSize should be greater than or equal to 3 for stars" )
    }
  }
  if( innerSymbol == c( "stars" ) ) {
    if( length( innerSize ) < 3 ) {
      stop( "Error! outerSize should be greater than or equal to 3 for stars" )
    }
  }
# types of color map and ring map
  if( is.null( colorMapType) ) stop( "colorMapType must be 'slope', 'pval', or 'blind'" )
  if( colorMapType != "pval" & colorMapType != "slope" & colorMapType  != "blind" ) stop( "wrong colorMapType. Must be 'slope', 'pval', or 'blind'" )
  if( !is.null( ringMapType ) && ( ringMapType  != "pval" & ringMapType  != "slope" & ringMapType  != "blind" ) ) stop( "wrong ringMapType. Must be 'slope', 'pval', or 'blind'" )
# init
  if( is.null( conormal ) ) {
    if( colorMapType == "pval" )  conormal <- 95
    if( colorMapType == "slope" ) conormal <- 0.5
  }
# construct patternmap
  evaltxt <- paste( vfinfo$tperimetry, "locmap$", vfinfo$tpattern, sep = "" )
  patternMap <- eval( parse( text=evaltxt ) )

# remove blindspot
  evaltxt <- paste( "vfsettings$", vfinfo$tpattern, "$bs", sep = "" )
  bspos <- eval( parse( text = evaltxt ) )
  patternMap <- patternMap[-bspos,]

  innerColor          <- as.data.frame( matrix( rep( c( 1, 1, 1 ), nrow( patternMap ) ), nrow( patternMap ) ) )
  names( innerColor ) <- c( "red", "green", "blue" )

# create a new window and plot data in it
# window rescale is set to fixed to ensure re-sizing window doesn't re-size the plot
  height <- width * yminmax / xminmax
  if( newWindow ) {
    device <- options( "device" )
    if( .Platform$OS.type == "unix" ) {
      if( Sys.info()["sysname"] == "Darwin" ) {
        options( device = "quartz" )
        dev.new( width = width, height = height, dpi = 85 )
      } else {
        options( device = "x11" )
        dev.new( width = width, height = height )
      }
    } else{
      options( device = "windows" )
      dev.new( width = width, height = height, rescale = "fixed" )
    }
    options( device = device )
  }

# init all color values to white
  pval <- 100 * pval
  pvalc <- rep( c( 100 ), length( pval ) )
  pvalc[which( pval <= visualFields::vfenv$nv$pmapsettings$cutoffs[1] )] <- visualFields::vfenv$nv$pmapsettings$cutoffs[1]
  for( i in 2:( length( visualFields::vfenv$nv$pmapsettings$cutoffs ) - 1) ) pvalc[which( pval > visualFields::vfenv$nv$pmapsettings$cutoffs[i-1] & pval <= visualFields::vfenv$nv$pmapsettings$cutoffs[i] )] <- visualFields::vfenv$nv$pmapsettings$cutoffs[i]

# get the conventional color scale
  if( colorMapType == "pval" ) {
    if( is.null( colorScale ) ) {
      colorScale  <- visualFields::vfenv$nv$pmapsettings
    }
    valForMapping <- pvalc
  }
# inform the color scale for slopes
  if( colorMapType == "slope" ) {
    if( is.null( colorScale ) ) {
      colorScale         <- NULL
      colorScale$cutoffs <- c( -1.5, -1.0, -0.5, 0.5, 1 )
      colorScale$red     <- c( 0.8914642, 0.9999847, 0.9999847, 0.9742432, 0.0000000 )
      colorScale$green   <- c( 0.0000000, 0.5706177, 0.9041748, 0.9355011, 0.9999847 )   
      colorScale$blue    <- c( 0.1622925, 0.1513214, 0.0000000, 0.9213409, 0.9999847 )
      colorScale         <- as.data.frame( colorScale )
    }
# map slope values to corresponding categories defined by colorScale$cutoffs
    slc <- NULL
    slc[c( 1:length( sl ) )] <- NA
    slc[which( sl <= colorScale$cutoffs[1] )] <- colorScale$cutoffs[1]
    slc[which( sl > colorScale$cutoffs[length( colorScale$cutoffs )] )] <- colorScale$cutoffs[length( colorScale$cutoffs )]
    for( i in 2:length( colorScale$cutoffs ) ) {
      slc[which( sl > colorScale$cutoffs[i-1] & sl <= colorScale$cutoffs[i] )] <- colorScale$cutoffs[i]
    }
    valForMapping <- slc
  }
# inform the color scale for years blind
  if( colorMapType == "blind" ) {
    if( is.null( colorScale ) ) {
      colorScale         <- NULL
      colorScale$cutoffs <- c( 5, 10, 15, 20 )
      colorScale$red     <- c( 0.8914642, 0.9999847, 0.9999847, 0.9742432 )
      colorScale$green   <- c( 0.0000000, 0.5706177, 0.9041748, 0.9355011 )   
      colorScale$blue    <- c( 0.1622925, 0.1513214, 0.0000000, 0.9213409 )
      colorScale         <- as.data.frame( colorScale )
    }
# calculate years to blindness
    sens                     <- as.numeric( vfinfo[visualFields::vfsettings$locini:ncol( vfinfo ) ] )
    sens[idxNotSeen]         <- 0
    yearstoblind             <- NULL
    yearstoblind             <- ( sens - imparedVision ) / -sl
    yearstoblind[which( yearstoblind < 0 & sl >= 0 )]    <- Inf
    yearstoblind[idxNotSeen] <- 0

# get category of years blind
    yearstoblindc <- NULL
    yearstoblindc[c( 1:length( yearstoblind ) )] <- NA
    yearstoblindc[which( yearstoblind <= colorScale$cutoffs[1] )] <- colorScale$cutoffs[1]
    yearstoblindc[which( yearstoblind > colorScale$cutoffs[length( colorScale$cutoffs )] )] <- colorScale$cutoffs[length( colorScale$cutoffs )]
    for( i in 2:length( colorScale$cutoffs ) ) {
      yearstoblindc[which( yearstoblind > colorScale$cutoffs[i-1] & yearstoblind <= colorScale$cutoffs[i] )] <- colorScale$cutoffs[i]
    }
    valForMapping <- yearstoblindc
  }

  outerColor <- vfcolormap( valForMapping, mapval = colorScale )

# init ring color values
  innerBorderThickness <- borderThickness
  outerBorderThickness <- borderThickness
  innerBorderColor     <- innerColor
  outerBorderColor     <- outerColor
  
# inform the color scale for slopes
  if( !is.null( ringMapType ) && ringMapType == "pval" ) {
    if( is.null( ringScale ) ) {
      ringScale             <- NULL
      ringScale$cutoffs     <- c( 0.5, 1, 5 )
      ringScale$innerCircle <- c( 1, 0, 1 )
      ringScale$outerCircle <- c( 1, 1, 0 )
      ringScale             <- as.data.frame( ringScale )
    }
    idx <- which( pvalc <= ringScale$cutoffs[1] )
    if( length( idx ) > 0 ) {
      if( ringScale$innerCircle[1] == 1 ) innerBorderColor[idx,] <- c( 0, 0, 0 )
      if( ringScale$outerCircle[1] == 1 ) outerBorderColor[idx,] <- c( 0, 0, 0 )
    }
    for( i in 2:nrow( ringScale ) ) {
      idx <- which( pvalc > ringScale$cutoffs[i-1] & pvalc <= ringScale$cutoffs[i] )
      if( length( idx ) > 0 ) {
        if( ringScale$innerCircle[i] == 1 ) innerBorderColor[idx,] <- c( 0, 0, 0 )
        if( ringScale$outerCircle[i] == 1 ) outerBorderColor[idx,] <- c( 0, 0, 0 )
      }
    }
  }

# inform the ring scale for slopes
  if( !is.null( ringMapType ) && ( ringMapType == "slope" & is.null( ringScale ) ) ) stop( "not implemented yet" )
# inform the ring scale for years blind
  if( !is.null( ringMapType ) && ( ringMapType == "blind" & is.null( ringScale ) ) ) stop( "not implemented yet" )
  
# if some are in range normal, then restore defaults. Only if colorMapType is slope
  if( length( rangeNormal ) > 0 ) {
    idxNormal <- which( sl >= rangeNormal[1] & sl <= rangeNormal[2] )
    if( length( idxNormal ) > 0 ) {
      innerColor[idxNormal,]       <- 1
      outerColor[idxNormal,]       <- colorScale[which( colorScale$cutoff == conormal ),c(2:4)]
      innerBorderColor[idxNormal,] <- innerColor[idxNormal,]
      outerBorderColor[idxNormal,] <- outerColor[idxNormal,]
    }
  }
# non-seen are black
  if( length( idxNotSeen ) > 0 ) {
    innerBorderColor[idxNotSeen,] <- 0
    outerBorderColor[idxNotSeen,] <- 0
    outerColor[idxNotSeen,] <- 0
    innerColor[idxNotSeen,] <- 0    
  }

# slope is in dB per 10 years
  sl <- round( 10 * sl )

  vfplotloc( sl, eye = vfinfo$seye, patternMap = patternMap,
             outerColor = outerColor, innerColor = innerColor,
             txtfont = txtfont, pointsize = pointsize,
             xminmax = xminmax, yminmax = yminmax,
             outerSymbol = outerSymbol, innerSymbol = innerSymbol,
             outerSize = outerSize, innerSize = innerSize,
             outerInch = outerInch, innerInch = innerInch,
             lengthLines = lengthLines, thicknessLines = thicknessLines,
             outerBorderColor = outerBorderColor, innerBorderColor = innerBorderColor,
             outerBorderThickness = outerBorderThickness, innerBorderThickness = innerBorderThickness )
}
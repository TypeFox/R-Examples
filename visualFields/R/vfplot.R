vfplot <- function( vf, plotType, notSeenAsBlack = TRUE, newWindow = FALSE,
                    txtfont = "serif", pointsize = 7, width = 6,
                    xminmax = 29, yminmax = 29,
                    outerSymbol = "circles", innerSymbol = "circles",
                    outerSize = 1, innerSize = 1,
                    outerInch = 0.14, innerInch = 0.08,
                    lengthLines = 4.25, thicknessLines = 2 ) {

# check that vf has only 1 entry
  if( nrow( vf ) > 1 ) {
    stop("Error! vf cannot have more than 1 rows")
  }

# Check that symbol input arguments are consistent
# sizes must be arrays, not matrices
  if( is.matrix( outerSize ) ) {
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

# construct the pattern string based on the pattern type
  evaltxt <- paste( "vfsettings$", vf$tpattern, "$locnum", sep = "" )
  loc_num <- eval( parse( text = evaltxt ) )

# construct patternmap
  evaltxt <- paste( vf$tperimetry, "locmap$", vf$tpattern, sep = "" )
  patternMap <- eval( parse( text=evaltxt ) )
# get bs
  if( plotType != "vf" ) {
    evaltxt <- paste( "vfsettings$", vf$tpattern, "$bs", sep = "" )
    bspos <- eval( parse( text = evaltxt ) )
  }

# read in the plotType and decide what to do
  if( plotType == "vf" ) {
    dev  <- vf
  }

# if plot type id 'td' then calculate total deviation and total deviation probability
  if( plotType == "td" ) {
    dev  <- tdval( vf )
    devP <- tdpmap( dev )
  }  
    
# if plot type id 'pd' then first calculate total deviation
# use the toal deviation to calculate probabilty deviation 
  if( plotType == "pd" ) {
    dev  <- tdval( vf )
    dev  <- pdval( dev )
    devP <- pdpmap( dev )
  }

# if plot type id "pdghr" then first calculate total deviation
# use the total deviation to calculate pattern deviation from global-sensitivity estimate
  if( plotType == "pdghr" ) {
    dev  <- tdval( vf )
    dev  <- pdvalghr( dev )
    devP <- pdpmapghr( dev )
  }

# getRGB will return a table with the red, green and blue intensity values 
# corresponding to pattern deviation at each location
  if( plotType == "vf" ) {
    plotColor  <- vfgrayscale( dev[,visualFields::vfsettings$locini:( visualFields::vfsettings$locini + loc_num - 1 )], age = vf$sage, pattern = vf$tpattern, algorithm = vf$talgorithm )
    cloneDev   <- as.character( round( dev[,visualFields::vfsettings$locini:( visualFields::vfsettings$locini + loc_num - 1 )] ) )
    cloneDev[which( dev[,visualFields::vfsettings$locini:( visualFields::vfsettings$locini + loc_num - 1 )] < 0 )] = "<0"
  }  else {
    plotColor  <- vfcolormap( as.numeric( devP[,visualFields::vfsettings$locini:( visualFields::vfsettings$locini + loc_num - 1 )] ) )
# exclude blind spot locations
    if( all( !is.na( bspos ) ) ) dev <- dev[,-( visualFields::vfsettings$locini + bspos - 1 )]
    if( notSeenAsBlack ) {
      idxblack <- which( vf[visualFields::vfsettings$locini:( visualFields::vfsettings$locini + loc_num - 1 )] <= 0)
      if( length( idxblack ) > 0 ) plotColor[idxblack,] <- 0
    }
    lenbs <- 0
    if( all( !is.na( bspos ) ) ) lenbs <- length( bspos )
    cloneDev <- as.character( round( dev[,visualFields::vfsettings$locini:( visualFields::vfsettings$locini + loc_num - lenbs - 1 )] ) )
    if( all( !is.na( bspos ) ) ) patternMap <- patternMap[-bspos,]
    if( all( !is.na( bspos ) ) ) plotColor  <- plotColor[-bspos,]
  }

  # if NA then plot all in black
  plotColor[is.na( plotColor )] <- 0

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
  vfplotloc( cloneDev, eye = vf$seye, patternMap = patternMap, outerColor = plotColor,
             bs = c( vf$sbsx, vf$sbsy ),
             txtfont = txtfont, pointsize = pointsize,
             xminmax = xminmax, yminmax = yminmax,
             outerSymbol = outerSymbol, innerSymbol = innerSymbol,
             outerSize = outerSize, innerSize = innerSize,
             outerInch = outerInch, innerInch = innerInch,
             lengthLines = lengthLines, thicknessLines = thicknessLines )

}
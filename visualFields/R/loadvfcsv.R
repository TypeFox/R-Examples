loadvfcsv <- function( filename, patternMap ) {
# loads CSV file with visual fields and converts the columns to the correct format
# the structure must be well-defined if something changes in the structure of the
# data we need to change this program
  
# read data from file
  vf <- read.csv( filename )

# remap locations to visualFields convention according to the pattern map.
# Format is always OD (RE) and ordering from top-left to bottom-right,
# reading row-wise. The structure vfsettings used here are part of the data
# of visualField package; it has information about how the VF-object is
# expected to be
  vf2 <- vf[,visualFields::vfsettings$locini:ncol( vf )]
  vf[,visualFields::vfsettings$locini:ncol( vf )] <- vf2[,patternMap$loc]
  remove( vf2 )
  
# convert all text columns into the correct class
  vf$id         <- as.character( vf$id )
  vf$tperimetry <- as.character( vf$tperimetry )
  vf$talgorithm <- as.character( vf$talgorithm )
  vf$tpattern   <- as.character( vf$tpattern )
  vf$tdate      <- as.Date( as.character( vf$tdate ), "%m/%d/%Y")
  vf$ttime      <- as.character( vf$ttime )
  vf$stype      <- as.character( vf$stype )
  vf$sage       <- as.numeric( vf$sage )
  vf$seye       <- as.character( vf$seye )
  vf$sduration  <- as.character( vf$sduration )
  vf$spause     <- as.character( vf$spause )

# convention is that RE and LE are OD and OS
  vf$seye[ which( vf$seye == "RE" ) ] <- "OD"
  vf$seye[ which( vf$seye == "LE" ) ] <- "OS"

  return( vf )
}

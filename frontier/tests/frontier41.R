library( frontier )

## *****************************
## Testing front41WriteInput

data( front41Data )
front41Data$logOutput  <- log( front41Data$output )
front41Data$logCapital <- log( front41Data$capital )
front41Data$logLabour  <- log( front41Data$labour )

insFile <- file()
dtaFile  <- file()

front41Ins <- front41WriteInput( front41Data, "firm", yName = "logOutput",
   xNames = c( "logCapital", "logLabour" ), insFile = insFile, dtaFile = dtaFile  )

print( front41Ins )

print( readLines( insFile ) )
print( readLines( dtaFile ) )

# irregular firm (cross section) identifier
set.seed( 20061705 )
front41Data$firm <- sample( c( 1:( nrow( front41Data ) + 20 ) ) )[ 1:nrow( front41Data ) ]

front41Ins <- front41WriteInput( front41Data, "firm", yName = "logOutput",
   xNames = c( "logCapital", "logLabour" ), insFile = insFile, dtaFile = dtaFile  )

print( front41Ins )

print( readLines( insFile ) )
print( readLines( dtaFile ) )

close( insFile )
close( dtaFile )


## *****************************
## Testing front41ReadOutput

outFile <- system.file( "front41/EG1.OUT", package = "frontier" )
sfa <- front41ReadOutput( outFile )
print( coef( sfa, which = "OLS" ) )
print( coef( sfa, which = "GRID" ) )
print( coef( sfa ) )
print( summary( sfa ) )
print( summary( sfa ), efficiencies = TRUE )
print( coef( summary( sfa ), which = "OLS" ) )
print( coef( summary( sfa ), which = "GRID" ) )
print( coef( summary( sfa ) ) )
print( vcov( sfa ) )
print( sfa )
print( sfa, efficiencies = TRUE )
class( sfa ) <- NULL
print( sfa )

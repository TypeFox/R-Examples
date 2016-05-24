
library( "soiltexture" )

#   Temporary file where the info will be exported:
f <- tempfile() 

#   Generate package information
soiltextureInfo( file = f, depends = TRUE, verbose = FALSE )

#   Read again the info (as for verbose = TRUE)
cat( readLines( f ), sep = "\n" ) 
unlink( f )

#   Also works with other packages
soiltextureInfo( packages = "sp" )


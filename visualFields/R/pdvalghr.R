pdvalghr <- function( td ) {
  pd <- td

  for( i in 1:nrow( pd ) ) {
# get how many locations we need to look at
    texteval <- paste( "vfsettings$", pd$tpattern[i], "$locnum", sep = "" )
    locnum <- eval( parse( text = texteval ) )
# get PD values from obtained gh
    pd[i,visualFields::vfsettings$locini:visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] <-
      pd[i,visualFields::vfsettings$locini:( visualFields::vfsettings$locini - 1 + locnum )] + ghranktd( td[i,] )$gh
  }

  return( pd )
}

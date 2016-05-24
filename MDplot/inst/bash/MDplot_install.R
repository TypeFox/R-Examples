install.packages( "./MDplot_0.31.tar.gz", dependencies = TRUE,
                  repos = NULL, type = "source", INSTALL_opts = c( "--no-lock" ) )
if( require( MDplot ) )
{
  print( "Installation of MDplot version 0.31 succeeded!" )
}
else
{
  print( "Installation of MDplot could not be completed successfully!" )
}
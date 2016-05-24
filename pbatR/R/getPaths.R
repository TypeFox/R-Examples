####################################################################
## Thomas Hoffmann                                                 #
## CREATED:  11/??/2005                                            #
## MODIFIED: 11/15/2005                                            #
##                                                                 #
## DESCRIPTION:                                                    #
##  Allows checking in user path directories for a certain file    #
##   (i.e. the unforgiving 'pbatdata.txt').                        #
####################################################################

#################################################################
## pathGet()                                                   ##
## Returns a vector of strings for each directory in the path. ##
#################################################################
pathGet <- function() {
  paths <- c();
  if( isWindows() ) {
    paths <- strsplit( as.character(Sys.getenv("PATH")), ";", fixed=TRUE )[[1]];
    paths <- paste( paths, "\\", sep="" );
  }else{
    # Good old linux/unix/mac?
    paths <- strsplit( as.character(Sys.getenv("PATH")), ":", fixed=TRUE )[[1]];
    paths <- paste( paths, "/", sep="" );
  }

  if( length(paths) < 1 )
    warning( "'PATH' has no length." );
  return( paths );
} ## DEBUGGED

#################################################################
## pathFindFile(...)                                           ##
## PARAM  fname  filename to find                              ##
## RETURN        string of full path, or "" if not found       ##
#################################################################
pathFindFile <- function( fname ) {
  if( file.exists(fname) )
    return( fname ); ## No reason to go further - in cwd

  ## Check the rest of the path
  posFnames <- paste( pathGet(), fname, sep="" );
  ##print( posFnames );
  for( i in 1:length(posFnames) ) {
    if( file.exists( posFnames[i] ) )
      return( posFnames[i] );
  }

  return( "" );
} ## DEBUGGED

#################################################################
##                                                             ##
##                                                             ##
##                                                             ##
##                                                             ##
##                                                             ##
##                                                             ##
#################################################################

## libname alteration for user-wide installation...
getVersion <- function( libname="" ) {
  ##print( "getVersion" );
  ##print( libname );

  if( libname=="" ) libname <- NULL;

  return(  installed.packages(lib.loc=libname)["pbatR","Version"]  );

  ## if the above doesn't work, this is the previous
  ##return(  installed.packages()["pbatR","Version"]  );
}

pbat.current <- function( libname="" ){
  cat( "Checking version of pbatR... " );
  fle <- NULL
  try( {
    ##filename <- "http://www.people.fas.harvard.edu/~tjhoffm/pbatRversion.txt";
    filename <- "http://sites.google.com/site/thomashoffmannproject/pbatRversion.txt";
    fle <- file( filename );
    lines <- readLines( fle, n=3 );
    curVersion <- lines[1];
    fixes <- lines[2];
    notes <- lines[3];
    close( fle );
    fle <- NULL

    usersVersion <- getVersion(libname)
    if( curVersion == usersVersion ) {
      cat( "version is current.\n" )
    }else if( curVersion < usersVersion ) {
      cat( "You appear to be running a pre-released version (potentially unstable). Upgrade from CRAN when you no longer see this message, even if it says current.\n" )
    }else{
      cat( "version is NOT CURRENT. Consider updating (see http://www.people.fas.harvard.edu/~tjhoffm/pbatR.html for details).\n" );
      if( nchar( fixes ) > 0 ) {
        cat( "The new version fixes:\n ", fixes, "\n", sep="" )
      }else{
        cat( "No version fixes specified.\n" )
      }

      ## this function has a few issues when user-installed -- NO - PATCHED!
      ##cat( "aside: your version may be current if this was not installed as super-user if you are on a linux machine, current version is ", curVersion, ".\n" );
    }

    ##if( nchar( notes ) > 0 )
    if(!is.na(notes) && is.character(notes) && nchar(notes) > 0)
      cat( "Notes:\n ", notes, "\n", sep="" );

    return( invisible() );
  }, silent=TRUE );
  if( !is.null(fle) ) try( close(fle) );
  cat( "version check FAILED.\n" );
  return( invisible() );
}


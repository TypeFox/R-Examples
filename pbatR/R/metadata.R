####################################################################
# Thomas Hoffmann                                                  #
# CREATED:     05/28/2006                                          #
#                                                                  #
# DESCRIPTION:                                                     #
#  Preserving information such as pbat location, job splitting info#
####################################################################

CLUSTER.TIME.DEFAULT <- 0; ## 15; ## 15 seconds - no; now it just batches it and you reload into R later if it's zero

####################################################################
# pbat.get.fname(...)                                              #
# DESCRIPTION: gets $HOME/.pbat.R                                  #
# RETURN: string of what was just described.                       #
####################################################################
pbat.get.fname <- function() {
  if( Sys.getenv("R_USER")=="" )
    return( "~/.pbat.R" );  # A *nix oddity...
  return( paste( Sys.getenv("R_USER"), "/.pbat.R", sep="" ) ); #win
}

####################################################################
# SET UP THE GLOBAL VARIABLES INTERFACE                            #
####################################################################
pbatenv <- new.env();
pbatenv.set <- function( x, value )
  assign( x, value, envir=pbatenv );
pbatenv.get <- function( x, mode="any" )
  get( x, envir=pbatenv, mode=mode, inherits=FALSE );

## Metadata format
## <pbat executable name>
## <mode={single,multi,cluster}>
## <splits (for multi/cluster mode)>
## <cluster command>

## Set up the default metadata
pbat.setmode.defaults <- function( save=TRUE ){
  pbatenv.set( "executable", "pbat" );
  pbatenv.set( "mode", "single" );
  pbatenv.set( "processors", 1 );
  pbatenv.set( "cluster", "batch -f" );  ## default used to be bsub, but I still don't have access to the goddamned server...
  pbatenv.set( "refresh", CLUSTER.TIME.DEFAULT );
  pbatenv.set( "version.check", 1 ); ## yes, do it!
  pbatenv.set( "wine", "" );
  pbatenv.set( "firsttime", 1 );  ## yes, it's our first...
}

## Loading in the metadata
pbat.loadMetadata <- function(){
  filename <- paste( pbat.get.fname(), "meta", sep="" );

  if( !file.exists(filename) )
    return(); ## keep the defaults

  file <- file( filename, "r" );
  tmp <- readLines(file);
  close(file);

  if( length(tmp)<5 ){ ## previous version was 5, so set that
    warning( "Potentially malformed metadata (or from a previous version of pbatR). Please check pbat.get() to make sure everything is in order, or ensure it is still set correctly in the GUI." );
    return();
  }

  ##if( length(tmp)==6 )
  pbatenv.set( "firsttime", 0 ); ## it's been done before if we got this far

  pbatenv.set( "executable", tmp[1] );
  pbatenv.set( "mode", tmp[2] );
  pbatenv.set( "processors", tmp[3] );
  pbatenv.set( "cluster", tmp[4] );
  pbatenv.set( "refresh", tmp[5] );
  if( length(tmp)>5 )
    pbatenv.set( "version.check", tmp[6] );  ## don't lose settings from prev version
  if( length(tmp)>6 )
    pbatenv.set( "wine", tmp[7] );
}

## Writing the metadata
pbat.writeMetadata <- function(){
  filename <- paste( pbat.get.fname(), "meta", sep="" );

  file <- file( filename, "w" );
  catn( pbatenv.get("executable"), file=file );
  catn( pbatenv.get("mode"), file=file );
  catn( pbatenv.get("processors"), file=file );
  catn( pbatenv.get("cluster"), file=file );
  catn( pbatenv.get("refresh"), file=file );
  catn( pbatenv.get("version.check"), file=file );
  catn( pbatenv.get("wine"), file=file );
  close(file);
}

## Setting the mode (aside from the executable)
pbat.setmode <- function( mode=NULL, jobs=NULL, clusterCommand=NULL, clusterRefresh=NULL ){
  if( !is.null(mode) && !is.na(mode) ) {
    if( isVecContained( mode, c("single","multi","cluster") ) ){
      pbatenv.set( "mode", mode );
      if( mode=="single" )
        jobs <- 1;  ## necessary override
    }else{
      stop("'mode' must be 'single', 'multi', or 'cluster'" );
    }
  }
  if( !is.null( jobs ) && !is.na(jobs) ) {
    processors <- floor( as.numeric(jobs) );
    if( processors<1 )
      stop( "The number of jobs must be > 0" );
    pbatenv.set( "processors", processors );
  }
  print( is.null(clusterCommand) )
  if( !is.null(clusterCommand) && !is.na(clusterCommand) )
    pbatenv.set( "cluster", clusterCommand );
  if( !is.null(clusterRefresh) && !is.na(clusterRefresh) ){
    refresh <- floor(as.numeric(clusterRefresh));
    if( refresh<0 ) stop( "Cluster 'refresh' must be nonnegative." );
    pbatenv.set( "refresh", refresh );
  }

  ## overrides
  pmode <- pbat.getmode();
  if( pmode$mode=='single' & pmode$jobs!=1 )
    pbatenv.set( "processors", 1 );
  if( pmode$jobs==1 & pmode$mode!='single' )
      pbatenv.set( "mode", 'single' );

  pbat.writeMetadata();
  return( pbat.getmode() );
}

## Getting the current mode
pbat.getmode <- function(){
  list <- list( executable=pbatenv.get("executable"),
                wine=pbatenv.get("wine"),
                mode=pbatenv.get("mode"),
                jobs=as.numeric(pbatenv.get("processors")),
                cluster=pbatenv.get("cluster"),
                refresh=as.numeric(pbatenv.get("refresh")) );
  return(list);
}

## OS functionality
#isWindows <- function()
#  return( Sys.info()["sysname"]=="Windows" );
## A more robust version of the above - I'm really not sure what to expect from vista
isWindows <- function()
  return( strfindf(tolower(Sys.info()["sysname"]),"windows")!=-1
          || strfindf(tolower(R.Version()$os),"windows")!=-1
          || strfindf(tolower(R.Version()$platform),"mingw32")!=-1 )

####################################################################
#                                                                  #
# PBAT PATH GET AND SET FOR THE EXECUTABLE.                        #
#                                                                  #
####################################################################

####################################################################
# pbat.get()                                                       #
# DESCRIPTION: Returns the current stored name of the pbat         #
#   executable.  NOTE: This is saved between sessions for a given  #
#   user.                                                          #
# RETURN  current stored name of the pbat executable.              #
####################################################################
pbat.get <- function() {
  return( pbat.getmode()$executable );
}


####################################################################
# pbat.set(...)                                                    #
# DESCRIPTION: Sets the name of the pbat executable.  NOTE: This is#
#   retained between sessions for a given user.                    #
# PARAM  executableStr  String of the pbat executable name, e.g.   #
#                         "c:/pbat/pbat25.exe"                     #
####################################################################
pbat.set <- function( executableStr="", CLEAR=FALSE ) {
  if( executableStr=="" && CLEAR==FALSE ) {
    if( isPackageLoaded("tcltk") ){
      ## now do the work!
      form <- tktoplevel();
      tkwm.title( form, "P2BAT - pbat.set()" );
      if( isWindows() ) {
        executableStr <- tclvalue(tkgetOpenFile(filetypes="{{Pbat Executable} {.exe}}"));
      } else{
        executableStr <- tclvalue(tkgetOpenFile(title="PBAT binary - NO SPACES in path")); # Unix exe's have no extension!  You've been in windows too long!
      }
      if( !nchar(executableStr) ) {
        tkdestroy(form);
        warning( "Pbat not set." );
        return(invisible());
      }
      if( spaceInFilename(executableStr) ) {
        ## Special case - tcl/tk has been loaded, so do a warning message
        pbatGUI.errorMessage( paste( "There can be no spaces in the pbat executable [", executableStr, "], please rename and remove them, and try again.", sep="" ) );
        #CLEAR <- TRUE; ## done later...
      }
      tkdestroy(form);
    }else{
      ## tcltk hasn't been loaded
      cat( "Either:\n" );
      cat( "1) Supply the full filename for pbat.\n" )
      cat( "2) load the tcl/tk library with the command:\n" )
      cat( "    library(tcltk)\n" )
      cat( "   and rerun this command for a graphical file choice (recommonded unless connecting via ssh without x support.\n" );
      return();
    }
  }

  if( spaceInFilename(executableStr) ) {
    warning( paste( "There can be no spaces in the pbat executable [", executableStr, "], please rename and remove them, and try again.", sep="" ) )
    CLEAR <- TRUE  ## and clear it out
  }

  if( executableStr!="" & file.exists(executableStr)==FALSE & CLEAR==FALSE )
    warning( paste("File may not exist.  If '",
                   executableStr,
                   "' is in your path, this is safe to ignore.  Make sure you are using '\\\\' or '/'.",
                   sep="") );

  if( CLEAR==TRUE )
    executableStr <- "";

  ## alteration for new metadata format
  pbatenv.set("executable",executableStr);
  pbat.writeMetadata();
}

pbat.getwine <- function() {
  return( pbat.getmode()$wine );
}

pbat.setwine <- function( wineStr="", CLEAR=FALSE ) {
  if( wineStr=="" && CLEAR==FALSE ) {
    if( isPackageLoaded("tcltk") ){
      form <- tktoplevel()
      tkwm.title( form, "P2BAT - pbat.setwine()" )
      if( isWindows() ) warning( "You are using windows. You should not be using wine. Rerun this function with CLEAR=TRUE." )
      wineStr = tclvalue(tkgetOpenFile(title="MAC/32-bit linux ONLY, NO SPACES"))
      if( !nchar(wineStr) ) {
        tkdestroy(form)
        warning( "pbat - wine not set." );
        return(invisible())
      }
      if( spaceInFilename(wineStr) ) {
        ## Special case - tcl/tk has been loaded, so do a warning message
        pbatGUI.errorMessage( paste( "There can be no spaces in the wine executable [", wineStr, "], please rename and remove them, and try again.", sep="" ) );
        #CLEAR <- TRUE; ## done later...
      }
       tkdestroy(form)
    }else{
      ## tcltk hasn't been loaded
      cat( "Either:\n" );
      cat( "1) Supply the full filename for wine, not wineHelper.\n" )
      cat( "2) load the tcl/tk library with the command:\n" )
      cat( "    library(tcltk)\n" )
      cat( "   and rerun this command for a graphical file choice (recommonded unless connecting via ssh without x support.\n" );
      return();
    }
  }

  if( spaceInFilename(wineStr) ) {
    warning( paste( "There can be no spaces in the wine executable [", wineStr, "], please rename and remove them, and try again.", sep="" ) )
    CLEAR <- TRUE  ## and clear it out
  }

   if( wineStr!="" & file.exists(wineStr)==FALSE & CLEAR==FALSE )
    warning( paste("Cannot find a file to where wine is set (",wineStr,").",sep="") )

  if( CLEAR==TRUE ) wineStr <- "";

  pbatenv.set("wine",wineStr)
  pbat.writeMetadata()
}

## New, setting version checking information
pbat.setVersionCheck <- function( check=TRUE )
{
  pbatenv.set( "version.check", as.integer(check) );
  pbat.writeMetadata();
}

## Deprecated for the most part
pbat.getNumProcesses <- function()
  return( pbat.getmode()$jobs );


## Loading in the metadata should be done on startup!
#.onLoad <- function(libname, pkgname){
.onAttach <- function(libname, pkgname){
  pbat.setmode.defaults( save=FALSE );
  pbat.loadMetadata();
  ##print( pbat.getmode() );

  ## even newer -- firsttime user help!
  if( pbatenv.get("firsttime")==1 ) {
    pbat.firsttime();
  }

  m <- pbat.getmode();
  msg = "##############################\n"
  msg = paste(msg, "# The current pbatR mode is:\n",sep='')
  msg = paste(msg, "#  PBAT executable:", m$executable, "\n", sep='');
  if( is.null(m$wine) || m$wine=="" ) {
    msg = paste(msg, "#  (wine is not being used)\n", sep='');
  }else{
    msg = paste(msg, "#  wine location:", m$wine, "\n", sep='');
  }
  msg = paste(msg, "#  mode:", m$mode, "\n", sep='');
  msg = paste(msg, "#  jobs:", m$jobs, "\n", sep='');
  msg = paste(msg, "#  cluster command:", m$cluster, "\n", sep='');
  msg = paste(msg, "#  cluster refresh:", m$refresh, "\n", sep='');
  msg = paste(msg, "##############################\n", sep='');

  ## newest - version check
  if( pbatenv.get("version.check")==1 ) {
    pbat.current(libname);
  }else{
    msg = paste(msg, "Warning: Version checking is turned off.  Issue the command:\n", sep='');
    msg = paste(msg, " pbat.setVersionCheck()\n", sep='');
    msg = paste(msg, "to turn version checking back on.\n", sep='');
  }

  packageStartupMessage(msg)

  ##print( libname )
  ##print( pkgname )
}

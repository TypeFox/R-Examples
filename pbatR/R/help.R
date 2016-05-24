catt <- function( line, endl=FALSE ) {
  ## chops a line into pieces so it displays a little nicer for the user
  CHOP.IT.OFF <- 60;

  while( nchar(line) > CHOP.IT.OFF ) {
    ## try to break it at a space if possible
    cur.chop <- CHOP.IT.OFF;
    while( substr(line,cur.chop,cur.chop)!=' ' && cur.chop > 1 )
      cur.chop <- cur.chop-1;

    cat( substr( line, 1, cur.chop ) );
    cat( "\n" );
    line <- substr( line, cur.chop+1, nchar(line) );
  }
  cat( line );
  if( endl ) cat( "\n" );
  return( invisible() )
}
catty <- function( line )
  catt( line, endl=TRUE );

pbat.help <- function( bug=FALSE, full=TRUE, ped=NULL, phe=NULL, lib.loc=NULL ) {
  cat( "-----------------------------------------------\n" );
  cat( date(), "P2BAT Report\n" );
  cat( "-----------------------------------------------\n" );
  cat( "OS:   ", version$os, "\n", sep="" );
  cat( "ARCH: ", version$arch, "\n", sep="" );
  cat( "R:    ", version$version.string, "\n", sep="" );
  cat( "-----------------------------------------------\n" );
  packageVersion <- getVersion();
  cat( "** P2BAT v",packageVersion, "**\n", sep="" );
  cat( "pbatR v", packageVersion, "\n", sep="" );
  pbat.current();  ## does the version check
  catt( "pbat v? \n" );
  catt( "[Run the binary/executable file by clicking on it in windows to find, or running it from the command prompt if that doesn't work in linux.]\n" );
  catt( "-----------------------------------------------\n" );
  cat( "pbat binary/executable: '", pbat.get(), "'\n", sep="" );
  if( file.exists( pbat.get() ) ) {
    catt( " [[binary exists]]\n" );
  }else{
    catt( " Ensure this file is in your path -- this warning indicates that it could not be found. If you get errors similar to 'pbat not found', see pbat.set with '?pbat.set' from the R command line, and try setting it to the full path to pbat.  If you are in unix, try from the command prompt (outside of R) the command 'which pbat' or 'which pbat32', e.g., to get the full path.\n" );
  }
  mode <- pbat.getmode();
  cat( "p2bat mode: mode='", mode$mode, "', jobs='", mode$jobs, "', cluster='", mode$cluster, "', refresh='", mode$refresh, "'\n", sep="" );
  catt( "-----------------------------------------------\n" );
  catt( "Assuming the above does not solve your problem, please describe, and include the text above if possible (don't forget to fill in 'pbat version' if you don't mind), and send to the package maintainer. Thanks.\n" );
  catt( "GUI [pbat()] / command line [pbat.m(...)] / other?\n" );
  catt( "function erroring on?\n" );
  catt( "Please describe:\n" );
  catt( "\n" );
  catt( "[For further information on this package, please type '?pbat'.]\n" );
  catt( "[The webpage [http://www.people.fas.harvard.edu/~tjhoffm/pbatR.html] also provides detailed installation instructions.]\n" );
  catt( "\n" );

  ## new further debug routines to give me some more information
  if( bug ) {
    extended.pbat.debug( full );

    if( !is.null(ped) && !is.null(phe) ) {
      catt( "\nCreating obfuscated files...\n" );
      if( !is.null(ped) ){
        oped <- obfuscate(ped);
        write.ped( "obfuscate.ped", oped );
      }
      if( !is.null(phe) ){
        ophe <- obfuscate(phe);
        write.phe( "obfuscate.phe", ophe );
      }
      catt( "\nConsider including the mangled files (./obfuscate.ped, ./obfuscate.phe) if requested along with your bug report. It is suggested to ensure the bug is reproducible with these files (perhaps you wish to obfuscate them yourself first), and ensure there is nothing being shared in them you are concerned about (you can open them with any test editor).\n" );
    }
  }
}


dir.by.date <- function( ..., PRINT=FALSE ){
  files <- dir( ... );
  dates <- rep( "", length(files) );
  for( i in 1:length(files) )
    dates[i] <- file.info(files[i])$ctime;
  files <- files[ order(dates) ];

  if( PRINT )
    print( data.frame( files=files, dates=dates ) ); ## Debug only

  return( files );
}
##dir.by.date.debug <- function(){
##  dir.by.date(pattern=".*.pdf");
##}
##dir.by.date.debug();
## so dir.by.date( ".*cmd.txt" ) ought to return all the pbat command files


extended.pbat.debug <- function( FULL=FALSE ) {
  cat( "#########################\n" );
  cat( "## pbat.help()          #\n" );
  cat( "#########################\n" );
  pbat.help();
  cat( "\n\n" );

  cat( "#########################\n" );
  cat( "## pbat command file    #\n" );
  cat( "#########################\n" );
  try( {
    commandfile <- dir.by.date( pattern=".*cmd.txt" );
    if( length(commandfile) < 1 ) {
      print( "There are no pbatcmd files here." );
      return();
    }
    commandfile <- commandfile[1];
    file <- file( commandfile );
    file.lines <- readLines( file );
    print( file.lines );
    close(file);
  } );
  cat( "\n\n" );

  ## Update -- resilient to multiple processor mode...!
  cat( "#########################\n" );
  cat( "## pbat logfile         #\n" );
  cat( "#########################\n" );
  try( {
    logfile.prefix <- unlist( strsplit( file.lines[1], " " ) )[2];
    logfiles <- dir.by.date( pattern=paste(logfile.prefix,".*",sep="") );
    for( logfile in logfiles ){
      print( paste( "LOGFILE", logfile ) );
      file <- file( logfile );
      print( readLines( file ) );
      close( file );
    }
  } );
  cat( "\n\n" );

  cat( "#########################\n" );
  cat( "## pbat csvfile         #\n" );
  cat( "#########################\n" );
  try( {
    logfile.prefix <- unlist( strsplit( file.lines[1], " " ) )[2];
    file <- file( paste( logfile.prefix, ".csv", sep="" ) );
    print( readLines( file ) );
    close( file );
  } )

  cat( "#########################\n" );
  cat( "## directory listing    #\n" );
  cat( "#########################\n" );
  try( {
    dir.by.date( PRINT=TRUE );
  } );

  if( FULL ) {
    cat( "#########################\n" );
    cat( "## running pbat directly#\n" );
    cat( "#########################\n" );
    try( {
      file.remove( "pbatstatus.txt" ); ## so we get what is directly related to this
      system( paste( pbat.get(), commandfile ), intern=TRUE );   ## pbat generally doesn't give any info when it crashes here...
      file <- file( "pbatstatus.txt" );  ## sometimes there's something here...
      print( readLines( file ) );
      close( file );
    } );

    ## and try printing it just one more time in case...
    cat( "#########################\n" );
    cat( "## pbat logfile (second)#\n" );
    cat( "#########################\n" );
    try( {
      logfile <- unlist( strsplit( file.lines[1], " " ) )[2];
      print( paste( "LOGFILE ", logfile ) );
      file <- file( logfile );
      print( readLines( file ) );
      close( file );
    } );
    cat( "\n\n" );
  }
}


pbat.firsttime <- function() {
  ## Runs on the first time this library is loaded...
  ## well, that's the idea -- we might just have to point people to this?
  cat( "\n" )
  catty( "****************** first time help - begin *****************" )
  catty( "Hello new user! Let me be the first to welcome you to P2BAT!" )
  catty( "This is a one-time startup message intended to help you setup pbatR on your machine.  To see this again, just enter the command 'pbat.firsttime()' anytime." )
  cat( "\n" )
  catty( "For analysis purposes, first ensure that you have been to" )
  catty( "  http://www.biostat.harvard.edu/~clange/default.htm" )
  catty( "and have accepted the PBAT license, downloaded PBAT to your computer, and extracted PBAT somewhere. Make a note of this location, you will need to set this in this location in this software." )
  catty( "[LINUX SPECIAL NOTE: if you downloaded the linux version, you also need the pbatdata.txt which does not seem to be supplied with this - go ahead and get the windows version as well and take this file from there and put it in this location. Also, the current version posted is a 64bit version only for linux, and will not run on 32bit linux.]" )
  cat( "\n" )
  catty( "Now, either start up the GUI with the R command:" )
  catty( "  pbat()" )
  catty( "AND" )
  catty( "make sure you set the executable to the pbat file you downloaded (it is the first option in the GUI)." )
  cat( "\n" )
  catty( "OR" )
  cat( "\n" )
  catty( "If you are on you are on your local machine, run:" )
  catty( " library( tcltk )" )
  catty( " pbat.set()" )
  cat( "\n" )
  catty( "OR" )
  cat( "\n" )
  catty( "Run " )
  catty( " pbat.set( \"path.to.your.pbat\" )." )
  catty( "Note that with this method in windows, use '\\\\' to separate directories rather than a single '\\', i.e. 'C:\\\\apps\\\\pbat', or just use the forward slash, i.e. 'C:/apps/pbat'." )
  cat( "\n" )
  catty( "Please also try the R command" )
  catty( " pbat.help()" )
  catty( "if you are having any problems." )
  cat( "\n" )
  catty( "You may wish to also check the P2BAT page at" )
  catty( " http://www.people.fas.harvard.edu/~tjhoffm/pbatR.html" )
  cat( "\n" )
  catty( "Lastly, to start the power GUI, enter for continuous/binary" )
  catty( " pbat.power()" )
  catty( " pbat.power(\"dichotomous\")" )
  catty( "****************** first time help - end *******************" )
}

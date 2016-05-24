####################################################################
# Thomas Hoffmann                                                  #
# CREATED:     06/??/2005                                          #
# MODIFIED:    06/??/2005                                          #
#                                                                  #
# DESCRIPTION:                                                     #
#  Loading back in the output from pbat.                           #
####################################################################

## New c kludge routines
kludgeConvertAwful <- function( csv.infilename, csv.outfilename ) {
  warning( "Kludge pbat input level 2 reached (output is unfixable, padded)." );

  .C( "kludgeConvertAwful", as.character(csv.infilename), as.character(csv.outfilename) );
}

kludgeConvert <- function( csv.infilename, csv.outfilename ) {
  warning( "Kludge pbat input level 1 reached (tries to fix output, should be okay?)." );

  status = as.integer(0);
  status <- .C( "kludgeConvert", as.character(csv.infilename), as.character(csv.outfilename), status )[[3]];
  print( status )

  #if( status != 1 )
  #  kludgeConvertAwful( csv.infilename, csv.outfilename );
}

# a <intersect> b
vectorIntersection <- function( a, b ) {
  remList <- c();
  for( i in 1:length(a) ) {
    if( sum(a[i]==b) < 1 )
      remList <- c(remList, i);
  }
  if( length(remList) > 0 )
    a <- a[-remList];
  return(a);
} # DEBUGGED

# a-b
vectorSubtraction <- function( a, b ) {
  # coding only altered a touch from 'vectorIntersection'
  remList <- c();
  for( i in 1:length(a) ) {
    if( sum(a[i]==b) > 0 )
      remList <- c(remList, i);
  }
  if( length(remList) > 0 )
    a <- a[-remList];
  return(a);
}

# gets the current pbat log file
getPbatlogs <- function() {
  strs <- dir(pattern="pbatlog.*"); # regular expressions
  datStrs <- dir(pattern="pbatlog.*dat");
  headerStrs <- dir(pattern="pbatlog.*header");

  return( vectorSubtraction( vectorSubtraction( strs, datStrs ), headerStrs ) );
} # DEBUGGED

# Idea: call getPbatlogs before running, and then after running
#        to get the new name.
getPbatlog <- function( beforeLogs, afterLogs ) {
  log <- vectorSubtraction( afterLogs, beforeLogs );
  if( length(log)!=1 ) {
    if( length(log)<1 )
      stop( "Pbat terminated before a log-file could be written." );
    stop( "Two possible logs were found - if you are running pbat twice simulataneously in the same directory, bad things happen." );
  }

  return(log);
} # DEBUGGED

strsplitFix2 <- function( x, split ) {
  if( length(x) > 1 ) stop( "strSplitFix(...) only works on a single string." );
  if( length(x)==0 || x=="" ) return("");
  res=unlist( strsplit( x, split, fixed=TRUE ) ); # split, return as vector of strings

  ## slow
  for( i in 1:length(res) ){
    ##print( res[i] );
    if( substring(res[i],1,1)==" " )
      res[i] <- substring(res[i],2);
    if( substring(res[i],strlen(res[i]))==" " )
      res[i] <- substring(res[i],1,strlen(res[i])-1);
  }
  return( res[res!=""] ); # eliminate any empty strings!
}


## _much_, __much__ faster version! stupid!!!
loadPbatlog <- function( log ){
  ##pbatlogfile <- log;
  callfile <- paste( log, ".call", sep="" );
  resultfile <- paste( log, ".csv", sep="" );
  ##apped <- 0
  .C( "launchPbatlog", log, callfile, resultfile, as.integer(0) );

  pbatCall <- NULL;  pbatData <- NULL;
  try(  pbatCall <- readLines( callfile )  );

  ## Addition pointed out to me by Dave!
  ## - If there are no results... empty file
  if( file.info(resultfile)$size == 0 ) {
    ## No results!
    warning( "Empty output. Generally this indicates the number of informative families in the markers specified is below your current 'min.info' threshhold (or pbat crashed)." );
    return( list( call=pbatCall, data=NA ) );
  }## end addi

  read <- FALSE;
  try(  { pbatData <- read.csv( resultfile, strip.white=TRUE );
          read <- TRUE; } );
  if( !read ) {
    kludgeLog <- paste( resultfile, ".kludge.csv", sep="" );
    kludgeConvert( resultfile, kludgeLog );
    try(  { pbatData <- read.csv( kludgeLog, strip.white=TRUE );
            read <- TRUE }  );
    if( !read ) {
      kludgeConvertAwful( resultfile, kludgeLog );
      try(  { pbatData <- read.csv( kludgeLog, strip.white=TRUE );
              read <- TRUE }  );

    }
    if( !read )
      warning( "Data could not be read in, despite kludges." );
  }

  return( list( call=pbatCall, data=pbatData ) );
}

## This is somewhat slow, but there is no logical reason why a piece of
##  loadPbatlog.bad won't work (slightly more efficient, but not all that much better).
## If this comes to be a real problem, then I will rewrite it in C code.
loadPbatlog.slow.but.good <- function( log ){
  pbatCall <- NULL; pbatData <- NULL;

  ## important to make sure the logfile actually exists - otherwise
  ##  something went drastically wrong...
  if( !file.exists(log) ) {
    ##stop( paste("Cannot load the pbat logfile '",log,"'; file does not exist",sep="") );

    ## 06/14/2006 alteration - potentially yet another alteration in the output format... either that or PBAT is crashing...
    print( "Nonexistant pbat output file; potentially safe to ignore if running a smaller analysis with multiple processes (the output changes in every PBAT release). Ensure that the output is proper length." );
    return( NULL );
  }
  ## The .header & .dat file format has been _dropped_ - I really haven't seen it anymore! so we'll just assume that it has been dropped.

  lines <- readLines( log );

  NUMLINES <- length(lines);

  ## 01/26/2006
  if( NUMLINES < 1 ) {
    print( "Empty pbat output file; safe to ignore if running a smaller analysis with multiple processes. Ensure that the output is proper length." );
    return( NULL );
  }

  ##find the first '&' symbol
  ## The following code is more robust than it need be, but from previous
  ##  experience, this output format changes lots b/n versions
  ##  and is the source of most breaks in the program!
  and.symbol <- -1;
  for( i in 1:NUMLINES ){
    if( !is.null(lines[i]) && lines[i]!="" ) {
      if( strfindf(lines[i], "&") != -1 ){
        if( and.symbol == -1 ) and.symbol <- i;
        break;
      }
    }
  }

  ## if we couldn't find the '&' symbol, we're really screwed
  ##  That means that there is no data there!
  if( and.symbol == -1 ) {
    if( pbat.getNumProcesses() < 2 ) {
      print( "ERROR: No data could be found in the file. The pbat output is as follows:" );
      print( lines );
    }
    print( "No output in the logfile - just batch commands. (1) Pbat may have crashed. (2) You may be doing a relatively small analysis, so that some processes had nothing to do (so completely safe to ignore in that circumstance). Ensure output is proper length." );  ## Remove the warnings for debug purposes (warnings are extremely hard to track, as most errors in R).
    return(NULL);  # I'm wondering if we can cut it into so many pieces some don't do anything
  }

  ## assuming that we found it

  ## extract the call
  if( and.symbol > 1 )
    pbatCall <- lines[1:(and.symbol-1)];

  ## check if the first line is a header
  dataNames <- NULL;
  ##unlist(strsplit(lines[and.symbol],"&",fixed=TRUE)) );
  firstLine <- strsplitFix2( lines[and.symbol], "&" );
  if( firstLine[1] == "Group" ){
    dataNames <- firstLine;
    and.symbol <- and.symbol + 1;
    ##print( dataNames ); ## DEBUG only
  }

  if( and.symbol>NUMLINES && length(dataNames)>0 ) {
    ## sometimes we get an empty header! Ouch!!! Why?? Why??? Abolutely horrible design!
    ##  why the he|| can't all of the output have a bloody label?!?!
    pbatData <- data.frame( matrix( NA, 1, length(dataNames) ) ); ## It's all we can do.
    names(pbatData) <- dataNames;
    ##print( pbatData );
    return( list( call=pbatCall, data=pbatData ) );
  }

  ## Now, we know that the output starts at the '&' symbol;

  ## The following mysteriously doesn't work, so we've got to try harder...
  ## I mean it should, and in fact, it does most of the time...
  ## Seems to screw up when the line is too long
  ## What the heck, I don't get it - so do something else that's probably
  ##  a little slower but will get it done.
  ##logfile <- file( log, open="r", blocking=FALSE );
  ##if( and.symbol>1 ) readLines(logfile, and.symbol);
  ##pbatData <- read.table( logfile, sep="&", header=FALSE, strip.white=TRUE );
  ##close(logfile);

  for( i in and.symbol:NUMLINES ){
    nextLine <- strsplitFix2( lines[i], "&" );
    ##print( nextLine );
    ##print( length( nextLine ) );

    pbatData <- rbind( pbatData, nextLine );
  }

  row.names(pbatData) <- 1:nrow(pbatData);

  pbatData <- data.frame( pbatData );
  if( !is.null(dataNames) ){
    if( length(dataNames) == ncol(pbatData) ) {
      names(pbatData) <- dataNames;
    }else{
      print( dataNames );
      warning( "Data Names do not match the data!" );
    }
  }

  ##print( pbatData[1,] );
  ## Lastly, perhaps try to reformat the data into numbers?
  if( !is.null(pbatData) ){
    for( i in 1:ncol(pbatData) ){
      suppressWarnings(
                       curcol <- as.numeric(as.character(pbatData[,i]))
                       );
      if( !is.na(sum(curcol)) )
        pbatData[,i] <- curcol;
    }
  }

  return( list( call=pbatCall, data=pbatData ) );
}

loadPbatlog.bad <- function( log ) {
  pbatCall <- NULL; pbatData <- NULL;

  if( !file.exists(log) )
    stop( paste("Cannot load the pbat logfile '",log,"'; file does not exist",sep="") );

  # If .header & .dat file exist
  if( file.exists(paste(log,".dat",sep="")) && file.exists(paste(log,".header",sep="")) ) {
    # First load in the data and the header
    header <- read.table( paste(log,".header",sep=""),
                          sep="&", comment.char="", header=TRUE );
    pbatData <- read.table( log, sep="&", header=FALSE );
    print(length(pbatData)) # They don't match - What???
    print(length(header))
    #names(pbatData) <- names(header);
    warning( "header and data do not match!!!" );

    # Now load in the call
    logfile <- file( paste(log,".dat",sep=""), open="r", blocking=FALSE );
    pbatCall <- readLines(logfile);
    NUMLINES <- length(pbatCall);
    close(logfile);
  }else {
    # .header & .dat don't exist

    # First get the number of lines to prevent an infinite loop.
    #  Yes, this is unnecessarily slow, but not enough to warrant concern,
    #  and R is being difficult this morning.

    logfile <- file(log, open="r", blocking=FALSE);
    tmp <- readLines(logfile);
    NUMLINES <- length(tmp);
    close(logfile);


    header <- TRUE; ## 01/25/2006
    addiLine <- NULL;
    if( NUMLINES>0 ) {

      ## Now, start reading in the input

      logfile <- file(log, open="r", blocking=FALSE);
      on.exit(close(logfile));

      MARKERSTR <- "Group&";

      ;# read the lines in from the log file, checking for the header...
      line <- readLines( logfile, n=1 );
      namesVector <- NULL;
      lastLine=-1;
      for( i in 1:NUMLINES ){
        if( substring(line,1,strlen(MARKERSTR))==MARKERSTR ) {
          namesVector <- make.names( unlist(strsplit(line,"&",fixed=TRUE)) );
          ##print( namesVector );
          break;
        }else if( strfindf(line,"&")!=-1 ){
          ## all added 01/25/2006 for erroneous multiple processes output (i.e. the second one doesn't work at all!)
          ##print(line); stop(i); ## DEBUG ONLY

          addiLine <- unlist(strsplit(line,"&",fixed=TRUE));
          namesVector <- "BAD"; ## less alteration of code
          header <- FALSE;
          break;
        }else{
          pbatCall <- c(pbatCall, line);
          line <- readLines( logfile, n=1 );
          ##print( line );
        }
        lastLine=i;
      }

      ##print( "LINE" );
      ##print( line );
      ##    line <- readLines( logfile, n=1 );
      ###    print( line );
      ##    line <- readLines( logfile, n=1 );
      ##    print( line );
      ##    line <- readLines( logfile, n=1 );
      ##    print( line );
      ##stop( "what the hell" );

      if( !is.null(namesVector) && lastLine<NUMLINES ) {
        ##print( "hi" );
        pbatData <- read.table( logfile, header=FALSE, sep="&" );
        ##print( "bye" );
        ##if( length(namesVector)!=length(pbatData) ) {
        if( length(namesVector)!=length(pbatData) && header==TRUE ) {
          warning( "Names vector is of improper length! I do not know what to do!" );
          ##print( "Names:" );
          ##print( namesVector );
        }else{
          ## 01/25/2006
          ##names(pbatData) <- namesVector;
          if( header ) {
            names(pbatData) <- namesVector;
          }else{
            warning( paste("Could not load in header for '",log,"' (bug workaround for multiple processes; safe to ignore).") );
            pbatData <- rbind( addiLine, pbatData );
            ## strange peculiarity
            if( pbatData[2,1]==999 )
              pbatData[2,1] <- -999; ## why does this get lost?
          }
        }
        ##print( namesVector ); # DEBUG ONLY
        ;#names( pbatData ) <- namesVector;
      } else if( lastLine>=NUMLINES ) {
        ## PBAT error - no headers! Try again. NEW 11/15/2005
        pbatData <- read.table( log, header=FALSE, sep="&" );
      }
    } else{
      warning( "No logfile exists." );
      pbatCall="";
      pbatData="";
    }
  }

  return( list( call=pbatCall, data=pbatData ) );
}

loadCurrentPbatLog <- function( beforeLogs ) {
  # Get the current logs
  afterLogs <- getPbatlogs();
  # Do the difference and find the string of the most recently
  #  run log!
  strLog <- getPbatlog( beforeLogs, afterLogs );
  # Load and return that log.
  return( loadPbatlog( strLog ) );
}

loadPbatlogExtended <- function( log ) {
  ##pbatlogfile <- log;
  callfile <- paste( log, ".call", sep="" );
  resultfile <- paste( log, ".csv", sep="" );
  ##append <- 0
  numProcesses <- pbat.getNumProcesses();
  if( numProcesses==1 ) return( loadPbatlog( log ) );
  .C( "launchPbatlogExtended", log, callfile, resultfile, as.integer(numProcesses) );

  pbatCall <- NULL; pbatData <- NULL;
  try( pbatCall <- readLines( callfile ) );

  ## Addition pointed out to me by Dave!
  ## - If there are no results... empty file
  if( file.info(resultfile)$size == 0 ) {
    ## No results!
    warning( "Empty output. Generally this indicates the number of informative families in the markers specified is below your current 'min.info' threshhold (or pbat crashed)." );
    return( list( call=pbatCall, data=NA ) );
  }## end

  read <- FALSE;
  try(  { pbatData <- read.csv( resultfile, strip.white=TRUE );
          read <- TRUE; } );
  if( !read ) {
    kludgeLog <- paste( resultfile, ".kludge.csv", sep="" );
    kludgeConvert( resultfile, kludgeLog );
    try(  { pbatData <- read.csv( kludgeLog, strip.white=TRUE );
            read <- TRUE }  );
    if( !read ) {
      kludgeConvertAwful( resultfile, kludgeLog );
      try(  { pbatData <- read.csv( kludgeLog, strip.white=TRUE );
              read <- TRUE }  );

    }
    if( !read )
      warning( "Data could not be read in, despite kludges." );
  }

  return( list( call=pbatCall, data=pbatData ) );
}

## Added 01/09/2006 - works with multiple processes
loadPbatlogExtended.slower <- function( log ) {
  numProcesses <- pbat.getNumProcesses();
  if( numProcesses == 1 )
    return( loadPbatlog(log) );

  res <- loadPbatlog(paste(log,"_1_",numProcesses,sep=""));
  ##print( res$call );
  for( i in 2:numProcesses ){
    res2 <- loadPbatlog(paste(log,"_",i,"_",numProcesses,sep=""));
    if( !is.null(res2) ) { ## warnings provided elsewhere
      ## 01/26/2005 update - more than one of the calls isn't informative
      ##res$call <- list(res$call,res2$call); ## tack them all together
      if( is.null( res$data ) ) {
        res$data <- res2$data;
      }else if( !is.null(res2$data) ) {
        names( res2$data ) <- names( res$data ); ## fixes problems w/ coersion
        res$data <- rbind( res$data, res2$data );  ## since error in mult proc
      }
    }
  }

  ## Newest 01/26/2006 - get rid of NA's - NA's introduced because of empty column headers
  ##  in multiple processing mode.
  res$data <- res$data[!is.na(res$data[,1]),];

  rownames(res$data) <- 1:nrow(res$data);

  return(res);
}

## Cluster mode alteration
loadPbatlogConcatenate <- function( log, filename, clean=FALSE ) {
  numProcesses <- pbat.getNumProcesses();
  if( numProcesses == 1 ) {
    if( !clean ){
      file.copy( from=log, to=filename );
    }else{
      file.rename( from=log, to=filename );
    }
  }

  ## rename/remove the first one
  firstlog <- paste(log,"_1_",numProcesses,sep="");
  if( !clean ){
    file.copy( from=firstlog, to=filename );
  }else{
    file.rename( from=firstlog, to=filename );
  }

  ## now reload the rest of the buggers
  for( i in 2:numProcesses ){
    nextlog <- paste(log,"_",i,"_",numProcesses,sep="");
    if( file.exists( nextlog ) ){
      file.append( filename, nextlog );
      if( clean ) file.remove( nextlog );
    }else{
      cat( "Warning, not all output files exist; PBAT may have crashed or not finished running. See also 'is.finished()'\n" );
    }
  }

  ## print a message that it's been written
  cat( "Output has been concatenated and left in '", filename, "'.\n", sep="" );

  ## and return nothing
  return(invisible());
}

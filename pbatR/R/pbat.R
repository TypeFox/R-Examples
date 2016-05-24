####################################################################
# Thomas Hoffmann                                                  #
# CREATED:     06/21/2005                                          #
# MODIFIED:    07/08/2005                                          #
#                                                                  #
# DESCRIPTION:                                                     #
#  The major pbat-interface commands.                              #
#   ( pbat.files, pbat.obj, all of the S4 methods )                #
####################################################################

####################################################################
#                                                                  #
# Making pbat be a class....                                       #
#                                                                  #
####################################################################
pbat <- function( x, ... )
  UseMethod( "pbat" );

summary.pbat <- function( object, ... ) {
  x <- object; ## need for R CMD check

  # print out the pretty call
  if( !is.null(x$call) ) {
    print( "Call:" );
    print( x$call );
  }

  # now, print out the results?
  catn( "Results:" );
  print( x$results );
}

plot.pbat <- function( x, ... ) {
  if( x$fbat=="logrank" && !is.null(x$rcode) && x$rcode!="" ) {
    pbat.logrank.replot(load=x$rcode);
  } else {
    stop( "No plot available (only for logrank)." );
  }
}

print.pbat <- function( x, ... ) {
  catn( "Class members:" );
  catn( " $call             The formula used." );
  catn( " $pbat.call        The batch-file commands sent to pbat." );
  catn( " $results          ** The results of pbat execution in a data.frame object. **" );
  catn( " $results.logfile  Filename of the raw results of pbat (may also have extension .hdr & .dat); esp. useful if pbat outputs an unknown format that fails to be read in (in case this happens in later releases of PBAT). Note this is _in the current working directory_." );
  catn( " $rcode            Name of the file containing plots for logrank." );
  catn( " $fbat             'pc','gee', or 'logrank'" );
  catn( " $commandfile      Name of the file cantaining pbat batchfile." );
  catn( "Type '?pbat' for more details." );
}

write.pbat <- function(x, filename, resultsOnly=FALSE) {
  if( str.file.extension(filename,extension='csv')==filename ) {
    write.pbat.csv(x, filename, resultsOnly);
    return(invisible());
  }

  ## update 05/14/2006
  if( resultsOnly ) {
    write.table( x$results, filename, row.names=FALSE, quote=FALSE );
    return(invisible());
  }

  ## updates 01/20/2006
  f <- file( filename, "w" );
  #cat( paste("PBAT ",as.character(x$fbat),"\n",sep=""), file=f );
  if( !is.null(x$call) ) {
    cat( "** FORMULA **\n", file=f );
    write( x$call, file=f );
    cat( "\n\n", file=f );
  }
  if( !is.null(x$pbat.call) ) {
    cat( "** PBAT BATCH FILE: **\n", file=f );
    write( x$pbat.call, file=f );
    cat( "\n\n", file=f );
  }
  if( !is.null(x$results) ) {
    cat( "** PBAT RESULTS **:\n", file=f );
    write.table( x$results, file=f, sep=" & ", row.names=FALSE, quote=FALSE );  ## sep consistent with Christoph
    cat( "\n\n", file=f );
  }

  if( !is.null(x$rcode) && x$rcode!="") {
    cat( "** R SOURCE CODE (all that follows): **\n", file=f );
    close(f);
    ##write( x$rcode, file=f );
    ##cat( "\n\n", file=f );
    file.append( filename, x$rcode );
  }else{
    close(f);
  }
  return(invisible());
}

write.pbat.csv <- function(x, filename, resultsOnly=FALSE) {
  filename <- str.file.extension(filename,extension='csv');

  quotify <- function( strList ) {
    for( i in 1:length(strList) )
      strList[i] <- paste("\"",strList[i],"\"");
    return( strList );
  }

  ## update 05/14/2006
  if( resultsOnly ) {
    write.csv( x$results, filename, row.names=FALSE, quote=FALSE );
    return( invisible() );
  }

  ## updates 01/20/2006
  f <- file( filename, "w" );
  #cat( paste("PBAT ",as.character(x$fbat),"\n",sep=""), file=f );
  if( !is.null(x$call) ) {
    cat( "** FORMULA **\n", file=f );
    write( x$call, file=f );
    cat( "\n\n", file=f ); ## all fine
  }
  if( !is.null(x$pbat.call) ) {
    cat( "** PBAT BATCH FILE: **\n", file=f );
    write( quotify(x$pbat.call), file=f );
    cat( "\n\n", file=f );
  }
  if( !is.null(x$results) ) {
    cat( "** PBAT RESULTS **:\n", file=f );
    write.table( x$results, file=f, sep=",", row.names=FALSE, quote=FALSE );
    cat( "\n\n", file=f );
  }

  if( !is.null(x$rcode) && x$rcode!="") {## major changes here
    cat( "** R SOURCE CODE (all that follows; in quotes for csv format): **\n", file=f );
    close(f);
    ##printFile2FileQuotesAppend( filename, x$rcode );
    printFile2FileQuotesAppend( x$rcode, filename );  ## syntax backwards!
  }else{
    close(f);
  }

  return(invisible());
}


####################################################################
#                                                                  #
# Functions that work on the files.                                #
#                                                                  #
####################################################################

## a few little extras...
pbatFilesFixNamesExtraSub <- function( mName )
  paste( unlist( strsplit( mName, " |\\.|\t" ) ), sep="", collapse="" )
pbatFilesFixNamesExtra <- function( names ) {
  res <- mapply( pbatFilesFixNamesExtraSub, names )
  names( res ) <- NULL
  return( res )
}

writePbatstatus <- function( str ) {
  f <- file("pbatstatus.txt", open="a")
  cat(str,"\n",sep="",file=f)
  close(f)
}

## 05/23/06 - MASSIVE alterations in this coding!
############################################################################
# \name{pbat.files}                                                        #
# \description{                                                            #
#   Typically this function will not be run by the user, and instead,      #
#   you will use the related \code{\link{pbat.pc.files}},                  #
#   \code{\link{pbat.gee.files}}, \code{\link{pbat.logrank.files}}         #
#   functions.                                                             #
#                                                                          #
#   This is mostly included for future expansion.                          #
# }                                                                        #
# \arguments{                                                              #
#   \item{pedfile}{Name of the pedigree file, as a string.  Extension is   #
#     not needed, and the phenotype file is assumed to have similar        #
#     filename unless otherwise specified (by setting                      #
#     'phefile="othername.phe"'). }                                        #
#   \item{tempPbatCommand.txt}{Name of the command file for pbat to read   #
#     in.}                                                                 #
#   \item{...}{Pbat options, as described in the function                  #
#     \code{\link{pbat.create.commandfile}}. }                             #
############################################################################
pbat.files <- function( pedfile, phefile,
                        fbat="gee",
                        commandfile="",
                        logrank.outfile="",
                        preds="", preds.order="",
                        max.pheno=1,  ## relevant to output loading!
                        LOAD.OUTPUT=TRUE,
                        ... )
{
  cat( "pbat.files pedfile", pedfile, "\n" );

  curTimeStamp = getTimeStamp();
  if( isTimeStamped(pedfile) ) {
    curTimeStamp = extractTimeStamp( pedfile );
  }else if( isTimeStamped(phefile) ) {
    curTimeStamp = extractTimeStamp( phefile );
  }
  ##print( "curTimeStamp" );
  ##print( curTimeStamp );
  ## It time-stamps correctly, even with symbolic

  if( commandfile=="" )
    commandfile <- paste("pbat",curTimeStamp,"cmd.txt",sep="");

  # Create the command file
  logfile <- pbat.create.commandfile( pedfile=pedfile, phefile=phefile, fbat=fbat, ## 05/31/06 fix
                                      commandfile=commandfile,
                                      preds=preds, preds.order=preds.order,
                                      max.pheno=max.pheno,
                                      ... );
  ##print( logfile ); ## DEBUG ONLY

  if( logrank.outfile=="" & fbat=="logrank" ) { # so we get the same timestamp :)
    logrank.outfile <- paste( "pbat",curTimeStamp,".R", sep="" );
  }

  # Kill the 'spluscode.txt' file
  if( file.exists( "spluscode.txt" ) )
    file.remove( "spluscode.txt" );

  # call the system 'pbat' command
  #TMPOUT <- paste( "pbat", curTimeStamp, "output.txt", sep="" );
  ## 01/09/2006 rewrite for multiple processes
  ## 01/18/2006 fix to allow spaces in windows
  ## 01/24/2006 Windows version of system spin-locks!! removing completely
  ## 05/23/2006 Altering for potential new clustering method...
  ## 09/20/2007 Altering for wine with Macs, deleting all possibility of spaces in filename
  mode <- pbat.getmode()
  numProcesses <- mode$jobs;
  CLUSTER.TIME <- mode$refresh;

  ## 09/20/2007 addition
  if( spaceInFilename(mode$wine) || spaceInFilename(mode$executable) || spaceInFilename(pedfile) || spaceInFilename(phefile) ) {
    msg <- paste("There can be no spaces (i.e. ' ') in the filename for 1) the pbat executable [",mode$executable,"] 2) the pedigree filename [",pedfile,"] 3) the phenotype filename [",phefile,"] or (4) the wine executable (mac/32-bit linux only) [",mode$wine,"]", sep="" )
    writePbatstatus(msg) ## neat! we can fool the GUI!!!
    stop(msg)
  }

  wineStr <- pbat.getwine();
  if( wineStr!="" ) wineStr <- paste( wineStr, " ", sep="" );

  ## now go through the modes
  if( mode$mode == "single" ){
    ## The original
    clearCommands()
    #if( isWindows() ) {
    #  addCommand( paste( "\"", pbat.get(), "\" \"", commandfile, "\"", sep="" ) #);
    #}else{
    #  ##print( "SINGLE Command" );
    #  ##print( paste( pbat.get(), commandfile ) );
    #  addCommand( paste( pbat.get(), commandfile ) );
    #}

    addCommand( paste( wineStr, pbat.get(), " ", commandfile, sep="" ) );
    runCommands();
  }else if( mode$mode != "cluster" ){
    ## The original multiple spawning method
    clearCommands();
    for( i in 1:numProcesses ) {
      #if( isWindows() ) {
      #  addCommand( paste( "\"", pbat.get(), "\" \"", commandfile, "\"",
      #                    " ", i, " ", numProcesses, sep="" ) );
      #}else{
      #  ##print( "MULTIPLE Command" );
      #  ##print( paste( pbat.get(), commandfile, i, numProcesses ) );
      #  addCommand( paste( pbat.get(), commandfile, i, numProcesses ) );
      #}

      addCommand( paste( wineStr, pbat.get(), " ", commandfile, " ", i, " ", numProcesses,  sep="" ) );
    }
    runCommands();
  }else{
    ## The bsub method for clusters... rather inefficient, but it's for clusters that don't support the above.

    clearCommands(); ## remember otherwise it spin-locks

    filenameSH <- rep( "", numProcesses );
    filenameTouch <- rep( "", numProcesses );
    finished <- 0;

    for( i in 1:numProcesses ) {
      ## set the filenames
      filenameSH[i] <- paste( 'pbatCluster', curTimeStamp, '.', i, '.sh', sep="" );
      filenameTouch[i] <- paste( 'pbatCluster', curTimeStamp, '.', i, '.touch', sep="" );

      ## create the shell file
      file <- file( filenameSH[i], 'w' );
      ##print( "CLUSTER Command" );
      ##print( paste( pbat.get(), commandfile, i, numProcesses ) );
      #catn( pbat.get(), commandfile, i, numProcesses, file=file );
      catn( pbat.getwine(), pbat.get(), commandfile, i, numProcesses, file=file );
      catn( 'touch', filenameTouch[i], file=file );  ## add an is.finished() command for some people?
      close(file);

      ## run the shell file (well, add it to the queue)
      addCommand( paste( mode$cluster, " ", wineStr, filenameSH[i], sep="" ) );
    }

    ## now run all of the shell files
    runCommands();

    ## wait for all the files to have been 'touched'
    if( CLUSTER.TIME>0 ){
      while( finished != numProcesses ) {
        if( file.exists(filenameTouch[finished+1]) ){
          ## The next file in line finished
          finished <- finished+1;
        }else{
          ## It isn't finished, so sleep so it doesn't eat up CPU time
          Sys.sleep( CLUSTER.TIME );
        }
      }
    }else{
      print( "Commands have been batched. When you quit, save your workspace [q(save='yes')], and restore it later. Check whether it has finished with is.finished(res) to be sure." );
      LOAD.OUTPUT <- FALSE;
    }

    ## Delete all the .sh and .touch files
    for( i in 1:numProcesses ) {
      file.remove( filenameSH[i] );
      file.remove( filenameTouch[i] );
    }
  }

  ##if( !file.exists( TMPOUT ) )
  ##  stop( "Either pbat execution was terminated, or pbat executable couldn't be found.  In the latter case, you need the 'pbat' software - see set.pbat() for more details and a web-link to download." ); # this might not work... seems to so far though :)
  ##printFile( TMPOUT ); # So we can see the pbat output... anyway to tell if erred?

  # if logrank, plot the picture
  if( fbat=="logrank" )
    pbat.logrank.replot(save=logrank.outfile);

  # future expansion, maybe load in some of the results (if gee/pc),
  #  and sort them by conditional power...

  # Get the getPbatlog()
  ##print( "loading logfile" );
  ##print( logfile );
  ##res <- loadPbatlog( logfile );
  if( LOAD.OUTPUT==TRUE ) {
    res <- loadPbatlogExtended( logfile ); ## 01/09/2006
  }else{
    res <- NULL;
  }

  pbatObj <- list();
  pbatObj$call <- NULL; # set by upper function
  pbatObj$pbat.call <- res$call; ## This isn't as good anymore...
  pbatObj$results <- res$data;
  pbatObj$rcode <- logrank.outfile;
  pbatObj$results.logfile <- logfile;
  if( fbat!="logrank" ) pbatObj$rcode <- "";
  pbatObj$fbat <- fbat;
  pbatObj$commandfile <- commandfile; ## Addition so it can be cleaned
  ##print( names(pbatObj) );
  class(pbatObj) <- c("pbat", "list");
  ##print( names(pbatObj) );

  if( CLUSTER.TIME==0 && mode$mode=='cluster' ){
    ## new additions for is.finished
    pbatObj$filenameTouch <- filenameTouch;
    pbatObj$filenameSH <- filenameSH;
  }

  ## New 06/10/2006 - addendum
  ## res$call isn't very good, so load in the `logfile' instead?
  cfile <- file( commandfile )
  pbatObj$pbat.call <- readLines( cfile )
  close( cfile ) ## wow, causes huge problems if we forget to close...
  ## End of new 06/10/2006

  ## New 09/11/2006 - another PBAT hack to try; this time we try to header
  ##  the data that comes out of PBAT
  ###if( !is.null( pbatObj$results ) &&
  ###   ncol(pbatObj$results)==14 && names(pbatObj$results)[1]=="C0" ){
  ###  ## at least this get's the defaults, but it's really not the best
  ###  names(pbatObj$results) <-
  ###    c( "Group", "snps", "haplotype",  "hap.freq", "model",
  ###      "X..info.fam",  "FBAT.Wilcoxon", "power", "FBAT.LOGRANK", "power.1",
  ###      "weighted.FBAT.LOGRANK", "power.2", "optimal.FBAT.LOGRANK", "power.3" );
  ###  warning( "Miscommunication with PBAT - column headers guessed." );
  ###}

  ## further modified 01/19/2006
  if( !is.na( pbatObj$results ) &&
      !is.null( pbatObj$results ) && names(pbatObj$results)[1]=="C0" ){
    if( fbat=="logrank" && ncol(pbatObj$results)==14 ){
      names(pbatObj$results) <-
        c( "Group", "snps", "haplotype",  "hap.freq", "model",
           "X..info.fam",  "FBAT.Wilcoxon", "power", "FBAT.LOGRANK", "power.1",
           "weighted.FBAT.LOGRANK", "power.2", "optimal.FBAT.LOGRANK", "power.3" );
    }else if( fbat=="gee" || fbat=="pc" ){
      ## The typical name guessing
      pre.names <- c("group","snps","haplotype","hap.freq","model","X..info.fam","FBAT","FBAT.GxE","power.FBAT","power.FBAT.GxE", "WALD.main", "WALD.GxE");
      ##post.names <- c("heritability")  ## Alteration 7/7/7
      post.names <- paste( "heritability", 1:max.pheno, sep="" )
      if( fbat=="pc" ) post.names <- c( post.names, paste( "fbatpcWeight", 1:max.pheno, sep="" ) )

      phe <- NULL; mid.names <- c("AffectionStatus");
      if( phefile!="" ) {
        phe <- read.phe( phefile );
        mid.names <- c( names(phe)[-c(1,2)], "AffectionStatus");
      }
      ## now we go further, for predictors of higher order...
      ## - the next three lines come into play depending on Christoph's answer
      if( !is.null(preds.order) &&
          length(preds.order)>=1 &&
          !( is.character(preds.order[1]) && preds.order[1]=="" ) ){
        new.order <- order( match(preds,names(phe)) );
        preds <- preds[new.order];
        preds.order <- preds.order[new.order];
        for( i in length(preds.order):1 )
          if( as.integer(preds.order[i]) > 1 )
            for( o in as.integer(preds.order[i]):2 )
              post.names <- c( paste(preds[i],".",o,sep=""), post.names );
      }
      ## and paste together the names!
      guessed.names <- c( pre.names, mid.names, post.names );
      if( length(guessed.names) == ncol(pbatObj$results) ){
        names(pbatObj$results) <- guessed.names;
      }else if( length(guessed.names) > ncol(pbatObj$results) ){
        ## try again -- short output is bloody different
        post.names <- c("formula","heritability");
        guessed.names <- c( pre.names, post.names ); ## leave out the middles
        if( length(guessed.names) == ncol(pbatObj$results) ){
          names(pbatObj$results) <- guessed.names;
        }else{
          warning("Guessed names length did not match actual column width.");
        }
      }else{
        ## indicates the guessed names fall short of the number
        ##  of columns that are actually there

        ## Generally this means that there is something here with
        ##  the mi(...), but I don't yet understand then how this
        ##  is titled...
        ## Also tends to have a column of 1's at the end of the
        ##  analysis sometimes.. haven't figured out that either

        if( length(guessed.names) + 9 <= ncol(pbatObj$results) ) {
          mi.names <- c("GxE.DAG.main.effect",
                        "GxE.DAG.main.effect.std",
                        "GxE.DAG.main.effect.p",
                        "GxE.DAG.itx.effect",
                        "GxE.DAG.itx.effect.std",
                        "GxE.DAG.itx.effect.p",
                        "GxE.FBATI",
                        "GxE.main.heritability",
                        "GxE.itx.heritability");
          names(pbatObj$results) <- c(guessed.names,
                                      rep("NA", ncol(pbatObj$results)-length(guessed.names)-length(mi.names)),
                                      mi.names);

        }else{
          names(pbatObj$results) <- c(guessed.names,
                                      rep("NA", ncol(pbatObj$results)-length(guessed.names)) );
        }
      }

      ## New addition 05/07/2007
      if( all(names(pbatObj$results)[(ncol(pbatObj$results)-2):ncol(pbatObj$results)] == "NA" ) ) {
        print( "Newest input fix." );
        newNames <- names(pbatObj$results);
        newNames <- c( newNames[1:4], "HW","freqParent","HWParents", newNames[5:(length(newNames)-3)] );
        names(pbatObj$results) <- newNames;
      }
    }
    warning( "Miscommunication with PBAT - column headers guessed." );
  }

  ## weN 09/11/2006

  if( is.na(pbatObj$results)
     || is.null(pbatObj$results)
     || nrow(pbatObj$results)==0 ){
    cat( "There are no results.  If you see anything to the effect of 'pbat: command not found', use pbat.set() and set the location of pbat if you have X or windows, otherwise use pbat.set('<full path to pbat>').  Note that the pbat you have will probably have a version number on it.\n" );
    cat( "Additionally try setting min.info=0, although note that these may be numerically unstable.\n" );
  }else{
    ## 5/17/07
    try( names( pbatObj$results ) <- pbatFilesFixNamesExtra( names( pbatObj$results ) ) )

    ## 6/03/08
    modelCol <- which( names(pbatObj$results) == "model" )
    if( length( modelCol ) > 0 )
      pbatObj$results[[modelCol]] <- factor( pbatObj$results[[modelCol]], levels=0:3, labels=c("a","d","r","h" ) )
  }

  return( pbatObj );
}


## CLUSTER mode alteration 1
## Check to see if the cluster mode has finished
is.finished <- function( pbatObj=NULL, clean=TRUE ){
  if( is.null(pbatObj) ) pbatObj <- pbat.last();

  if( is.null(pbatObj$filenameTouch) || length(pbatObj)==0 )
    return( NA );

  ## see if all the filenames have been touched
  for( i in 1:length(pbatObj$filenameTouch) ){
    if( !file.exists(pbatObj$filenameTouch[i]) )
      return(FALSE);
  }

  ## and clean up if they've finished
  if( clean ){
    for( i in 1:length(pbatObj$filenameTouch) ){
      file.remove( pbatObj$filenameSH[i] );
      file.remove( pbatObj$filenameTouch[i] );
    }
  }
  return( TRUE );
}

## CLUSTER mode alteration 2
## Reloading in the output file
pbat.load <- function( pbatObj=NULL ){
  if( is.null(pbatObj) ) pbatObj <- pbat.last();

  pbatObj$results <- loadPbatlogExtended( pbatObj$results.logfile );
  return( pbatObj );
}

## CLUSTER mode alteration 3
pbat.concatenate <- function( pbatObj=NULL, filename="myResults.txt", clean=FALSE ){
  if( is.null(pbatObj) ) pbatObj <- pbat.last();

  loadPbatlogConcatenate( pbatObj$results.logfile, filename, clean );
}

####################################################################
#                                                                  #
# Functions on the 'phe' and 'ped' objects.                        #
#                                                                  #
####################################################################


#################################
## Kludge for affection status ##
affectionPhe <- function( ped, trait="affected", offset=0.0 ) {
  ## load the dataset into memory potentially
  if( is.sym(ped) ) {
    if( !is.cped(ped) ) {
      ped <- as.ped( ped, clearSym=TRUE )
    }else{
      ped <- as.cped( ped, clearSym=TRUE )
    }
  }
  ped$AffectionStatus[ped$AffectionStatus==0] <- NA

  ## create a phe file with the affectionStatus
  phe <- as.phe( data.frame( pid=ped$pid, id=ped$id, affected=as.integer(ped$AffectionStatus==2) ) )
  names(phe)[3] <- trait
  phe[[3]] <- phe[[3]] - offset

  return( phe )
}

## Ordering on phe and ped here is opposite to that of
##  pbat.files! Not so good! But we potentially break
##  others coding if we change it now.
####################################################################
# pbat.obj(..)           <EXTERNAL>                                #
# DESCRIPTION: takes a phe object and a ped "object" {not          #
#   explicitly defined classes, but as described in write.phe()    #
#   and write.ped() }                                              #
# PARAM phe         phe "object" as described in write.phe()       #
#       ped         ped "object" as described in write.ped()       #
#       fileprefix  prefix of the output datafile (phe & ped must  #
#                                                  match)          #
# (PARAM) ...  pbat options, as referenced to in the function      #
#                pbat.files().                                     #
####################################################################
pbat.obj <- function( phe, ped, file.prefix, phenos="", offset="gee", LOAD.OUTPUT=TRUE, ... ) {
  #cat("entered pbat.obj") ## debug hell

  #write.phe( paste( file.prefix, ".phe", sep="" ), phe );
  #write.ped( paste( file.prefix, ".ped", sep="" ), ped );
  #return( pbat.files( file.prefix, ... ) );

  ## Common error by people - check for it and more helpful message
  ##if( !(is.ped(ped) || is.pedlist(ped) ) || !is.phe(phe) ){
  ##  stop( "`phe' must be a phenotype object, and `ped' must be a pedigree object.  A common mistake is to pass the pedigree as the phenotype and vice versa." );
  ##}
  ## Update 02/25/2007 - phe can be empty
  if( !( is.ped(ped) || is.pedlist(ped) || is.cped(ped) ) )
    stop( "'ped' must be a pedigree object or cnv pedigree object. A common mistake is to switch the order of the phenotype and pedigree object when passing them to the function." );
  if( !is.phe(phe) )
    warning( "'phe' object is either of a wrong class, or unspecified (when you are just using AffectionStatus safe to ignore)." );

  ## Write out files to disk if necessary
  if( !is.sym(ped) ) {
    if( !is.cped(ped) ) {
      write.ped( paste( file.prefix, ".ped", sep="" ), ped );
      pedname <- file.prefix;
    }else{
      ## ??? WHAT ???
      write.cped( paste( file.prefix, ".cped", sep="" ), ped );
      pedname <- paste( file.prefix, ".cped", sep="" )
    }
  }else{
    pedname <- get.sym( ped );
  }
  ##cat( "pbatObj pedname", pedname, "\n" );

  ## new, nasty little kludge begin
  if( is.element("AffectionStatus",phenos) ){
    if( is.pped(ped) )
      stop( "AffectionStatus does not currently work with compressed pedigree files (pped files)." )

    ## we'll let the user specify an offset if they like
    newOffset = 0.0
    if( is.numeric( offset ) ) {
      newOffset = offset ## we do it manually
      offset="none"  ## and then tell pbat not to do it
    }

    if( is.null(phe) ) {
      ## sweet, then we can just do it
      phe <- affectionPhe( ped, offset=newOffset )
    }else{
      ## ugh... we're going to have to merge...
      if( is.element("affected",names(phe)) )
        stop( "You cannot use 'AffectionStatus' when the phenotype file has an element called 'affected'. Try renaming it." )

      ## make sure the phe is in memory
      if( is.sym(phe) )
        phe <- read.phe( get.sym(phe), sym=FALSE )

      phe2 <- affectionPhe( ped, offset=newOffset )
      phe <- as.phe( merge( phe, phe2, by.x=c("pid","id"), by.y=c("pid","id"), all.y=TRUE ) )
    }

    ## and replace it in the call
    phenos[phenos=="AffectionStatus"] <- "affected"
  }
  ## new, nasty little kludge end

  phename <- "";  ## don't always need a phefile
  if( !is.null(phe) && class(phe)=="phe" ) {
    if( !is.sym(phe) ) {
      write.phe( paste( file.prefix, ".phe", sep="" ), phe );
      phename <- file.prefix;
    }else{
      phename <- get.sym( phe );
    }
  }

  ##cat( "pbatObj (about to go pbat.files) pedname", pedname, "\n" );
  ## run the command
  res <- pbat.files( pedname, phename, phenos=phenos, offset=offset, LOAD.OUTPUT=LOAD.OUTPUT, ... );

  ## take note of what was symbollic (for the clean routine)
  ##  Nevermind - the clean routine doesn't need this!
  #if( is.sym(ped) )
  #  res$pedSym <- TRUE;
  #if( is.sym(phe) )
  #  res$pheSym <- TRUE;

  #if( CLEAN & LOAD.OUTPUT )
  #  pbat.clean( res ); ## doesn't delete the logrank stuff

  ## and return the result
  return( res );
}

####################################################################
# pbat.logrank.replot()                                            #
# DESCRIPTION: (Re)plots the survival graphs from running either   #
#   pbat.logrank(...) or pbat.logrank.files(...)                   #
# PARAM  save  filename to copy the current plotting commands to   #
# PARAM  load  filename to load previously 'save'd plots           #
#               [ or just call source(load) ]                      #
####################################################################
pbat.logrank.replot <- function( save="", load="" ) {
  if( save!="" && load!="" )
    stop( "You can only save or load at a time!" );

  if( save!="" ) {
    file.copy( "spluscode.txt", str.file.extension(save,extension=".R") );
  }else if(load!="" ) {
    source( str.file.extension(load,extension="R") );
  }else{
    source( "spluscode.txt" ); # call on load or otherwise
  }
}

####################################################################
# printFile(...)                                                   #
# DESCRIPTION: Like running cat <filename> from unix prompt,       #
#   but from within R.                                             #
# PARAM filename  Name of the file to print out.                   #
####################################################################
printFile <- function( filename ) {
  file = file(filename, "rt", blocking=FALSE );
  ##on.exit(close(file));
  ##print( readLines(file) ); ## Alteration 11/15/2005
  lines <- readLines(file);
  if( length(lines)>=1 ) {
    for( i in 1:length(lines) )
      cat( lines[i], "\n" );
  }
  close(file);
}

####################################################################
# printFile(...)                                                   #
# DESCRIPTION: Like running cat <filename> from unix prompt,       #
#   but from within R.                                             #
# PARAM filename  Name of the file to print out.                   #
####################################################################
printFile2FileQuotesAppend <- function( filename, filenameAppend ) {
  file = file(filename, "rt", blocking=FALSE );
  fileA = file(filenameAppend, "at", blocking=FALSE );
  ##on.exit(close(file));
  ##on.exit(close(fileA));
  ##print( readLines(file) ); ## Alteration 11/15/2005
  lines <- readLines(file);
  if( length(lines)>=1 ) {
    for( i in 1:length(lines) )
      cat( "\"", lines[i], "\"", "\n", sep="", file=fileA );
  }
  close(file); ## should fix
  close(fileA);
}

####################################################################
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
#                                                                  #
####################################################################

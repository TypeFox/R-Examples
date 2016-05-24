##################################################################
## Thomas Hoffmann                                               #
## CREATED:   2005                                               #
## MODIFIED: 05/23/06                                            #
##                                                               #
## DESCRIPTION:                                                  #
##  'ped' class files (pedigree files)                           #
##################################################################

## WARNING: HARD CODED CONSTANT in this file for the maximum
##  before going pure symbolic.

## See phe.R for explanation of symbollic modification

## Provide a wrapper to load in full dataset
fread.ped <- function( filename, ... )
  return( read.ped( filename, sym=FALSE, ... ) )

##################################################################
## S3 Methods for 'ped' class                                   ##
##################################################################

## There are two different ways to store it.. as ped and a pedlist...

ped <- function( x, ... )
  UseMethod( "ped" );

pedlist <- function( x, ... )
  UseMethod( "pedlist" );

## Modified 09/07/2006 for the pped
is.ped <- function( obj, pure.ped=FALSE ) {
  if( !pure.ped ){
    if( sum( class(obj)=="ped" ) == 1 )
      return(TRUE);
    return(FALSE);
  }else{
    return( is.ped( obj, pure.ped=FALSE ) & !is.pped(obj) );
  }
}

is.pedlist <- function( obj ) {
  if( sum( class(obj)=="pedlist" ) == 1 )
    return(TRUE);
  return(FALSE);
}

## S3 methods
print.ped <- function( x, ... ) {
  if( is.sym( x ) ) {
    if( length(x)==0 ){
      catn( "( pure SYMBOLIC reference -", get.sym(x), ")" )
    }else{
      cat( "Names: " );
      catn( names(x) );
      catn( "( SYMBOLIC reference -", get.sym(x), ")" )
    }
  }else{
    print.data.frame(x);
  }
}
print.pedlist <- function( x, ... ) {
  if( is.sym( x ) ) {
    if( length(x)==0 ){
      catn( "( pure SYMBOLIC reference -", get.sym(x), ")" )
    }else{
      cat( "Names: " );
      catn( names(x) );
      catn( "( SYMBOLIC reference -", get.sym(x), ")" )
    }
  }else{
    for( i in 1:length(x) ){
      cat( names(x)[i], "\n\n", sep="" );
      print( x[[i]] );
    }
  }
}
sort.ped <- function( x, decreasing=FALSE, ... ) {
  if( !is.sym(x) )
    return( x[ order(x$pid, x$id, decreasing=decreasing), ] )
  stop( "Data is symbolic, i.e. data not fully read into R. Try loading in with read.ped(...,sym=FALSE)) if you really want to do this." )
}


####################################################################
# read.ped(...)   <EXTERNAL>                                       #
# DESCRIPTION: Reads in the .ped pedigree file, as described in    #
#  the pbat literature.                                            #
# PARAM: filename  filename to open; does not need .ped extension. #
#        format=="ped":      RETURN dataframe with markers         #
#                                   subscripted as .a, .b          #
#                                   (class 'ped')                  #
#        format=="pedlist":   RETURN list with markers being lists #
#                                    with members $a, $b           #
#                                    (class 'pedlist')             #
#        sym                 Toggles symbolic reading (sets attr ) #
####################################################################
read.ped <- function( filename, format="ped", lowercase=TRUE, sym=TRUE, max=100, ... ) {
  #--------------------------------------------------------------------
  # http://www.biostat.harvard.edu/~clange/default.htm                -
  #  * Following the FBAT convention, PBAT pedigree files have the    -
  #    extension *.ped.                                               -
  #  * The first line of the PBAT pedigree file contains the names    -
  #    of the markers.                                                -
  #  * From the second line on, the pedigree file contains the        -
  #    pedigree data. Each line stands for one individual/subject.    -
  #    Each line of the pedigree file starts with the pedigree id,    -
  #    followed by the individual/subject id, the id of the father,   -
  #    the id of the mother, the individual's sex and affection       -
  #    status. After this information, for each marker, both marker   -
  #    alleles are listed. The order of the markers has to correspond -
  #    to the order of the marker names in the first line of the file.-
  #                                                                   -
  # * In the pedigree file, missing alleles, unknown affection status -
  #   and unknown sex are coded as 0?                                 -
  # * If an individual's mother or father is referred to in the       -
  #   pedigree file, but does not have its own entity in the pedigree -
  #   file, PBAT assumes that her/his marker alleles are missing.     -
  # (*) Missing values in the phenotype file have to be coded either as
  #    '.' or '-'.  Coding missing phenotypes as 0 will lead          -
  #     to incorrect results.                                         -
  #--------------------------------------------------------------------

  #warning( "DO WE NEED TO DO ANYTHING WITH MISSINGNESS???" );

  filename <- str.file.extension( filename, ".ped" );
  if( spaceInFilename(filename) & sym==TRUE )  ## sym==TRUE added 01/14/2008
    stop( spaceInFilenameError(filename) ) ## added 05/17/2007
  ped <- read.badheader( filename, na.strings="", lowercase=lowercase, onlyHeader=sym, max=max, ... ); # 0 is NA only for censor & sex
  firstNames <- c( "pid", "id", "idfath", "idmoth", "sex", "AffectionStatus" );

  if( sym ){
    ## overrides other settings
    if( length(ped$header) < max ){
      pedlist <- data.frame( matrix( 0, nrow=1, ncol=length(firstNames)+length(ped$header) ) );
      names( pedlist ) <- c( firstNames, ped$header );
      class( pedlist ) <- "pedlist";
      return( set.sym( pedlist, filename ) );
    }else{
      ## We hit the max on the load.. it needs to be pure symbolic
      pedlist <- data.frame(); ## make it empty
      class( pedlist ) <- 'pedlist';
      return( set.sym( pedlist, filename ) );
    }
  }

  if( format=="ped" ) {
    # Then we tack on extensions to each of the markers, so, for example,
    #  suppose we had marker m7, then the dataframe would have m7.a, m7.b;
    #  note, however, that these are unphased.

    tack.extension <- function( strvec, extension ) {
      for( i in 1:length(strvec) ){
        strvec[i] <- paste( strvec[i], extension, sep="" );
      }
      return( strvec );
    }
    alt.vecs <- function( v1, v2 ) {
      if( length(v1)!=length(v2) ) stop("alt.vecs: length must be the same!");
      v <- rep(0,2*length(v1));
      v[seq(1,2*length(v1),by=2)] <- v1;
      v[seq(2,2*length(v1),by=2)] <- v2;
      return(v);
    }

    markerNames <- alt.vecs( tack.extension(ped$header,".a"),
                             tack.extension(ped$header,".b")    );
    names( ped$table ) <- make.names(  c( firstNames, markerNames )  );
    class(ped$table) <- c("ped","data.frame");
    return( ped$table );
  }
  #else{

  # keep each of the markers the same, but each marker will be a
  #  list of 'a' and 'b' for the two.
  pedlist <- list();
  for( i in 1:length(firstNames) )
    pedlist[[i]] <- ped$table[[i]];
  names( pedlist ) <- firstNames;

  ENTRYPOINT <- length(firstNames)+1;
  for( i in 1:length(ped$header) ) {
    idx1=ENTRYPOINT+(i-1)*2;
    idx2=idx1+1;
    newMarker <- list( a=ped$table[[idx1]], b=ped$table[[idx2]] );
    pedlist[[length(pedlist)+1]] <- newMarker;
    names(pedlist)[length(pedlist)] <- ped$header[i];
  }
  class(pedlist) <- c("pedlist","list"); #OUCH!
  return(pedlist);
} ## VERIFIED ## (with 'total' dataset)

as.ped <- function( x,
                    pid="pid", id="id", idfath="idfath",
                    idmoth="idmoth", sex="sex", affection="AffectionStatus",
                    clearSym=FALSE )
{
  if( is.sym(x) ){
    if( clearSym==TRUE )
      return( read.ped( get.sym(x), sym=FALSE ) );
    return( x );
  }

  if( is.ped(x) )
    return(x);

  if( is.pedlist(x) ) {
    # convert it
    tmpdf = x;
    remList <- c();
    for( i in 1:length(x) ) {
      if( is.list(x[[i]]) ){
        remList <- c(remList, i); # to be removed later...
        tmpdf[[length(tmpdf)+1]] <- x[[i]][[1]];
        tmpdf[[length(tmpdf)+1]] <- x[[i]][[2]];
        name1 <- paste( names(x)[i], ".a", sep="" );
        name2 <- paste( names(x)[i], ".b", sep="" );
        jj <- length(tmpdf);
        names(tmpdf)[(jj-1):jj] <- c(name1,name2);
      }
    }
    df <- data.frame( tmpdf[-remList] );
    class(df) <- c("ped", "data.frame");
    return( df )
  }

  if( is.data.frame(x) ) {
    # The we just need to ensure the proper ordering...

    # ensure proper ordering
    idpedCol <- x[pid];
    idsubCol <- x[id];
    idfathCol <- x[idfath];
    idmothCol <- x[idmoth];
    idsexCol <- x[sex];
    idAffectionCol <- x[affection];
    df <- dfr.r( x, c(pid,id,idfath,idmoth,sex,affection) );
    df <- cbind( idpedCol, idsubCol, idfathCol, idmothCol,
                 idsexCol, idAffectionCol, df );
    names(df)[1:6] <- c("pid","id","idfath","idmoth","sex","AffectionStatus");
    class(df) <- c("ped", "data.frame" );
    return( df );
  }
  stop( "parameter 'x' must be of class 'ped', 'pedlist', or 'data.frame'." );
}

rem.dot.a <- function( strVec ){
  for( i in 1:length(strVec) ){
    dotloc <- strlen(strVec[i])-1;
    if( substring( strVec[i], dotloc, dotloc )!="." ) #@$%!!!
      stop( "write.ped: malformed header" );
    strVec[i] <- substring( strVec[i], 1, strlen(strVec[i])-2 ); #@$%!!!

    #if( substring( strVec[i], dotloc, dotloc )=="." ) #@$%!!!
    #  strVec[i] <- substring( strVec[i], 1, strlen(strVec[i])-2 ); #@$%!!!
  }
  return( strVec );
}

####################################################################
# write.ped(...)  <EXTERNAL>                                       #
# Writes a pedigree file out from the proper data format.          #
# PARAM file  string or connection for file output                 #
####################################################################
write.ped <- function( file, ped ) {
  # assuming its in a dataframe format...
  # also assuming that "marker.a" and "marker.b" are next to each
  #  other in the data frame...

  if( is.character(file) ){
    file <- str.file.extension(file,".ped");
  }

  if( is.sym(cped) )
    stop( "ped object is symbolic -- it was not really read into R. Thus you did not modify it, and so there is no point to doing this." )

  if( is.pedlist(ped) )
    ped <- as.ped(ped); # convert to 'ped' instead of 'pedlist'

  if( !is.ped(ped) )
    stop( "Can only write objects of class 'ped'.  See as.ped(...)" );

  #markerNames <- unique( rem.dot.a( names()[7:length(
  header <- unique( rem.dot.a( names(ped)[7:length(ped)] ) );

  # and dump it to file!
  write.badheader( file, ped, header );
}

as.pedlist <- function( x,
                        pid="pid", id="id", idfath="idfath",
                        idmoth="idmoth", sex="sex", affection="AffectionStatus",
                        clearSym=FALSE ) {
  if( is.sym(x) ){
    if( clearSym==TRUE )
      return( read.ped( get.sym(x), sym=FALSE, format="pedlist" ) );
    return( x );
  }

  if( is.pedlist(x) )
    return(x);

  if( is.data.frame(x) )
    x <- as.ped( x, pid, id, idfath, idmoth, sex, affection );

  if( !is.ped(x) )
    stop( "'x' must be of class 'ped' or 'data.frame'." );

  header <- unique( rem.dot.a( names(x)[-c(1:6)] ) );
  firstNames <- c( "pid", "id", "idfath", "idmoth", "sex", "AffectionStatus" );

  # keep each of the markers the same, but each marker will be a
  #  list of 'a' and 'b' for the two.
  pedlist <- list();
  for( i in 1:length(firstNames) )
    pedlist[[i]] <- x[[i]];
  names( pedlist ) <- firstNames;

  ENTRYPOINT <- length(firstNames)+1;
  for( i in 1:length(header) ) {
    idx1=ENTRYPOINT+(i-1)*2;
    idx2=idx1+1;
    newMarker <- list( a=x[[idx1]], b=x[[idx2]] );
    pedlist[[length(pedlist)+1]] <- newMarker;
    names(pedlist)[length(pedlist)] <- header[i];
  }
  class(pedlist) <- c("pedlist","list");
  return(pedlist);

}

## pped stuff
is.pped <- function( obj ){
  if( is.sym(obj) & file.extension(get.sym(obj))=="pped" )
    return(TRUE);
  return(FALSE);
}
read.pped <- function( filename, max=100 ){
  ## mirrors the symbolic piece of read.ped
  filename <- str.file.extension( filename, "pped" );
  if( !file.exists(filename) )
    stop( paste( "Cannot open '", filename, "' - file does not exist.", sep="" ) );

  ## addition to load in the names of the compressed format if possible
  file <- file( filename, open="r" );
  numNames <- as.numeric( readLines( file, n=1 ) );
  if( numNames < max ) {
    firstNames <- c( "pid", "id", "idfath", "idmoth", "sex", "AffectionStatus" );
    new.names <- readLines( file, n=numNames );
    close( file );
    pedlist <- data.frame( matrix( 0, nrow=1, ncol=length(firstNames)+length(new.names) ) );
    names( pedlist ) <- c( firstNames, new.names );
    class( pedlist ) <- "pedlist";
    return( set.sym( pedlist, filename ) );
  }
  close( file );

  ## otherwise back to being pure symbolic...
  pedlist <- data.frame();
  class( pedlist ) <- 'pedlist';
  return( set.sym( pedlist, filename ) );
}
as.pped <- function( ped, ppedname="" ){
  ## Get the filename
  kill <- FALSE;

  ## we _do_ need the 'pbatdata.txt' file, despite it doing
  ##  absolutely nothing!  ah well, probably need it later anyway
  checkAndGetPbatdata();

  pedname <- "killme.ped";
  if( is.sym(ped) ) {
    pedname <- get.sym( ped );
    if( ppedname=="" ){
      ppedname <- file.strip.extension( pedname ); ## get rid of the .ped
      ppedname <- paste( ppedname, ".pped", sep="" );
    }
  }else{
    kill <- TRUE;
    ## pedname is set
    if( ppedname=="" )
      ppedname <- "pped.pped";
  }

  ## make sure the file doesn't exist
  if( file.exists(ppedname) )
    stop( paste( "The file '", ppedname, "' already exists.", sep="" ) );

  ## do we need to write a temporary pedigree file?
  if( kill )
    write.ped( pedname, ped );

  ## create the batchfile to convert
  pbatfile <- file( "killme.txt", open="w" );
  cat( "pedfile ", pedname, "\n", sep="", file=pbatfile );
  cat( "xwriteped ", ppedname, "\n", sep="", file=pbatfile );
  close( pbatfile );
  #print( paste( pbat.get(), "killme.txt" ) ); ## debug only

  ##if( isWindows() ) {
  ##  system( paste( "\"", pbat.get(), "\" killme.txt", sep="") ); ## assuming don't need pbatdata.txt?
  ##}else{
  ##  system( paste( pbat.get(), "killme.txt" ) );
  ##}
  ## 09/25/2007 update
  wineStr <- pbat.getwine();
  if( wineStr != "" ) wineStr <- paste( wineStr, " ", sep="" );
  system( paste( wineStr, pbat.get(), " ", "killme.txt", sep="" ) );

  ## and kill the temp files
  file.remove( "killme.txt" );
  ## important to be very safe with this
  if( pedname == "killme.ped" & kill==TRUE )
    file.remove( "killme.ped" );

  print( pedname );
  print( ppedname );

  return( read.pped( ppedname ) );
}

## DEPRECATED: see plot.pedigree.R for an awesome rewrite
## NEW! Plotting routines
# plotPed <- function( ped, sink=NULL ) {
#   require( kinship ) ## replaced from 'library' for codetools...
# 
#   ## is it symbolic? it can't be for these routines...
#   if( is.sym(ped) )
#     ped <- as.ped( ped, clearSym=TRUE )
# 
#   ## move it to their format
#   if( any( ped$sex==0 ) )
#     ped$sex[ped$sex==0] <- 3;
# 
#   ## Huh? the documentation on this package doesn't make much sense...
#   #ped$affection <- 0
#   #ped$affection[ped$AffectionStatus==2] <- 1
#   ped$affection <- ped$AffectionStatus
# 
#   ## If sink = filename, sink each plot to a file!
#   ## See if we should sink it to file
#   if( !is.null(sink) ) {
#     pdf( sink );
#   }else{
#     par(ask=TRUE);
#   }
# 
#   for( pid in unique(ped$pid) ) {
#     ## pull out the pedigree piece
#     subPed <- ped[ ped$pid==pid, ]
#     ## fix it so it's their program happy
#     pedigr <- pedigree(id = subPed$id, dadid=subPed$idfath, momid=subPed$idmoth, sex=subPed$sex, affected=subPed$affection)
#     print( pedigr )
#     print( str(pedigr) )
# 
#     SUCCESS <- FALSE;  ## sometimes it fails...
#     try( {
#       plot( pedigr );
#       title( pid );
#       SUCCESS <- TRUE;
#     } );
#     if( !SUCCESS )
#       print( paste( "Plotting pedigree", pid, "failed." ) );
#   }
# 
#   ## Close off the filename if necessary
#   if( !is.null(sink) )
#     dev.off()
# }


## 12/07/2008
ped.markerNames <- function( ped ) {
  if( is.sym(ped) ) {
    if( length(ped)==0 )
      stop( "Pedigree is completely symbolic, try reading it in with 'fread.ped' instead of 'read.ped'.")
    return( names(ped)[7:length(ped)] )
  }
  
  n <- names(ped)[seq(from=7, to=ncol(ped), by=2)]
  return(  substring( n, 1, nchar(n)-2 )  )
}

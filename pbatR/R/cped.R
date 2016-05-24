## Thomas Hoffmann
## Created:  04/19/2008
## Modified: same
##
## Reads in cped files (copy number variation).
## Again -- still allow the symbolic loading

## Tom, try to code missing CNVs by -1234.0. That should work. Best, Christoph

## We have one cped file where there are the usual 6 columns
## as in a pedigree file, followed by 12 columns, which is
## data for two cped markers. It is read in with
##  cped 2 3
## This indicates that each marker has three intensities,
## the second is selected for analysis.
## So, we must have three intensities, and each intensity
## has two values for each marker. So what might make sense,
## for the naming convention to be consistent with before,
## Suppose we have our marker is 'c1' and 'c2'
##   c1.1 c1.2 c1.3 c2.1 c2.2 c2.3.a
## class 'cped'

## Did Christoph ever fix that damned AffectionStatus bug? That's going to haunt us here...

## What does the
##  channing 1
## option do?

## Provide a wrapper to load in full dataset
fread.cped <- function( filename, ... )
  return( read.cped( filename, sym=FALSE, ... ) )

## S3 Integration
cped <- function( x, ... )
  UseMethod( "cped" )

## taken from is.phe
is.cped <- function( obj ) {
  if( sum( class(obj) == "cped" ) == 1 )
    return( TRUE )
  return( FALSE )
}## DEBUGGED

## essentially taken from print.ped
print.cped <- function( x, ... ) {
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
}## DEBUGGED

sort.cped <- function( x, decreasing=FALSE, ... ) {
  if( !is.sym(x) )
    return( x[ order(x$pid, x$id, decreasing=decreasing), ] )
  stop( "Data is symbolic, i.e. data not fully read into R. Try loading in with read.cped(...,sym=FALSE)) if you really want to do this." )
}## DEBUGGED

## EXTERNAL
read.cped <- function( filename, lowercase=TRUE, sym=TRUE, max=100, ... ) {
  filename <- str.file.extension( filename, ".cped" )
  if( spaceInFilename(filename) && sym==TRUE )
    stop( spaceInFilenameError(filename) )
  cped <- read.badheader( filename, na.strings=c("-",".","NA"), lowercase=lowercase, onlyHeader=sym, max=max, ... );
  firstNames <- c( "pid", "id", "idfath", "idmoth", "sex", "AffectionStatus" )

  if( sym ) {
    cpedlist <- NULL
    if( length(cped$header) < max ) { ## Stupid bug! cped, not ped!!! How on earth did this _ever_ work???
      cpedlist <- data.frame( matrix( 0, nrow=1, ncol=length(firstNames)+length(cped$header) ) )
      names( cpedlist ) <- c( firstNames, cped$header )
    }else{
      ## We hit the max on the load... needs to be _pure_ symbolic
      cpedlist <- data.frame()
    }
    class( cpedlist ) <- c("cped","data.frame")
    return( set.sym( cpedlist, filename ) );
  }

  NC <- ncol(cped$table) - 6
  ##cat( "read.cped NC" , NC, "\n" )
  cnvNames <- rep("",NC)
  #numCnvs <- length(cped$header)
  numIntensity <- NC / length(cped$header)
  ##suffixes <- c(".a",".b")## another mistake by C! ARGH!!!!

  for( i in 1:numIntensity )
    cnvNames[ seq( from=i, to=NC, by=numIntensity ) ] <- paste( cped$header, ".", i, sep="" )

  names( cped$table ) <- make.names( c( firstNames, cnvNames ) )
  class( cped$table ) <- c("cped","data.frame")

  attr( cped$table, "numIntensity" ) <- numIntensity  ## Doesn't really do that much

  #print( cped$table ) ## DEBUG ONLY

  ## NEW! Need to recode missing data...
  for( i in 7:ncol(cped$table) ) {
    wh <- which(cped$table[[i]]==-1234.0)
    if( length(wh) > 0 )
      cped$table[[i]][wh] <- NA
  }
  ## End of recode...

  return( cped$table )
}## DEBUGGED

## as.cped
as.cped <- function( x,
                     pid="pid", id="id", idfath="idfath",
                     idmoth="idmoth", sex="sex", affection="AffectionStatus",
                     clearSym=FALSE ) {
  ## usage for clearing the sym-ness of it
  if( is.sym(x) ) {
    if( clearSym==TRUE )
      return( read.cped( get.sym(x), sym=FALSE ) )
    return( x )
  }

  ## usage if already desired object
  if( is.cped(x) )
    return( x )

  if( is.data.frame(x) ) {
    ## Then we just need to ensure the proper ordering...

    idpedCol <- x[pid]
    idsubCol <- x[id]
    idfathCol <- x[idfath]
    idmothCol <- x[idmoth]
    idsexCol <- x[sex]
    idAffectionCol <- x[affection]
    df <- dfr.r( x, c(pid,id,idfath,idmoth,sex,affection) ) ## removes them
    df <- cbind( idpedCol, idsubCol, idfathCol, idmothCol, idsexCol, idAffectionCol, df )
    names(df)[1:6] <- c("pid","id","idfath","idmoth","sex","AffectionStatus")
    class(df) <- c("cped","data.frame")
    return( df )
  }

  stop( "parameter 'x' must be of class 'cped', or 'data.frame'." );
}## DEBUGGED (by the example!)

## Sigh... this is backwards for all other
##  R routines... it kind of drives me crazy
##  but we probably shouldn't change
##  the convention now?
write.cped <- function( file, cped ) {
  ## if character, enforce the extension
  if( is.character(file) )
    file <- str.file.extension( file, ".cped" )

  if( is.sym(cped) )
    stop( "cped object is symbolic -- it was not really read into R. Thus you did not modify it, and so there is no point to doing this." )

  ## tack off the .a, and then the .1, .2, .3 pieces
  ##header <- unique( rem.dot.a( rem.dot.a( names(cped)[7:ncol(cped)] ) ) )
  header <- unique( rem.dot.a( names(cped)[7:ncol(cped)] ) )

  ## and dump it to file
  write.badheader( file, cped, header, na="-1234.0" )  ## na update for cped missing data!
}## DEBUGGED

## plotting it?
plotCPed <- function( cped, sink=NULL ) {
  if( is.sym(cped) )
    cped <- as.cped( cped, clearSym=TRUE )

  ## Want to piggyback on other routine..
  ##  this is a complete hack!
  class( cped ) <- c("ped", "data.frame")
  plotPed( cped, sink=sink )
}## DEBUGGED

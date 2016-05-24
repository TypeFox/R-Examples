####################################################################
# Thomas Hoffmann                                                  #
# CREATED:  2005                                                   #
# MODIFIED: 05/23/06                                               #
#                                                                  #
# DESCRIPTION:                                                     #
#  'phe' object files.                                             #
####################################################################

## The idea of symbolicially loading files was thought of after
##  all of this infrastructure was laid.  We need the names(...)
##  function to still do what it used to, or we need to modify a
##  lot of things. Thus this is slightly complicated in that
##  we set the 'sym' attribute of an object to indicate its
##  just 'pointing' to that file if you will. Otherwise the file is
##  redundantly written back to disk.

## Provide a wrapper to load in full dataset
fread.phe <- function( filename, ... )
  return( read.phe( filename, sym=FALSE, ... ) )

## On to the usual code:

phe <- function( x, ... )
  UseMethod( "phe" );

is.phe <- function( obj ) {
  if( sum( class(obj)=="phe" ) == 1 )
    return(TRUE);
  return(FALSE);
}

print.phe <- function( x, ... ) {
  if( is.sym( x ) ) {
    cat( "Names: " );
    catn( names(x) );
    catn( "( SYMBOLIC reference -", get.sym(x), ")" )
  }else{
    print.data.frame(x);
  }
}
sort.phe <- function( x, decreasing=FALSE, ... ) {
  if( !is.sym(x) )
    return( x[ order(x$pid, x$id, decreasing=decreasing), ] )
  stop( "Not symbolic, i.e. data not read into R. Try loading in with read.phe(...,sym=FALSE)) if you really want to do this." )
}


####################################################################
# read.phe(...)   <EXTERNAL>                                       #
# DESCRIPTION: Reads in the .phe phenotype file, as described in   #
#  the pbat literature.                                            #
# PARAM: filename    filename to open; does not need .phe extension#
#        na.strings  strings to represend NA; '-','.' are defaults #
#        ...         options for read.table, do NOT put in         #
#                     header=TRUE (will cause an error - incorrect!#
#        sym         Toggles symbolic reading (sets attr )         #
# RETURN: dataframe representing the .phe file                     #
####################################################################
read.phe <- function( filename, na.strings=c("-",".","NA"), lowercase=TRUE, sym=TRUE, ... ) {

  #-----------------------------------------------------------------
  # From PBAT_overview.doc (I only have a temporary internet add...-
  # Format of the phenotype file (xbat.phe):                       -
  #  First line:  names of all markers in the sequence of the      -
  #               genotype data                                    -
  #  Remaining lines:  pid id trait_1 trait_2, with pid and id the -
  #                     pedigree and individual ID, respectively.  -
  #    Trait_i refers to the i-th trait value of the individual.   -
  #  Note:  Missing values are coded as '-' or '.' in fbat         -
  #-----------------------------------------------------------------

  # according to other documentation, na.strings can be '-' or '.'
  filename <- str.file.extension( filename, ".phe" );
  if( spaceInFilename(filename) & sym==TRUE )  ## sym==TRUE added 01/14/2008
    stop( spaceInFilenameError(filename) ) ## added 05/17/2007
  phe <- read.badheader( filename, na.strings=na.strings, lowercase=lowercase, onlyHeader=sym, max=-1, ... );

  if( sym ){
    phe2 <- data.frame( matrix( 0, nrow=1, ncol=2+length(phe$header) ) );
    names( phe2 ) <- c( "pid","id", phe$header );
    class( phe2 ) <- 'phe';
    return( set.sym( phe2, filename ) );
  }

  names( phe$table ) <- make.names(  c( "pid","id", phe$header )  );

  if( lowercase )
    names( phe$table ) <- tolower( names( phe$table ) );

  class(phe$table) <- c("phe","data.frame");
  return( phe$table );
} ### VERIFIED ### on the total dataset; unknown's for now



####################################################################
# as.phe(...)  <EXTERNAL>                                          #
# Creates a 'phe' (phenotype) class object for use with pbat routs.#
# PARAM df    dataframe with the data                              #
#       pid   string for the column header for 'pid'               #
#       id    string for the column header for 'id'                #
####################################################################
as.phe <- function( df, pid="pid", id="id" ) {
  if( is.ped(df) )
    stop( "Cannot write out the pedigree file as a phenotype file.  Did you mix up 'phe' and 'ped' objects?" );

  # ensure proper ordering...
  pidCol <- df[pid]; idCol <- df[id];
  df <- dfr.r( df, c(pid,id) );
  #header <- names(df);
  df <- cbind( pidCol, idCol, df );
  names(df)[1:2] <- c("pid","id");

  class(df) <- c("phe","data.frame");

  return(df);
}

####################################################################
# write.phe(...)  <EXTERNAL>                                       #
# Writes a phenotype file out from proper data format.             #
# PARAM file  string representing filename                         #
#       df    dataframe with the data                              #
#       pid   string for the column header for 'pid'               #
#       id    string for the column header for 'id'                #
####################################################################
write.phe <- function( file, phe ) {
  if( !is.phe(phe) )
    warning( "'phe' parameter is not of class phe; errors may occur in writing this object.");

  # tack on extension, if not already there
  if( !is.null(file) && is.character(file) ){
    file <- str.file.extension(file,".phe");
  }

  ##write.badheader( file, phe, names(phe)[-c(1,2)] );
  write.badheader( file, phe, names(phe)[-c(1,2)], na="-" ); ## 01/19/2006 bugfix, 06/14/07 update so fbat friendly!
}

####################################################################
# Thomas Hoffmann                                                  #
# CREATED:     06/07/2005                                          #
# RE-CREATED:  06/20/2005                                          #
# MODIFIED:    05/28/2006                                          #
#                                                                  #
# DESCRIPTION:                                                     #
#  Methods for setting and getting the pbat binary                 #
#   executable filename, and methods for reading/writing           #
#   both .phe/.ped files ( read.badheader )                        #
####################################################################

####################################################################
#                                                                  #
# .phe & .ped FILE READING AND WRITING                             #
#                                                                  #
####################################################################

## New 05/17/07
spaceInFilename <- function( str ) {
  if( is.null(str) || is.na(str) ) return( FALSE );
  return( length( unlist( strsplit( str, " " ) ) ) > 1 )
}
spaceInFilenameError <- function( str )
  return( paste( "The filename specified, '", str, "', contains spaces (i.e. the file itself contains spaces or the pathname contains spaces). This is not supported in P2BAT. Please rename the file or move it to a path that does not contain spaces. This is regretably unfixable.", sep="" ) )

####################################################################
# read.badheader(...)                                              #
# DESCRIPTION: Use when the number of headers does not match the   #
#              data.                                               #
# WARNING: DO NOT do read.badheader( ...., header=T )              #
# PARAM: file         file or filename                             #
#        sep          seperator                                    #
#        lowercase    makes all the headers lowercase (if onlyHeader=F)
#        onlyHeader   for symbolic - only reads in the header      #
#        ...   other parameters to 'read.table'; header=T not sup. #
# RETURN $header  headers read from the file                       #
#        $table   table read from the file (read.table(header=F))  #
####################################################################
read.badheader <- function( file, sep="", lowercase=TRUE, onlyHeader=FALSE, max=100, ... ) {
  # The idea behind this is we read in the header, and then piggy-back
  #  onto the read.table function!
  # Some of this coding was taken from read.table...

  # Open the file (unless it's already a file)
  if (is.character(file)) {
    file <- file(file, "r")
    on.exit(close(file))
  }

  # If it's not a symbolic (onlyHeader=T), then max shouldn't apply
  if( onlyHeader==FALSE ) max <- -1;

  # Read in the header line.
  header <- scan(file, what="", sep=sep,
                 nlines=1, quiet=TRUE, skip=0, strip.white=TRUE,
                 blank.lines.skip=TRUE, comment.char="#", nmax=max );
  if( lowercase & onlyHeader==FALSE )  ## when symbolic, the case cannot be changed; is this the best place to make that change??
    header <- tolower(header); # convenience

  ## If we've reached the max, then it should be pure symbolic (no info whatsoever)
  ##if( max>0 && length(header)==max ){
    ##cat( 'SNPs number >', max, 'so no consistency checks will be done by pbatR on the SNPs (but still on the rest of the formula, and then PBAT will check to some extent later.' );
    ##header <- NULL;
  ##}

  # Read in the rest of the table, prevent that error described above
  #  in the "WARNING" section.
  table <- NULL;
  if( !onlyHeader )
    table <- read.table( file, sep=sep, header=FALSE, ... );

  return( list( header=header, table=table) );
} ### VERIFIED ###

####################################################################
# write.badheader(...)                                             #
# DESCRIPTION: Use when the number of headers does not match the   #
#              data; writes to disk.                               #
# PARAM: file       file (connection)  or filename                 #
#        dataframe  data frame to write to disk                    #
#        header     header to write                                #
# (PARAM):  sep        seperator                                   #
#           col.names  T/F whether to print column names           #
#           row.names  T/F whether to print row names              #
#           ...        other parameters to 'read.table'            #
####################################################################
write.badheader <- function( file, dataframe, header,
                             col.names=FALSE, row.names=FALSE,
                             sep=" ", ... ) {
  # like read.badheader - piggy-back on write.table :)

  # Open the file (unless it's already a file)
  if (is.character(file)) {
    file <- file(file, "w")
    on.exit(close(file))
  }

  # write the header
  cat( header, file=file, sep=sep );
  cat( "\n", file=file );

  # write the table
  write.table( dataframe, file, col.names=col.names, row.names=row.names, sep=sep, quote=FALSE, ... ); ## 05/14/07 update ...
} ## VERIFIED ##


####################################################################
# dfr.r(...)                                                       #
# Removes columns from a dataframe, searching by column names.     #
# PARAM df      data frame                                         #
#       strVec  vector of strings representing column names to     #
#                eliminate                                         #
# RETURNS dataframe with columns removed                           #
####################################################################
#####################################
## data frame remove, using strings #
#####################################
dfr <- function( df, str ) {
  if( sum(names(df)==str) != 1 ) {
    warning("dfr: no headers match, or more than one header matches");
    return( df );
  }
  return( df[-which(names(df)==str)] );
}
dfr.r <- function( df, strVec ) {
  for( i in 1:length(strVec) )
    df <- dfr( df, strVec[i] );

  return( df );
}

as.phe <- function( ... ) {
  return( write.phe(NULL,...) );
}
as.ped <- function( ... ) {
  return( write.ped(NULL,...) );
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

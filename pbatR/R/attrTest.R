####################################################################
# Thomas Hoffmann                                                  #
# CREATED:  05/??/06                                               #
#                                                                  #
# DESCRIPTION:                                                     #
#  Contains attribute setting routines for symbolic loading.       #
#  Contains 'expandPath' to aid with converting relative filenames #
#   to absolute filenames.                                         #
####################################################################

###############################################################
## expandPath(...)                                            #
## Expands the path of a filename (in case the use cwd's).    #
###############################################################
expandPath <- function( filename ){
  if( (strfindf(filename,'/')==-1 & strfindf(filename,'\\')==-1 )
     | substring(filename,1,2)=='..' ) {
    if( !isWindows() ) {
      filename <- paste( getwd(), '/', filename, sep='' );
    }else{
      filename <- paste( getwd(), '\\', filename, sep='' );
    }
  }
  return( filename );
}


###############################################################
## is.sym(...)                                                #
## Returns if an object was set to be symbolic.               #
###############################################################
is.sym <- function( x )
  return( !is.null(attr(x,'sym')) );

###############################################################
## set.sym(...)                                               #
## Sets an object to be symbolic.                             #
## WARNING: you must do                                       #
##   x <- set.sym(x)                                          #
###############################################################
set.sym <- function( x, filename, expand.path=TRUE ){
  if( expand.path )
    filename <- expandPath( filename );
  attr(x,'sym') <- filename;
  return( x );
}

###############################################################
## get.sym(...)                                               #
## Returns the filename of the symbolic object.               #
###############################################################
get.sym <- function( x )
  return( attr(x,'sym') )

## Leftover debugging
#df <- data.frame( c(0,1), c(3,4) )
#is.sym( df )  ## should be false
#df <- set.sym( df, 'attrTest.R' )
#is.sym( df )  ## Now should be true
#get.sym( df ) ## Should return the value
#
#df <- set.sym( df, '/home/notexist.R' )
#get.sym( df ) ## Should return the value
#
#df <- set.sym( df, '../notexist.R' )
#get.sym( df ) ## Should return the value
#
#df <- set.sym( df, '/home/../notexist.R' )  ## This is malformed to some extent, but still works just fine!
#get.sym( df ) ## Should return the value



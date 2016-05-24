####################################################################
# Thomas Hoffmann                                                  #
# CREATED:  some time ago                                          #
# MODIFIED: 01/25/2006                                             #
#                                                                  #
# DESCRIPTION:                                                     #
#  Contains string functions.                                      #
####################################################################

##################################################
# I can't find the right string functions in R;  #
#  so I just created a few real quick to aid us. #
##################################################

strlen <- function( s ) { return(nchar(s)) };

strfindf <- function( s, p, start=1 ) {
  if( length(s)!=1 || length(p)!=1 )
    return(-1) ## 08/11/20008

  if( strlen(s)<strlen(p) )
    return(-1); ## 01/25/2006
  ls <- strlen(s);  lp <- strlen(p);
  if( ls < lp ) return(-1);
  if( ls-lp+1<start ) return(-1); ## 01/26/2006
  for( i in start:(ls-lp+1) ){
    if( p==substring(s,i,i+lp-1) ) #@$%!!!
      return(i);
  }
  return(-1); # not found
}

# 10/06/2005
# find - from behind!
strfindb <- function( s, p, start=strlen(s)-strlen(p) ) {
  if( length(s)!=1 || length(p)!=1 )
    return(-1) ## 08/11/2008

  if( strlen(s)<strlen(p) )
    return(-1); ## 01/25/2006
  lp <- strlen(p);
  for( i in seq(start,1,by=-1) ) {
    if( p==substring(s,i,i+lp-1) ) #@$%!!!
      return(i);
  }
  return(-1); # not found
}
# get's the path of a file / executable
str.getpath <- function( s ) {
  # two delimiters
  loc1 <- strfindb(s,"/");
  loc2 <- strfindb(s,"\\");
  if( loc2 > loc1 )
    loc1=loc2;

  # make sure it's a path (and not just cwd)
  if( loc1 < 2 )
    return("");

  # return the substring
  return( substring( s, 1, loc1-1 ) );
}


str.extract.after <- function( s, p='~' ) {
  loc <- strfindf( s,p );
  if( loc == -1 ) return("");
  return( substring( s, loc+strlen(p) ) );
}

str.extract.afterb <- function( s, p='~' ) {
  loc <- strfindb( s,p );
  if( loc == -1 ) return("");
  return( substring( s, loc+strlen(p) ) );
}

str.file.extension <- function( s, extension='txt' ){
  if( substring(extension,1,1)!="." ) #@$%!!!
    extension = paste(".",extension,sep="");

  if( strlen(s) > strlen(extension) ){
    if( substring(s,strlen(s)-strlen(extension)+1) == extension )
      return( s );
  }
  return( paste(s,extension,sep="") );
}

######################################################
## strReplace(...)                                   #
## Replaces occurences of one string with another.   #
## HOWEVER - if it's at the beginning or the end it  #
##  will get lost.                                   #
## PARAM  str   string to replace                    #
##        find  string to look for                   #
##        repl  string to replace with               #
######################################################
strReplace <- function( str, find, repl )  {
  strList <- strsplit( str, find )[[1]];

  if( length(strList) <= 1 )
    return(str);

  strRet <- strList[1];
  ##print( strRet );
  for( i in 2:length(strList) ){
    strRet <- paste( strRet, strList[i], sep=repl );
    ##print( strRet );
  }

  return(strRet);
}


#substring( string, start, stop=... )
#substr( string, start, stop )
# WARNING: substr() does not exist in s-plus, but it does the same
#  thing as substring() with the third parameter.. strange R;
#  probably for backwards compatibility??

catn <- function( ..., file="" ) {
  cat( ..., file=file );
  cat( '\n', file=file );
}


## Additions 09/07/2006

## returns the extension of a file, _assuming_ that there is one
file.extension <- function( filename ){
  loc <- strfindb( filename, "." );
  if( loc!=-1 ) {
    ## then it was found, so extract and return
    return( substring( filename, first=loc+1 ) );
  }
  return( "" );
} ## Debugged
#file.extension( "badger.pped" );
#file.extension( "buck.badger.pped" );
#file.extension( "badger.ped" );
#file.extension( "badger" );

file.strip.extension <- function( filename ){
  loc <- strfindb( filename, "." );
  if( loc!=-1 ){
    ## then it was found, so extract and return
    return( substring( filename, first=1, last=loc-1 ) );
  }
  return( "" );
} ## Debugged
#file.strip.extension( "badger.pped" );
#file.strip.extension( "buck.badger.pped" );
#file.strip.extension( "badger.ped" );
#file.strip.extension( "badger" );

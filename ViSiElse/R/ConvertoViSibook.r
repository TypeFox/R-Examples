#' ConvertoViSibook convert a \code{data.frame} in \code{ViSibook} object.
#' 
#' @title Function \code{ConvertoViSibook}
#' @rdname ConvertoViSibook
#' @aliases ConvertoViSibook
#' @export ConvertoViSibook
#' @param x a dataframe. \code{x} should contains at least the columns \strong{vars, label, typeA, showorder, deb, fin }.
#' Optionally other characteristics can be filled :
#'  \strong{GZDebn,  GZFin, Repetition, BZBeforeDeb, BZBeforeFin, BZAfterDeb, BZAfterFin, BZLong , BZLtype }.
#' @return a ViSibook object.
#' @seealso See \code{\linkS4class{ViSibook}} to get the definitions of the columns 
#' and see  \code{\link{plot-ViSigrid-method}} for examples.
ConvertoViSibook <- function(x ) {  
  if (is.na( match( "vars" , colnames( x ) ) ) ) { stop( " ConvertoViSibook : colname \"vars\" not found \n " ) }
  if (is.na( match( "label" , colnames( x ) ) ) ) { stop( " ConvertoViSibook : colname \"label\" not found \n " ) }
  if (is.na( match( "typeA" , colnames( x ) ) ) ) { stop( " ConvertoViSibook : colname \"typeA\" not found \n " ) }
  if (is.na( match( "showorder" , colnames( x ) ) ) ) { stop( " ConvertoViSibook : colnames \"showorder\" not found " ) }
  if (is.na( match( "deb" , colnames( x ) ) ) ) { stop( " ConvertoViSibook : colnames \"deb\" not found " ) }
  if (is.na( match( "fin" , colnames( x ) ) ) ) { stop( " ConvertoViSibook : colnames \"fin\" not found " ) }
  return( ViSibook( 
    vars = x[ , match( "vars" , colnames( x ) ) ] ,
    label = x[ , match( "label" , colnames( x ) ) ] ,
    typeA = x[ , match( "typeA" , colnames( x ) ) ] ,
    showorder = x[ , match( "showorder" , colnames( x ) ) ] ,
    deb = x[ , match( "deb" , colnames( x ) ) ] ,
    fin = x[ , match( "fin" , colnames( x ) ) ] ,
    GZDeb = switch( as.character( is.na( match( "GZDeb" , colnames( x ) ) ) ) , 
                    "TRUE" = vector() , 
                    "FALSE" = x[ , match( "GZDeb" ,colnames( x ) ) ] ) ,
    GZFin = switch( as.character( is.na( match( "GZFin" , colnames( x ) ) ) ) , "TRUE" = vector() , "FALSE" = x[ , match( "GZFin" , colnames( x ) ) ] ) ,
    Repetition = switch( as.character( is.na( match( "Repetition" , colnames( x ) ) ) ) , "TRUE" = vector() , "FALSE" = x[ , match( "Repetition" , colnames( x ) ) ] ) ,
    BZBeforeDeb = switch( as.character( is.na( match( "BZBeforeDeb" , colnames( x ) ) ) ) , "TRUE" = vector() , "FALSE" = x[ , match( "BZBeforeDeb" , colnames( x ) ) ] ),
    BZBeforeFin = switch( as.character( is.na( match( "BZBeforeFin", colnames( x ) ) ) ) , "TRUE" = vector() , "FALSE" = x[ ,match( "BZBeforeFin" , colnames( x ) ) ] ) ,
    BZAfterDeb = switch( as.character(is.na( match( "BZAfterDeb" , colnames( x ) ) ) ) , "TRUE" = vector() , "FALSE" = x[ , match( "BZAfterDeb" , colnames( x) ) ] ) ,
    BZAfterFin = switch( as.character(is.na( match( "BZAfterFin" , colnames( x ) ) ) ) , "TRUE" = vector() , "FALSE" = x[ , match( "BZAfterFin" , colnames( x ) ) ] ) ,
    BZLong = switch( as.character(is.na( match( "BZLong" , colnames( x ) ) ) ) , "TRUE" = vector() , "FALSE" = x[ , match( "BZLong" , colnames( x) ) ] ) ,
    BZLtype = switch( as.character(is.na( match( "BZLtype" , colnames( x ) ) ) ) , "TRUE" = vector() , "FALSE" = x[ , match( "BZLtype" , colnames( x ) ) ] ) ,
    NAMES = colnames(x)
  ) )
}

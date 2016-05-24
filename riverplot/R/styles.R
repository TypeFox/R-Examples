#' @rdname riverplot-styles
#' @export 
default.style <- function( ) {

    ret <- list(
      nodestyle= "regular",
      col=       "grey",
      srt=       "90",
      lty=       1,
      textcol=   "black",
      edgecol=   "gradient",
      edgestyle= "sin"
      )
  
  class( ret ) <- c( class( ret ), "riverplotStyle" )
  return( ret )
}

#' @rdname riverplot-styles
#' @export 
updateRiverplotStyle <- function( style, master ) getstyle( style, master )

# function for updating styles. s is filled up with default values if these
# values are empty. If update.missing is TRUE, update also these fields
# which are missing from the global default style.
getstyle   <- function( s, defaults= NULL, update.missing= FALSE ) {

  if( is.null( s ) ) s <- list( )
  
  if( is.null( defaults ) ) defaults <- default.style()

  for( n in names( defaults ) ) {
    if( is.null( s[[n]] ) ) s[[n]] <- defaults[[n]]
  }

  if( update.missing ) {
    defaults <- default.style()
    for( n in names( defaults ) ) {
      if( is.null( s[[n]] ) ) s[[n]] <- defaults[[n]]
    }
  }

  class( s ) <- c( class( s ), "riverplotStyle" )

  return( s )
}


# checks whether attr for id in styles is equal to value
isStyle <- function( styles, id, attr, value ) {

  if( is.null( styles ) ) return( FALSE )
  if( is.null( styles[[id]] ) ) return( FALSE )
  if( is.null( styles[[id]][[attr]] ) ) return( FALSE )
  if( styles[[id]][[attr]] %in% value ) return( TRUE ) 
  #printf( "isStyle: %s", styles[[id]][[attr]] )

  FALSE
}

# copy attribute from id.from to id.to
copyattr <- function( styles, id.from, id.to, attr ) {
  
  val <- getattr( styles, id.from, attr )
  styles <- setattr( styles, id.to, attr, val )
  return( styles )

}

setattr <- function( styles, id, attr, value ) {

  if( is.null( styles ) ) styles <- list()
  if( is.null( styles[[id]] ) ) styles[[id]] <- list()
  styles[[id]][[attr]] <- value

  return( styles )
}

# return attribute for id in styles. If NULL, return the default
getattr <- function( styles, id, attr ) {

  def <- TRUE

  if( is.null( styles ) || 
      is.null( styles[[id]] ) ||
      is.null( styles[[id]][[attr]] ) ) 
    tmp <- default.style()
  else
    tmp <- styles[[id]]

  return( tmp[[attr]] )
}

# merges styles s1 and s2, overwriting s1 if IDs are repeated
# if s1 is NULL, it will be created.
# if s2 is NULL, it will be ignored
mergestyles <- function( s1, s2 ) {

  if( is.null( s1 ) ) s1 <- list()
  if( ! is.null( s2 ) ) {
    for( n in names( s2 ) ) {
      s1[[n]] <- s2[[n]]
    }
  }

  return( s1 )
}

# As of 2.8.1, we are still not able to write S4 methods for rbind(). Our
# original plan was to use combine() as in biobase, but this seems to be causing
# clashes for the generic function when both packages are loaded. Improved S4
# support for functions whose first argument is ... seems to be on the way, so
# we'll just do S3 for now, and move up to S4 once it's implemented. Note that
# as of 2.8.1, rbind.data.frame() still exists and is used for data frames.

# Update! S3 method dispatch for rbind() is non-standard (see its documentation)
# and it produced unexpected dispatch to the matrix method when presented with a
# mix of Intervals and Intervals_full objects. As a consequence, we switched to
# c(), which uses standard S3 dispatch.

c.Intervals <- function( ... ) {
  args <- list(...)
  # Drop NULL arguments
  if ( any( sapply( args, is.null ) ) )
    args <- args[ !sapply( args, is.null ) ]
  # Check if we should just return a list
  classes <- sapply( args, class )
  if ( !all( classes %in% c( "Intervals", "Intervals_full" ) ) )
    return( list( ... ) )
  same_class <- all( classes == "Intervals" )
  # We are in fact dealing with intervals only...
  if ( !all( sapply( args, type ) == type( args[[1]] ) ) )
    stop( "All arguments should have the same 'type' slot." )
  # Check for common closure
  same_closed <- all( sapply( args[-1], function(x) identical( closed( args[[1]] ), closed( x ) ) ) )
  # Coerce up if necessary
  if ( !same_class || ( type( args[[1]] ) == "R" & !same_closed ) ) {
    warning( "Coercion to 'Intervals_full' required.", call. = FALSE )
    return( do.call( c, lapply( args, as, "Intervals_full" ) ) )
  }
  # Convert to common closure for Z
  if ( type( args[[1]] ) == "Z" & !same_closed )
    args <- lapply(
                   args,
                   adjust_closure,
                   close_left = closed( args[[1]] )[1],
                   close_right = closed( args[[1]] )[2]
                   )
  result <- args[[1]]
  result@.Data <- do.call( rbind, lapply( args, function(x) x@.Data ) )
  return( result )
}

c.Intervals_full <- function( ... ) {
  args <- list(...)
  if ( any( sapply( args, is.null ) ) )
    args <- args[ !sapply( args, is.null ) ]
  classes <- sapply( args, class )
  if ( !all( classes %in% c( "Intervals", "Intervals_full" ) ) )
    return( list( ... ) )
  if ( !all( sapply( args, type ) == type( args[[1]] ) ) )
    stop( "All arguments should have the same 'type' slot." )
  if ( !all( classes == "Intervals_full" ) ) {
    warning( "Coercion to 'Intervals_full' required.", call. = FALSE )
    args <- lapply( args, as, "Intervals_full" )
  }
  result <- args[[1]]
  result@.Data <- do.call( rbind, lapply( args, function(y) y@.Data ) )
  closed(result) <- do.call( rbind, lapply( args, closed ) )
  return(result)
}

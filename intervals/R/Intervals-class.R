# We define two classes for two-column interval endpoint matrices. The basic class
# has a two-element boolean vector indicating whether endpoints are closed or
# not. The full class caries a matrix with one boolean per endpoint, permitting
# full control.




######## Class definitions

#### Virtual

# (R 2.6.2) If I add the "VIRTUAL" tag to the representation, I am not able to
# extend this class! I believe this arises because I am already extending a base
# class, but I am not sure. This tag would be more appropriate, but I leave it
# off... 

setClass(
         "Intervals_virtual",
         representation( type = "character" ),
         prototype(
                   matrix( 0, 0, 2 ),
                   type = "R"
                   ),
         contains = "matrix",
         validity = function( object ) {
           # Check main matrix
           if ( !is.double( object@.Data ) || ncol( object@.Data ) != 2  )
             return( "The 'Intervals' classes are based on two-column, numeric matrices." )
           # Check 'type' string
           if ( length( object@type ) != 1 || !( object@type %in% c( "Z", "R" ) ) )
             return( "The 'type' slot should be 'Z' or 'R'." )
           # For type 'Z', check for integral endpoints
           if ( object@type == "Z" && !all( object@.Data[ is.finite( object@.Data ) ] %% 1 == 0 ) )
             return( "Non-integer-valued endpoints not permitted for type 'Z'." )
           # Check for valid intervals
           if ( any( object@.Data[,2] < object@.Data[,1], na.rm = TRUE ) )
             return( "One or more intervals with second endpoint before first." )
           return( TRUE )
         }
         )

setMethod(
          "initialize",
          signature( "Intervals_virtual" ),
          function( .Object, .Data, ... ) {
            if ( missing( .Data ) ) callNextMethod( .Object, ... )
            else {
              if ( is.data.frame( .Data ) )
                .Data <- as.matrix( .Data )
              if ( !is.matrix( .Data ) )
                .Data <- matrix( .Data, ncol = 2 )
              if ( is.integer( .Data ) ) {
                # warning( "Converting endpoints from 'integer' to 'numeric' data type. See class documentation.", call. = FALSE )
                .Data <- matrix( as.numeric( .Data ), nrow( .Data ), ncol( .Data ) )
              }
              callNextMethod( .Object, .Data, ... )
            }
          }
          )

#### Intervals

# Common endpoint closure state for all intervals

setClass(
         "Intervals",
         representation( closed = "logical" ),
         prototype( closed = c( TRUE, TRUE ) ),
         contains = "Intervals_virtual",
         validity = function( object ) {
           # Check 'closed' slot
           if ( length( object@closed ) != 2 || any( is.na( object@closed ) ) )
             return( "The 'closed' slot should be a logical of length 2. NA values are not permitted." )
           return( TRUE )
         }
         )

setMethod(
          "initialize",
          signature( "Intervals" ),
          function( .Object, .Data, closed, ... ) {
            if ( missing( .Data ) )
              callNextMethod( .Object, ... )
            else {
              if ( missing( closed ) )
                callNextMethod( .Object, .Data, ... )
              else {
                if ( length( closed ) == 1 ) closed <- c( closed, closed )
                callNextMethod( .Object, .Data, closed = closed, ... )
              }
            }
          }
          )

#### Intervals_full

# Full control of endpoint closure. Note that if the 'closed' slot is omitted,
# we use an 'initialize' method to create an appropriately sized matrix of TRUE
# values. We also permit vector input, with recycling, for the 'closed' slot.

setClass(
         "Intervals_full",
         representation( closed = "matrix" ),
         prototype( closed = matrix( TRUE, 0, 2 ) ),
         contains = "Intervals_virtual",
         validity = function( object ) {
           # Check 'closed' slot
           if (
               !is.logical( object@closed ) ||
               dim( object@.Data ) != dim( object@closed ) ||
               any( is.na( object@closed ) )
               )
             return( "The 'closed' slot should be a logical matrix with the same dimensions as the main endpoints matrix. NA values are not permitted." )
           return( TRUE )
         }
         )

setMethod(
          "initialize",
          signature( "Intervals_full" ),
          function( .Object, .Data, closed, ... ) {
            if ( missing( .Data ) )
              callNextMethod( .Object, ... )
            else {
              if ( !is.matrix( .Data ) )
                .Data <- matrix( .Data, ncol = 2 )
              if ( missing( closed ) )
                closed <- matrix( TRUE, nrow( .Data ), 2 )
              if ( is.vector( closed ) ) {
                if ( length( closed ) > 2 )
                  stop( "The 'closed' argument should be a matrix, or a vector of length 1 or 2." )
                closed <- matrix(
                                 if ( nrow( .Data ) == 0 ) logical() else closed,
                                 nrow( .Data ),
                                 2,
                                 byrow = TRUE
                                 )
              }
              callNextMethod( .Object, .Data, closed = closed, ... )
            }
          }
          )




######## Constructor functions

Intervals <- function( ... )
  new( "Intervals", ... )

Intervals_full <- function( ... )
  new( "Intervals_full", ... )




######## Subsetting

setMethod(
          "[",
          signature( "Intervals" ),
          function( x, i, j, ..., drop ) {
            if ( missing(i) ) i <- rep( TRUE, nrow(x) )
            if ( missing(j) ) {
              # Preserve class. Note that both [i,] and [i] syntax subset rows.
              if ( any( is.na( i ) ) )
                warning( "NA indices encountered.", call. = FALSE )
              x@.Data <- x@.Data[i,,drop=FALSE]
              return( x )
            }
            else return( x@.Data[i,j] )
          }
          )

# Note: row name handling for matices is completely inadequate for our purposes
# here. Lots of obviously desirable behavior (like setting rownames(x)[i] when x
# doesn't yet have rownames) tends to produce errors.

setMethod(
          "[<-",
          signature( x = "Intervals", i = "ANY", j = "missing", value = "Intervals_virtual" ),
          function( x, i, j, value ) {
            #### Error checking
            if ( type(x) != type(value) )
              stop( "Types do not match (Z vs. R)." )
            if ( is.character(i) ) i <- match( i, rownames( x ) )
            if ( any( is.na( i ) ) )
              stop( "Cannot assign to NA indices or row names which do not exist." )
            n <- length( (1:nrow(x))[i] )
            if ( n != nrow( value ) )
              stop( "Replacement object is of the wrong size." )
            #### Coerce up?
            coerce_x <- FALSE
            if ( is( value, "Intervals_full" ) ) coerce_x <- TRUE
            else {
              matches <- all( closed(x) == closed(value) )
              if ( !matches ) {
                if ( type(x) == "Z" ) value <- adjust_closure( value, closed(x)[1], closed(x)[2] )
                else coerce_x <- TRUE
              }
            }
            if ( coerce_x ) {
              warning( "Coercion to 'Intervals_full' required.", call. = FALSE )
              x <- as( x, "Intervals_full" )
              x[ i, ] <- value
              return(x)
            }            
            #### Data
            x@.Data[i,] <- value@.Data
            #### Rownames
            has_names_x <- !is.null( rownames(x) )
            has_names_value <- !is.null( rownames(value) )
            if ( has_names_x ) {
              if ( has_names_value ) rownames(x)[i] <- rownames(value)
              else rownames(x)[i] <- ""
            }
            else {
              if ( has_names_value ) {
                rownames(x) <- rep( "", nrow(x) )
                rownames(x)[i] <- rownames(value)
              }
            }
            return(x)
          }
          )

setMethod(
          "[",
          signature( "Intervals_full" ),
          function( x, i, j, ..., drop ) {
            if ( missing(i) ) i <- rep( TRUE, nrow(x) )
            if ( missing(j) ) {
              # Preserve class. Note that both [i,] and [i] syntax subset rows.
              if ( is.character(i) ) i <- match( i, rownames( x ) )
              if ( any( is.na( i ) ) )
                warning( "NA indices encountered.", call. = FALSE )
              x@.Data <- x@.Data[i,,drop=FALSE]
              x@closed <- x@closed[i,,drop=FALSE]
              # We may have NAs in closed if present in i, if any(is.na(i)) == TRUE
              x@closed[ is.na(x@closed) ] <- TRUE
              return( x )
            }
            else return( x@.Data[i,j] )
          }
          )

setMethod(
          "[<-",
          signature( x = "Intervals_full", i = "ANY", j = "missing", value = "Intervals_virtual" ),
          function( x, i, j, value ) {
            #### Error checking
            if ( type(x) != type(value) )
              stop( "Types do not match (Z vs. R)." )
            if ( is.character(i) ) i <- match( i, rownames( x ) )
            if ( any( is.na( i ) ) )
              stop( "Cannot assign to NA indices or row names which do not exist." )
            n <- length( (1:nrow(x))[i] )
            if ( n != nrow( value ) )
              stop( "Replacement object is of the wrong size." )            
            #### Data
            x@.Data[i,] <- value@.Data
            if ( is( value, "Intervals" ) )
              x@closed[i,] <- matrix( value@closed, n, 2, byrow = TRUE )
            else
              x@closed[i,] <- value@closed
            #### Rownames
            has_names_x <- !is.null( rownames(x) )
            has_names_value <- !is.null( rownames(value) )
            if ( has_names_x ) {
              if ( has_names_value ) rownames(x)[i] <- rownames(value)
              else rownames(x)[i] <- ""
            }
            else {
              if ( has_names_value ) {
                rownames(x) <- rep( "", nrow(x) )
                rownames(x)[i] <- rownames(value)
              }
            }
            return(x)
          }
          )




######## Coercion

setMethod(
          "coerce",
          signature( from = "Intervals", to = "Intervals_full" ),
          function( from, to, strict ) {
            new(
                "Intervals_full",
                from@.Data,
                type = type( from ),
                closed = cbind(
                  rep( closed(from)[1], nrow( from ) ),
                  rep( closed(from)[2], nrow( from ) )
                  )
                )
          }
          )

setMethod(
          "coerce",
          signature( from = "Intervals_full", to = "Intervals" ),
          function( from, to, strict ) {
            if ( nrow( from ) == 0 ) new_closed <- rep( TRUE, 2 )
            else new_closed <- closed( from )[1,]
            if ( !all( t( closed( from ) ) == new_closed ) )
              stop( "Intervals do not all have the same endpoint closure." )
            new(
                "Intervals",
                from@.Data,
                type = type( from ),
                closed = new_closed
                )
          }
          )

setMethod(
          "coerce",
          signature( from = "Intervals_virtual", to = "character" ),
          function( from, to, strict ) {
            if ( nrow( from ) == 0 )
              return( character() )
            else {
              cl <- closed( from )
              # So we only write main code once
              if ( is( from, "Intervals" ) )
                cl <- matrix( cl, nrow(from), 2, byrow = TRUE )
              result <- paste(
                              ifelse( cl[,1], "[", "(" ),
                              from[,1], ", ", from[,2],
                              ifelse( cl[,2], "]", ")" ),
                              sep = ""
                              )
              names( result ) <- rownames( from )
              return( result )
            }
          }
          )




######## Union

setClassUnion( "Intervals_virtual_or_numeric", c( "Intervals_virtual", "numeric" ) )

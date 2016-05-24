setGeneric( "interval_complement", def = function(x, ...) standardGeneric( "interval_complement" ) )

setMethod(
          "interval_complement",
          signature( "Intervals_virtual" ),
          function(x, check_valid = TRUE) {
            # Sort and clean up
            x <- reduce( x, check_valid )
            # When the data type of the endpoints matrix is integer,
            # complications arise from maximum/minimum representable integer
            # values. For Intervals objects, the complement will often
            # need to include the minimal/maximal represenable integer, but one
            # or both endpoints may be open. In such cases, we adjust close with
            # a warning. For the moment, we force numeric endpoints throughout
            # the package.
            endpoints <-
              if ( nrow(x) == 0 )
                matrix( c( -Inf, Inf ), 1 )
              else
                rbind(
                      if ( is.finite( x[1,1] ) ) c( -Inf, x[1,1] ) else NULL,
                      cbind( x[-nrow(x),2], x[-1,1] ),
                      if ( is.finite( x[nrow(x),2] ) ) c( x[nrow(x),2], Inf ) else NULL
                      )
            closed <-
              if ( class(x) == "Intervals" )
                # Note that we ignore closure for non-finite endpoints.
                !closed(x)[2:1]
              else
                if ( nrow(x) == 0 ) TRUE                  
                else                    
                  !rbind(
                         if ( is.finite( x[1,1] ) ) c( TRUE, closed(x)[1,1] ) else NULL,
                         cbind( closed(x)[-nrow(x),2], closed(x)[-1,1] ),
                         if ( is.finite( x[nrow(x),2] ) ) c( closed(x)[nrow(x),2], TRUE ) else NULL
                         )
            new(
                class(x),
                endpoints,
                type = type(x),
                closed = closed
                )
          }
          )

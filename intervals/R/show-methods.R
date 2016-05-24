setMethod(
          "show",
          signature( "Intervals_virtual" ),
          function( object ) {
            cat(
                "Object of class ",
                class( object ),
                "\n",
                nrow( object ),
                " interval",
                ifelse( nrow( object ) == 1, "", "s" ),
                " over ",
                type(object),
                ":\n",
                sep = ""
                )
            ints <- as( object, "character")
            if ( !is.null( rownames( object ) ) ) {
              fmt <- sprintf( "%%%is", max( nchar( rownames( object ) ) ) )
              ints <- paste( sprintf( fmt, rownames( object ) ), ints )
            }
            cat( ints, sep = "\n" )
          }
          )

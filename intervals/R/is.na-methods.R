setMethod(
          "is.na",
          signature( x = "Intervals_virtual" ),
          function(x) is.na( x[,1] ) | is.na( x[,2] )
          )

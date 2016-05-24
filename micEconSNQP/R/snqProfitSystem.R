## -------- system formulas -------------------
snqProfitSystem <- function( nNetput, nFix, form = 0, profit = FALSE ) {
   system <- list()
   for(i in 1:nNetput) {
      system[[i]] <- paste("q",as.character(i)," ~",sep="")
      for(j in 1:nNetput) {
         if(j==1) {
            system[[i]] <- paste(system[[i]]," pr", as.character(j),sep="")
         } else {
            system[[i]] <- paste(system[[i]]," + pr", as.character(j),sep="")
         }
      }
      for(j in 1:nNetput) {
         for(k in 1:nNetput) {
            system[[i]] <- paste(system[[i]], " + pq", as.character(i), ".",
                                 as.character(j), ".", as.character(k), sep="")
         }
      }
      if( nFix > 0 ) {
         for(j in 1:nFix) {
            system[[i]] <- paste(system[[i]]," + f", as.character(j),sep="")
         }
         for(j in 1:nFix) {
            for(k in 1:nFix) {
               system[[i]] <- paste(system[[i]], " + fq", as.character(i), ".",
                                    as.character(j), ".", as.character(k), sep="")
            }
         }
      }
      system[[i]] <- as.formula(system[[i]])
   }
   if( profit ) {
      system[[ nNetput  + 1 ]] <- "profit ~ -1"
      for( j in 1:nNetput ) {
         system[[ nNetput  + 1 ]] <- paste( system[[ nNetput  + 1 ]],
            " + tp", j, sep = "" )
      }
      for( j in 1:nNetput ) {
         for( k in 1:nNetput ) {
            system[[ nNetput  + 1 ]] <- paste( system[[ nNetput  + 1 ]],
               " + tpq", j, ".", k, sep = "" )
         }
      }
      if( nFix > 0 ) {
         for( j in 1:nNetput ) {
            for( k in 1:nFix ) {
               system[[ nNetput  + 1 ]] <- paste( system[[ nNetput  + 1 ]],
                  " + tpf", j, ".", k, sep = "" )
            }
         }
         if( form == 0 ) {
            for( j in 1:nFix ) {
               for( k in 1:nFix ) {
                  system[[ nNetput  + 1 ]] <- paste( system[[ nNetput  + 1 ]],
                     " + tfq", j, ".", k, sep = "" )
               }
            }
         } else {
            for( j in 1:nNetput ) {
               for( k in 1:nFix ) {
                  for( l in 1:nFix ) {
                     system[[ nNetput  + 1 ]] <- paste( system[[ nNetput  + 1 ]],
                        " + tfq", j, ".", k, ".", l, sep = "" )
                  }
               }
            }
         }
      }
      system[[ nNetput  + 1 ]] <- as.formula(system[[ nNetput  + 1 ]])
   }
   return( system )
}

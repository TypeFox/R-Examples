writeMps <- function( file, cvec, bvec, Amat, name="LP" ) {
   nCon <- length(bvec)
   nVar <- length(cvec)

   if( is.null( names(bvec) ) ) {
      blab <- rep("",nCon)
      for(i in 1:nCon) {
         blab[i] <- paste("R_",as.character(i))
      }
   } else {
      blab <- names(bvec)
      for(i in 1:nCon) {
         blab[i] <- gsub(" ","",blab[i])
         if( nchar( blab[i] ) > 8 ) {
            blab[i] <- substr( blab[i], 1, 8 )
         }
         j <- 2
         while( i>1 & blab[i] %in% blab[1:(i-1)])  {
            blab[i] <- paste( substr(blab[i], 1, 7-nchar(as.character(j))),
                              "_", as.character(j), sep="" )
            j <- j+1
         }
      }
   }

   if( is.null( names(cvec) ) ) {
      clab <- rep("",nVar)
      for(i in 1:nVar) {
         clab[i] <- paste("C_",as.character(i))
      }
   } else {
      clab <- names(cvec)
      for(i in 1:nVar) {
         clab[i] <- gsub(" ","",clab[i])
         if( nchar( clab[i] ) > 8 ) {
            clab[i] <- substr( clab[i], 1, 8 )
         }
         j <- 2
         while( i>1 & clab[i] %in% clab[1:(i-1)])  {
            clab[i] <- paste( substr(clab[i], 1, 7-nchar(as.character(j))),
                              "_", as.character(j), sep="" )
            j <- j+1
         }
      }
   }

   write( paste("NAME          ",name,sep=""), file )

   write( "ROWS", file, append=TRUE )
   write( " N  obj", file, append=TRUE )
   for(i in 1:nCon) {
      write( paste(" L  ",blab[i],sep="" ), file, append=TRUE )
   }

   write( "COLUMNS", file, append=TRUE )
   for(i in 1:nVar) {
      line <- paste("    ",clab[i], sep="" )
      line <- paste( line, paste( rep( " ", 14-nchar(line)), collapse=""), "obj", sep="")
      line <- paste( line, paste( rep( " ", 36-nchar(line) - nchar(
         as.character( signif( cvec[i], 10 )))), collapse=""),
         as.character( signif( cvec[i], 10)), sep="")
      write( line, file, append=TRUE )
      for(j in 1:nCon) {
         if( Amat[j,i] != 0 ) {
            line <- paste("    ",clab[i], sep="" )
            line <- paste( line, paste( rep( " ", 14-nchar(line)),collapse=""),
               blab[j], sep="")
            line <- paste( line, paste( rep( " ", 36 - nchar(line) - nchar(
                      as.character( signif( Amat[j,i], 10 )))), collapse=""),
                      as.character( signif( Amat[j,i], 10 )), sep="")
            write( line, file, append=TRUE )
         }
      }
   }

   write( "RHS", file, append=TRUE )
   for(i in 1:nCon) {
      line <- paste("    RHS       ",blab[i], sep="" )
      line <- paste( line, paste( rep( " ", 36-nchar(line) - nchar(
               as.character( signif( bvec[i], 10 )))), collapse=""),
               as.character( signif( bvec[i], 10 )), sep="")
      write( line, file, append=TRUE )
   }

   write( "ENDATA", file, append=TRUE )
}

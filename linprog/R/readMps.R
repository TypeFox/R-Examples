readMps <- function( file, solve=FALSE, maximum=FALSE ) {

   mps <- readLines(file)
   i <- 1

   ## Name
   while( substr( mps[i], 1, 4 ) != "NAME" & i < length(mps) ) {
      i <- i + 1
   }
   if( substr( mps[i], 1, 4 ) == "NAME" ) {
      name <- substr( mps[i], 15, nchar( mps[i] ) )
   } else {
      stop( "MPS file must have a line starting with 'NAME'" )
   }

   ## Rows / Constraints
   while( substr( mps[i], 1, 4 ) != "ROWS" & i < length(mps) ) {
      i <- i + 1
   }
   if( substr( mps[i], 1, 4 ) != "ROWS" ) stop( "MPS file must have a line starting with 'ROWS'" )

   objname <- NULL
   bvec <- NULL   # constraints
   svec <- NULL   # sign of the constraints

   i <- i + 1
   while( substr( mps[i], 1, 7 ) != "COLUMNS" & i < length(mps)) {
      sign <- substr( mps[i], 2,2 )
      if( sign == "E" ) stop("Equaltity constraints are not implemented yet.")
      if( !( sign %in% c("E", "L", "G", "N" ) ) ) {
         sign <- substr( mps[i], 3,3 )
      }
      if( sign %in% c("E", "L", "G" ) ) {
         rname <- strsplit( mps[i], " " )[[1]][length( strsplit( mps[i], " " )[[1]] )]
         bvec <- c( bvec, 0 )
         svec <- c( svec, sign )
         names(bvec)[length( bvec ) ]  <-  rname
         names(svec)[length( svec ) ]  <-  rname
      } else {
         if( sign == "N" ) {
            if( is.null( objname ) ) {
               objname <- strsplit( mps[i], " " )[[1]][length( strsplit( mps[i], " " )[[1]] )]
            }
         } else {
            stop("the 2nd or 3rd column of the rows section must be 'N', 'E', 'L' or 'G'")
         }
      }
      i <- i + 1
   }

   ## Columns
   if( substr( mps[i], 1, 7 ) != "COLUMNS" )
      stop( "MPS file must have a line starting with 'COLUMS'" )
   cvec <- NULL
   Amat <- matrix(0, length(bvec), 0 )
   rownames(Amat) <- names(bvec)
   i <- i + 1
   while( substr( mps[i], 1, 3 ) != "RHS" & i < length(mps)) {
      temp <- strsplit( mps[i], " " )[[1]]
      temp <- temp[ temp != "" ]
      if( !(temp[1] %in% colnames(Amat) ) ) {
         cvec <- c( cvec, 0 )
         Amat <- cbind( Amat, rep( 0, nrow(Amat) ) )
         names(cvec)[length(cvec)]  <- temp[1]
         colnames(Amat)[ncol(Amat)] <- temp[1]
      }
      for( j in 1:((length(temp)-1)/2) ) {
         if( temp[ 2*j ] == objname ) {
            cvec[ temp[ 1 ] ] <- as.numeric( temp[ 2*j + 1 ] )
         } else {
            if( temp[ 2*j ] %in% names(bvec) ) {
               Amat[ temp[ 2*j ], temp[ 1 ] ] <- as.numeric( temp[ 2*j + 1 ] )
            } else {
               stop( paste( "Constraint name '",temp[ 2*j ],"' is not defined", sep="") )
            }
         }
      }
      i <- i + 1
   }

   ## Restriction values
   if( substr( mps[i], 1, 3 ) != "RHS" ) stop( "MPS file must have a line starting with 'RHS'" )
   i <- i + 1
   while( substr( mps[i], 1, 6 ) != "BOUNDS" & substr( mps[i], 1, 6 ) != "ENDATA" & i < length(mps)) {
      temp <- strsplit( mps[i], " " )[[1]]
      temp <- temp[ temp != "" ]
      for( j in 1:((length(temp)-1)/2) ) {
         if( temp[ 2*j ] %in% names(bvec) ) {
            bvec[ temp[ 2*j ] ] <- as.numeric( temp[ 2*j + 1 ] )
         } else {
            stop( paste( "Constraint name '",temp[ 2*j ],"' is not defined", sep="") )
         }
      }
      i <- i + 1
   }

   ## Bounds
   if( substr( mps[i], 1, 6 ) == "BOUNDS" ) {
      i <- i + 1
      while( substr( mps[i], 1, 6 ) != "ENDATA" & i <= length(mps)) {
         temp <- strsplit( mps[i], " " )[[1]]
         temp <- temp[ temp != "" ]
         if( temp[ 3 ] %in% colnames(Amat) ) {
            if(temp[1] == "UP") {
               svec <- c( svec, "L" )
               bvec <- c( bvec, as.numeric(temp[ 4 ]) )
               Amat <- rbind( Amat, rep( 0, ncol(Amat) ) )
               Amat[ nrow(Amat), temp[3] ] <- 1
               names( svec )[length(svec)] <- paste(temp[1], temp[3], sep="" )
               names( bvec )[length(bvec)] <- paste(temp[1], temp[3], sep="" )
               rownames( Amat )[nrow(Amat)] <- paste(temp[1], temp[3], sep="" )
            } else {
               if( temp[1] %in% c("LO","FX","FR") ) {
                  stop("'LO', 'FX', and 'FR' Bounds are not implemented yet")
               } else {
                  stop(" A 'BOUND' line must start with 'UP', 'LO', 'FX' or 'FR'")
               }
            }
         } else {
            stop( paste( "Variable name '",temp[ 3 ],"' is not defined", sep="") )
         }
         i <- i + 1
      }
   }

   if( substr( mps[i], 1, 6 ) != "ENDATA" )
      stop( "MPS file must have a line starting with 'ENDDATA'" )

   ## Changing 'Greater' constraints to 'Lower' constraints
   for( j in 1:length(svec) ) {
      if(svec[j] == "G" ) {
         bvec[ j ]  <- - bvec[ j ]
         Amat[ j, ] <- - Amat[ j, ]
      }
   }

   res <- NULL
   if(solve) {
      res <- solveLP(cvec,bvec,Amat,maximum)
   }

   result <- list( name=name, cvec=cvec, bvec=bvec, Amat=Amat, res=res )
   return( result )
}

.aidsRestr <- function( nGoods, nShifter = 0, hom = TRUE, sym = TRUE,
      LA = TRUE, restrict.regMat = FALSE ) {

   if( sym && !hom ) {
      hom <- TRUE  # symmetry implies homogeneity
      warning( "symmetry implies homogeneity: imposing additionally homogeniety" )
   }
   nExogEq <- nGoods + 2 + nShifter # number of exog. variables per equation
   if( restrict.regMat ) {
      nExog <- ( nGoods - 1 ) * ( nExogEq )
      restr <- diag( nExog )
      delCols <- NULL
      if( hom ) {
         for( i in 1:( nGoods - 1 ) ) {
            delCol <- ( i - 1 ) * ( nExogEq ) + 2 + nGoods
            addCols <- ( ( i - 1 ) * ( nExogEq ) + 3 ):(
               ( i - 1 ) * ( nExogEq ) + 2 + ( nGoods - 1 ) )
            restr[ delCol, addCols ] <- -1
            delCols <- c( delCols, delCol )
         }
      }
      if( sym ) {
         for( i in 1:( nGoods - 2 ) ) {
            for( j in ( i + 1 ):( nGoods - 1 ) ) {
               delCol <- ( j - 1 ) * ( nExogEq ) + 2 + i
               addCol <- ( i - 1 ) * ( nExogEq ) + 2 + j
               restr[ , addCol ] <- restr[ , addCol ] + restr[ , delCol ]
               delCols <- c( delCols, delCol )
            }
         }
      }
      if( hom || sym ) {
         restr <- restr[ , -delCols ]
      } else {
         restr <- NULL
      }
      if( !is.null( restr ) ) {
         rownames( restr ) <- .aidsCoefNamesEst( nGoods = nGoods,
            nShifter = nShifter, hom = FALSE, sym = FALSE )
         colnames( restr ) <- .aidsCoefNamesEst( nGoods = nGoods,
            nShifter = nShifter, hom = hom, sym = sym )
      }
   } else {
      restr <- NULL
      if( LA ) {
         if( hom ) {
            restr <- matrix( 0, nGoods - 1, ( nGoods - 1 ) * ( nExogEq ) )
            for( i in 1:( nGoods - 1 ) ) {
               for( j in 1:nGoods ) {
                  restr[ i, ( i - 1 ) * ( nExogEq ) + 2 + j ] <- 1
               }
            }
         }
         if( sym ) {
            restr <- rbind( restr, matrix( 0, ( nGoods - 1 ) * ( nGoods - 2 ) / 2,
               ( nGoods - 1 ) * ( nExogEq ) ) )
            k <- 0
            for( i in 1:( nGoods - 2 ) ) {
               for( j in ( i + 1 ):( nGoods - 1 ) ) {
                  k <- k + 1
                  restr[ nGoods - 1 + k, ( i - 1 ) * ( nExogEq ) +
                     2 + j ] <-  1
                  restr[ nGoods - 1 + k, ( j - 1 ) * ( nExogEq ) +
                     2 + i ] <- -1
               }
            }
         }
      }
      if( !is.null( restr ) ) {
         colnames( restr ) <- .aidsCoefNamesEst( nGoods = nGoods,
            nShifter = nShifter, hom = FALSE, sym = FALSE )
      }
   }
   return( restr )
}

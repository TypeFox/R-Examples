solveLP <- function( cvec, bvec, Amat, maximum=FALSE,
               const.dir = rep( "<=", length( bvec ) ),
               maxiter=1000, zero=1e-9, tol=1e-6, dualtol = tol,
               lpSolve=FALSE, solve.dual=FALSE, verbose = 0 )
{

   result <- list()  # list for results that will be returned
   result$status <- 0

   rdigits <- -round( log10( zero ) )

   nVar <- length(cvec)  # number of variables
   nCon <- length(bvec)  # number of constraints

   if( !all.equal( dim( Amat ), c( nCon, nVar ) ) == TRUE ) {
      stop( paste( "Matrix A must have as many rows as constraints (=elements",
         "of vector b) and as many columns as variables (=elements of vector c).\n" ) )
   }
   if( length( const.dir ) != nCon ) {
      stop( paste( "'const.dir' must have as the elements as constraints",
         "(=elements of vector b).\n" ) )
   }
   if( sum( const.dir == ">=" | const.dir == ">" | const.dir == "=" |
         const.dir == "==" | const.dir == "<=" | const.dir == "<" ) < nCon ) {
       stop( "'const.dir' may only contain '>=', '>', '=', '==', '<=' or '<'" )
   }
   if( any( const.dir %in% c( "=", "==" ) ) && ( ! lpSolve ) ) {
      warning( "solveLP() might return incorrect results",
         " if the model includes equality constraints",
         " and argument 'lpSolve' is 'FALSE';",
         " please check if solveLP() returns the same results",
         " with argument 'lpSolve' equal to 'TRUE';",
         " more information on this bug available at",
         " linprog's R-Forge site" )
   }

   ## Labels
   if( is.null(names(cvec))) {
      clab <- as.character(1:nVar)
   } else {
      clab <- names(cvec)
   }
   if( is.null(names(bvec))) {
      blab <- as.character(1:nCon)
   } else {
      blab <- names(bvec)
   }
   const.dir2 <- rep( 0, nCon )
   const.dir2[ const.dir == ">=" | const.dir == ">" ] <-  1
   const.dir2[ const.dir == "<=" | const.dir == "<" ] <- -1

   ## lpSolve
   if( lpSolve ) {
      library( lpSolve )
      if( maximum ) {
         direction <- "max"
      } else {
         direction <- "min"
      }
      lpres <- lp( direction = direction, cvec, Amat, const.dir, bvec )
      result$lpStatus <- lpres$status
      if( result$lpStatus == 0 ) {
         if( min( lpres$solution ) < -tol ) {
            result$lpStatus <- 7
         } else if( max( ( bvec - c( Amat %*% lpres$solution ) ) *
               const.dir2 ) > tol ) {
            result$lpStatus <- 3
         }
      }

      if( result$lpStatus != 0 ) {
         result$status <- 1
      } else  {
         result$solution  <- lpres$solution
         names( result$solution ) <- clab
         result$opt    <- lpres$objval

         ## Results: Constraints
         result$con <- data.frame( actual=NA, dir=const.dir, bvec=bvec, free=NA )
         result$con$actual <- round( c( Amat %*% result$solution ), digits=rdigits )
         names( result$con$actual ) <- blab
         result$con$free   <- round( result$con$bvec - result$con$actual, digits=rdigits )
         result$con$free[ const.dir2 == 1 ] <- -result$con$free[ const.dir2 == 1 ]
         result$con$free[ const.dir2 == 0 ] <- -abs( result$con$free[ const.dir2 == 0 ] )
      }

   } else {
      ## Simplex algorithm
      iter1 <- 0
      iter2 <- 0

      ## Slack Variables
      for(i in 1:nCon) clab <- c( clab, paste("S", blab[i] ) )
      cvec2 <- c( cvec, rep( 0, nCon ) )

      ## Tableau ( Basic Variables, Slack,Variables, P0, Z-C )
      Tab <- rbind( cbind( -Amat * const.dir2, diag( 1, nCon, nCon ), -bvec * const.dir2 ),
                 c( cvec2 * (-1)^maximum, 0 ) )
      rownames(Tab) <- c( blab, "Z-C" )
      colnames(Tab) <- c( clab, "P0" )
      if( verbose >= 3 ) {
         print("initial Tableau")
         print(Tab)
      }

      ## searching for feasible solution for starting
      # basis: Zero Solution ( Basic Variables = Slack Variables )
      basis <- c( (nVar+1) : (nVar+nCon) )
      if(min(Tab[ 1:nCon, nVar+nCon+1]) < 0 ) {
         Tab2 <- Tab
         Tab2 <- rbind( Tab2, rep(0, ncol(Tab2) ) )
         rownames(Tab2)[nCon+2] <- "M Z-C"     # additional artificial 'Z-C' row
         basis2 <- basis
         nArt   <- 0       # number of artificial variables
         for(i in 1:nCon) {
            if(Tab[ i, nVar+nCon+1] < 0 ) {
               Tab2[ i, ] <- -Tab2[ i, ]
               Tab2 <- cbind( Tab2[ , 1:(nVar+nCon+nArt) ], rep(0,nCon+2),
                              Tab2[ , (nVar+nCon+nArt+1)] )
               nArt <- nArt + 1
               colnames(Tab2)[ nVar+nCon+nArt ] <- paste("M", rownames(Tab2)[i] )
               Tab2[ i, nVar+nCon+nArt ] <- 1
               Tab2[ nCon+2, nVar+nCon+nArt ] <- 1
               # put artificial variables in basis
               rownames(Tab2)[ i ] <- paste("M", rownames(Tab2)[i] )
               basis2[i] <- nVar+nCon+nArt
            }
         }
         for(i in 1:nCon) {    # artificial objective function (Z-C)
            if(Tab[ i, nVar+nCon+1] < 0 ) {
               Tab2[nCon+2, 1:(nVar+nCon+nArt)] <- Tab2[nCon+2, 1:(nVar+nCon+nArt)] -
                                                   Tab2[ i , 1:(nVar+nCon+nArt)]
            }
         }
         for(i in 1:nCon) {    # value of artificial objective function
            Tab2[nCon+2, nVar+nCon+nArt+1 ] <- Tab2[nCon+2, nVar+nCon+nArt+1 ] -
                         Tab2[i, nVar+nCon+nArt+1] * Tab2[nCon+2, basis[i] ]
         }
         colnames(Tab2)[ nVar+nCon+nArt+1 ] <- "P0"
         if( verbose >= 3 ) {
            print("initial Tableau for Phase 1")
            print(Tab2)
         }

         ## Simplex algorithm (Phase 1)
         while( min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ] ) < -zero & iter1 < maxiter) {
            iter1 <- iter1 + 1
           ## Pivot
            Tab[ abs(Tab) < zero ] <- 0
#            pcolumn <- which.min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ]) # Pivot column
            decval <- array( NA, nVar+nCon )
            for( pcolumnt in 1:(nVar+nCon+nArt) ) {
               if( Tab2[ nCon+2, pcolumnt ] < 0 ) {
                  rwerte  <- Tab2[ 1:nCon, nVar+nCon+nArt+1 ] / Tab2[ 1:nCon , pcolumnt ]
                       # R-values
                  rwerte[ Tab2[1:nCon, pcolumnt ] <= 0 ] <- max(rwerte,na.rm=TRUE)+1
                  prow  <- which.min( rwerte )    # Pivot row
                  if( length( rwerte[ !is.na(rwerte) & is.finite(rwerte) ] ) >= 1 ) {
                     decval[ pcolumnt ] <- Tab2[ nCon+2, pcolumnt ] *
                           min( rwerte[ !is.na(rwerte) & is.finite(rwerte) ] )
                  }
               }
            }
            if( min( decval, na.rm=TRUE ) < -zero ) {
               pcolumn <- which.min( decval ) # Pivot column
            } else {
               pcolumn <- which.min( Tab2[ nCon+2, 1:(nVar+nCon+nArt) ]) # Pivot column
            }
            rwerte  <- Tab2[ 1:nCon , nVar+nCon+nArt+1 ] / Tab2[ 1:nCon , pcolumn ] # R-values
            rwerte[ Tab2[1:nCon, pcolumn ] <= 0 ] <- max(rwerte, na.rm=TRUE)+1
            prow  <- which.min( rwerte )    # Pivot row
            if( verbose >=2 ) {
               cat( paste( "\nPivot Column:", as.character(pcolumn),
                           "(",colnames(Tab2)[pcolumn],")\n" ) )
               cat( paste( "Pivot Row:", as.character( prow ),
                  "(", rownames(Tab2)[prow], ")\n\n" ) )
            }

            ## New Basis
            basis[prow] <- pcolumn
            rownames(Tab2)[prow] <- colnames(Tab2)[pcolumn]

            ## New Tableau
            Tab2[ prow, ] <- Tab2[ prow, ] / Tab2[ prow, pcolumn ]
            for( i in 1:(nCon+2) ) {
               if( i != prow ) {
                  Tab2[ i, ] <- Tab2[ i, ] - Tab2[ prow, ] * Tab2[ i, pcolumn ]
               }
            }
            if( verbose >= 4 ) print(Tab2)
         }
         if(iter1 >= maxiter ) warning("Simplex algorithm (phase 1) did not reach optimum.")
         Tab <- cbind( Tab2[ 1:(nCon+1), 1:(nCon+nVar) ],
            Tab2[ 1:(nCon+1), nVar+nCon+nArt+1 ] )
         if( verbose >= 3 ) {
            print("New starting Tableau for Phase II")
            print(Tab)
         }
      }
      ## Simplex algorithm (Phase 2)
      while( min( Tab[ nCon+1, 1:(nVar+nCon) ] ) < -zero  & iter2 < maxiter ) {
         iter2 <- iter2 + 1
         ## Pivot
         Tab[ abs(Tab) < zero ] <- 0
#         pcolumn <- which.min( Tab[ nCon+1, 1:(nVar+nCon) ]) # Pivot column
         decval <- array( NA, nVar+nCon )
         for( pcolumnt in 1:(nVar+nCon) ) {
            if( Tab[ nCon+1, pcolumnt ] < 0 ) {
               rwerte  <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon , pcolumnt ]
                    # R-values
               rwerte[ Tab[1:nCon, pcolumnt ] <= 0 ] <- max(rwerte,na.rm=TRUE)+1
               prow  <- which.min( rwerte )    # Pivot row
               if( length( rwerte[ !is.na(rwerte) & is.finite(rwerte) ] ) >= 1 ) {
                  decval[ pcolumnt ] <- Tab[ nCon+1, pcolumnt ] *
                        min( rwerte[ !is.na(rwerte) & is.finite(rwerte) ] )
               }
            }
         }
         if( min( decval, na.rm=TRUE ) < -zero ) {
            pcolumn <- which.min( decval ) # Pivot column
         } else {
            pcolumn <- which.min( Tab[ nCon+1, 1:(nVar+nCon) ]) # Pivot column
         }
         rwerte  <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon , pcolumn ]     # R-values
         rwerte[ Tab[1:nCon, pcolumn ] <= 0 ] <- max(rwerte,na.rm=TRUE)+1
         prow  <- which.min( rwerte )    # Pivot row
         if( verbose >= 2 ) {
            cat( paste( "\nPivot Column:", as.character(pcolumn),
                        "(",colnames(Tab)[pcolumn],")\n" ) )
            cat( paste( "Pivot Row:", as.character( prow ) ,
               "(",rownames(Tab)[prow],")\n\n") )
         }

         ## New Basis
         basis[prow] <- pcolumn
         rownames(Tab)[prow] <- colnames(Tab)[pcolumn]

         ## New Tableau
         Tab[ prow, ] <- Tab[ prow, ] / Tab[ prow, pcolumn ]
         for( i in 1:(nCon+1) ) {
            if( i != prow ) {
               Tab[ i, ] <- Tab[ i, ] - Tab[ prow, ] * Tab[ i, pcolumn ]
            }
         }
         if( verbose >= 4 ) print(Tab)
      }
      if(iter2 >= maxiter ) warning("Simplex algorithm (phase 2) did not reach optimum.")
      ## Results: Basic Variables
      basvar <- matrix( NA, nCon, 1 )
      colnames(basvar) <- c("opt")
      rownames(basvar) <- rep("a",nCon)
      for( i in 1:nCon ) {
         rownames(basvar)[i] <- clab[sort(basis)[i]]
         basvar[i,1] <- Tab[ which(basis==sort(basis)[i]), nVar+nCon+1 ]
      }

      ## Results: All Variables (Including Slack Variables)
      allvar <- data.frame( opt=rep( NA, nVar+nCon ), cvec=cvec2, min.c=NA,
                              max.c=NA, marg=NA, marg.reg=NA )
      rownames(allvar) <- clab
      for( i in 1:(nVar+nCon) ) {
         if(i %in% basis ) {
            allvar$opt[ i ] <- Tab[ which(basis==i), nVar+nCon+1 ]
            ## Stability of Basic Variables
            quot <- Tab[ nCon+1, 1:(nVar+nCon) ] / Tab[ which(basis==i), 1:(nVar+nCon) ]
            if( maximum ) {
               if(max(quot[!is.na(quot)]) > 0 ) {
                  suppressWarnings(
                     allvar$min.c[ i ] <- cvec2[ i ] - min(quot[quot>0 & !is.na(quot)])
                  )
               }
               if(min(quot[!is.na(quot) & is.finite(quot)]) < 0 ) {
                  if(max(quot[quot<0 & !is.na(quot)]) > -1e14 ) {
                     allvar$max.c[ i ] <- cvec2[ i ] - max(quot[quot<0 & !is.na(quot)])
                  } else {
                     allvar$max.c[ i ] <- Inf
                  }
               } else {
                  allvar$max.c[ i ] <- Inf
               }
            } else {
               if(max(quot[!is.na(quot)]) > 0 ) {
                  suppressWarnings(
                     allvar$max.c[ i ] <- cvec2[ i ] + min(quot[quot>0 & !is.na(quot)])
                  )
               }
               if(min(quot[!is.na(quot)]) < 0 ) {
                  if(max(quot[quot<0 & !is.na(quot)]) > -1e14 ){
                     allvar$min.c[ i ] <- -cvec2[ i ] + max(quot[quot<0 & !is.na(quot)])
                  } else {
                     allvar$min.c[ i ] <- NA
                  }
               } else {
                  allvar$min.c[ i ] <- NA
               }
            }
         } else {
             allvar$opt[ i ] <- 0
             if( i <= nVar ) {
                if( maximum ) {
                   allvar$min.c[ i ] <- -Inf
                   allvar$max.c[ i ] <- Tab[ nCon+1, i ] + cvec2[i]
                } else {
                   allvar$min.c[ i ] <- 99#-Tab[ nCon+1, i ] - cvec2[i]
                   allvar$max.c[ i ] <- 77#Inf
                }
             }
         }
         allvar$cvec[ i ] <- cvec2[ i ]

         # marginal contribution to objective function (Shadow prices)
         if( !( ( i %in% basis ) & ( i <= nVar ) ) ) {
            allvar$marg[ i ] <- Tab[ nCon+1, i ] * (-1)^maximum
            if( !( i %in% basis ) & ( i > nVar ) ) {
               if( maximum ) {
                  allvar$max.c[ i ] <- Tab[ nCon+1, i ] #* (-1)^maximum
                  allvar$min.c[ i ] <- -Inf
               } else {
                  allvar$min.c[ i ] <- -Tab[ nCon+1, i ]
                  allvar$max.c[ i ] <-  Inf
               }
            }
            quot <- Tab[ 1:nCon , nVar+nCon+1 ] / Tab[ 1:nCon, i ]
               suppressWarnings(
               if( !( i %in% basis) ) {
                  allvar$marg.reg[ i ] <- min(quot[quot>0 & !is.na(quot)])
               } else {
                  allvar$marg.reg[ i ] <- NA
               }
            )
         }
      }
      allvar$min.c[ allvar$min.c >  1e16 ] <-  Inf
      allvar$min.c[ allvar$min.c < -1e16 ] <- -Inf

      ## Results: Constraints
      con <- data.frame( actual=NA, dir=const.dir, bvec=bvec, free=NA, dual=NA, dual.reg=NA )
      names( con$actual ) <- blab
      for(i in 1: nCon) {
         if( (i+nVar) %in% basis ) {
            con$actual[ i ] <- round( bvec[i] + Tab[ which((i+nVar)==basis),
                                 nVar+nCon+1 ] * const.dir2[ i ], digits=rdigits )
         } else {
            con$actual[ i ] <- round( bvec[i], digits=rdigits )
         }
         if( -allvar$opt[ i+nVar ] == 0 ) {
            con$dual[ i ]     <- allvar$marg[ i+nVar ] * (-1)^maximum
            con$dual.reg[ i ] <- allvar$marg.reg[ i+nVar ]
         } else {
            con$dual[ i ]     <- 0
            con$dual.reg[ i ] <- allvar$opt[ i+nVar ]
         }
      }
      con$free <- round( con$bvec - con$actual, digits=rdigits )
      con$free[ const.dir2 == 1 ] <- -con$free[ const.dir2 == 1 ]
      con$free[ const.dir2 == 0 ] <- -abs( con$free[ const.dir2 == 0 ] )

      result$opt      <- round( -Tab[ nCon+1, nCon+nVar+1 ], digits=rdigits ) * (-1)^maximum
      result$iter1    <- iter1
      result$iter2    <- iter2
      result$allvar   <- round( allvar, digits=rdigits )
      result$basvar   <- round( basvar, digits=rdigits )
      result$solution <- result$allvar$opt[ 1 : nVar ]
      names( result$solution ) <- clab[ 1: nVar ]

      result$con      <- con
      if( verbose >= 1 ) result$Tab <- Tab
      if( iter1 >= maxiter ) result$status <- 4
      if( iter2 >= maxiter ) result$status <- 5
   }

   if( result$status == 0 ) {
      if( min ( result$con$free ) < - tol ) {
         result$status <- 3
      }
   }

   ## solving the dual problem
   if( solve.dual && result$status == 0 ) {
      if( any( const.dir2 == 0 ) ) {
         stop( paste( "At the moment the dual problem can not be solved",
            "with equality constraints" ) )
      }

      if( maximum ) {
         const.dir.dual <- rep(">=",nVar)
      } else {
         const.dir.dual <- rep("<=",nVar)
      }

      result$con$dual.p <- result$con$dual
      dualres <- solveLP( cvec = bvec * const.dir2 * (-1)^maximum, bvec = cvec,
         Amat = t( Amat * const.dir2 ) * (-1)^maximum, maximum = !maximum,
         const.dir = const.dir.dual, maxiter = maxiter, zero = zero,
         tol = dualtol, lpSolve = lpSolve, verbose = verbose )
      result$dualStatus <- dualres$status
      if( result$dualStatus == 0 ) {
         result$con$dual <- dualres$solution
      } else {
         result$status <- 2
      }
   }

   ## List of Results
   result$maximum  <- maximum
   result$lpSolve  <- lpSolve
   result$solve.dual <- solve.dual
   result$maxiter    <- maxiter
   class(result)   <- "solveLP"
   result
}

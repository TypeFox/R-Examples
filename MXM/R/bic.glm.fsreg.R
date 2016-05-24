bic.glm.fsreg <- function( target, dataset, robust = FALSE, tol = 0, ncores = 1, maxit = 100 ) {

  p <- ncol(dataset)  ## number of variables
  bico <- numeric( p )
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  tool <- NULL
  oiko <- NULL
  info <- matrix( 0, ncol = 2 )
  # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )


 #########
  ## if it is binomial or poisson regression
  #########

  if ( length( unique(target) ) == 2 || sum( round(target) - target ) == 0 ) {

   if ( length( unique(target) ) == 2 ) {
      oiko <- "binomial"  ## binomial regression
    } else {
      oiko <- "poisson"  ## poisson regression
    }
    
    runtime <- proc.time()

    if ( robust == FALSE ) {
      ini = glm( target ~ 1, family = oiko )$deviance + log(n)  ## initial BIC
    } else {
      ini = robust::glmRob( target ~ 1, family = oiko, maxit = maxit )$deviance + log(n)  ## initial BIC
    }

    if (ncores <= 1) {
      if ( robust == FALSE ) {  ## Non robust
        for (i in 1:p) {
          mi <- glm( target ~ dataset[, i], family = oiko )
          bico[i] <- BIC( mi )
        }

      } else {  ## Robust
        for (i in 1:p) {
          mi <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit )
          bico[i] <- mi$deviance + length( coef( mi ) ) * log(n)
        }
      }
      mat <- cbind(1:p, bico)

      if( any( is.na(mat) ) ) {
        mat[ which( is.na(mat) ) ] = ini
      }

    } else {
      if ( robust == FALSE ) {  ## Non robust
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        bico <- numeric(p)
        mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
          ww <- glm( target ~ dataset[, i], family = oiko )
          bico[i] <- BIC( ww )
        }
        stopCluster(cl)

      } else {  ## Robust
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        bico <- numeric(p)
        mod <- foreach( i = 1:p, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
          ww <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit )
          bico[i] <- ww$deviance + length( coef( ww ) ) * log(n)
        }
        stopCluster(cl)
      }

      mat <- cbind(1:p, mod)

      if ( any( is.na(mat) ) ) {
        mat[ which( is.na(mat) ) ] = ini
      }
      
    }

    colnames(mat) <- c("variable", "BIC")
    rownames(mat) <- 1:p
    sel <- which.min( mat[, 2] )
    sela <- sel

    if ( mat[sel, 2] < ini ) {

      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ]
      if ( !is.matrix(mat) ) {
        mat <- matrix(mat, ncol = 2) 
      }
      mat <- mat[ order( mat[, 2] ), ]

      if ( robust == FALSE ) {
        mi = glm( target ~ dataset[, sel], family = oiko )
        tool[1] <- BIC( mi )
      } else {
        mi = robust::glmRob( target ~ dataset[, sel], family = oiko, maxit = maxit )
        tool[1] <- mi$deviance + length( coef( mi ) ) * log(n)
        if ( is.na(tool[1]) )  tool[1] <- ini
      }

      moda[[ 1 ]] <- mi
    }


    ######
    ###     k equals 2
    ######


    if ( length(moda) > 0  &  nrow(mat) > 0 ) {

      k <- 2
      pn <- p - k  + 1
      mod <- list()

      if ( ncores <= 1 ) {
        bico <- numeric( pn )
        if ( robust == FALSE ) {  ## Non robust
          for ( i in 1:pn ) {
            ma <- glm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), family = oiko )
            bico[i] <- BIC( ma )
          }

        } else {  ## Robust
          for ( i in 1:pn ) {
            ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), family = oiko, maxit = maxit )
            bico[i] <- ma$deviance + length( coef( ma ) ) * log(n)
          }
        }

        mat[, 2] <- bico

      } else {

        if ( robust == FALSE ) {  ## Non robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- glm( target ~ dataset[, sel ] + dataset[, mat[i, 1] ], family = oiko )
            bico[i] <- BIC( ww )
          }

          stopCluster(cl)

        } else {  ## Robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
            ww <- glmRob( target ~ dataset[, sel ] + dataset[, mat[i, 1] ], family = oiko, maxit = maxit )
            bico[i] <- ww$deviance + length( coef( ww ) ) * log(n)
          }

          stopCluster(cl)
        }

        mat[, 2] <- mod

      }

      ina <- which.min( mat[, 2] )
      sel <- mat[ina, 1]

      if ( mat[ina, 2] >= tool[1] ) {
        info <- rbind( info,  c( -10, 1e300 ) )

      } else {
        tool[2] <- mat[ina, 2]
        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina , ]
        if ( !is.matrix(mat) ) {
          mat <- matrix(mat, ncol = 2) 
        }
        mat <- mat[ order( mat[, 2] ), ]

        if ( robust == FALSE ) {
          mi = glm( target ~., data =as.data.frame( dataset[, sela] ), family = oiko )
          tool[2] <- BIC( mi )
        } else {
          mi = robust::glmRob( target ~ ~., data =as.data.frame( dataset[, sela] ), family = oiko, maxit = maxit )
          tool[2] <- mi$deviance + length( coef( mi ) ) * log(n)
        }
        moda[[ 2 ]] <- mi

      }

     if ( sum( info[2, ] -  c( -10, 1e300 ) ) == 0 ) {
       info <- matrix( info[1, ], ncol = 3)
     }

   }


      #########
      ####      k is greater than 2
      #########


    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( ( k < n - 10 ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) ) {

        k <- k + 1
        pn <- p - k + 1

        if (ncores <= 1) {
          if ( robust == FALSE ) {  ## Non robust
            for ( i in 1:pn ) {
              ma <- ma <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko )
              mat[i, 2] <- BIC( ma )
            }

          } else {  ## Robust
            for ( i in 1:pn ) {
              ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
              mat[i, 2] <- ma$deviance + length( coef( ma ) ) * log(n)
            }
          }

        } else {
          if ( robust == FALSE ) {  ## Non robust
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            bico <- numeric(pn)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko )
              bico[i] <- BIC( ww )
            }

           stopCluster(cl)

          } else {  ## Robust
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            bico <- numeric(pn)
            mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
              ww <- glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
              bico[i] = ww$deviance + length( coef( ww ) ) * log(n)
            }

            stopCluster(cl)
          }
          mat[, 2] <- mod

        }

        ina <- which.min( mat[, 2] )
        sel <- mat[ina, 1]

        if ( mat[ina, 2] >= tool[k - 1] ) {
          info <- rbind( info,  c( -10, 1e300 ) )
          tool[k] <- 1e+300

        } else {

          tool[k] = mat[ina, 2]
          info <- rbind(info, mat[ina, ] )
          sela <- info[, 1]
          mat <- mat[-ina , ]
          if ( !is.matrix(mat) ) {
            mat <- matrix(mat, ncol = 2) 
          }
          mat <- mat[ order( mat[, 2] ), ]

          if ( robust == FALSE ) {
            ma = glm( target ~., data =as.data.frame( dataset[, sela] ), family = oiko )
            tool[k] <- BIC( mi )
          } else {
            ma = robust::glmRob( target ~., data =as.data.frame( dataset[, sela] ), family = oiko, maxit = maxit )
            tool[k] <- ma$deviance + length( coef( ma ) ) * log(n)
          }
          moda[[ k ]] <- ma

        }

      }

    }

    runtime <- proc.time() - runtime


    ##############
    ##############
    ##          ##
    ##          ##  else it is linear regression
    ##          ##
    ##############
    ##############


  } else {

  
    runtime <- proc.time()

      if ( robust == FALSE ) {  ## Non robust
        mi = lm( target ~ 1 )
        ini <- BIC( mi )

      } else {  ## Robust
         mi = MASS::rlm( target ~ dataset[, i], maxit = 2000 )
         ini = BIC(mi)
      }

    if (ncores <= 1) {
      if ( robust == FALSE ) {  ## Non robust
       for (i in 1:p) {
         mi = lm( target ~ dataset[, i] )
         bico[i] <- BIC( mi )
       }

      } else {  ## Robust
        for (i in 1:p) {
         mi = MASS::rlm( target ~ dataset[, i], maxit = 2000 )
         bico[i] = BIC(mi)
        }
      }

     mat <- cbind(1:p, bico)

    } else {
       if ( robust == FALSE ) {  ## Non robust
         cl <- makePSOCKcluster(ncores)
         registerDoParallel(cl)
         bico <- numeric(p)
         mod <- foreach( i = 1:p, .combine = rbind ) %dopar% {
           ww <- lm( target ~ dataset[, i] )
           bico[i] <- BIC( ww )
         }
         stopCluster(cl)

       } else {  ## Robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(p)
          mod <- foreach( i = 1:p, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
            ww <- rlm( target ~ dataset[, i], maxit = 2000 )
            bico[i] = BIC(ww)
          }
          stopCluster(cl)
       }
      mat <- cbind(1:p, mod)

    }

    colnames(mat) <- c("variable", "BIC")
    rownames(mat) <- 1:p
    sel <- which.min( mat[, 2] )
    sela <- sel

    if ( mat[sel, 2] < ini ) {
      if ( robust == FALSE ) {
        mi = lm( target ~ dataset[, sel])
        tool[1] <- BIC( mi )
      } else {
        mi = MASS::rlm( target ~ dataset[, sel], maxit = 2000 )
         tool[1] <- BIC(mi)
      }
      moda[[ 1 ]] <- mi
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ]
      if ( !is.matrix(mat) ) {
        mat <- matrix(mat, ncol = 2) 
      }
      mat <- mat[ order( mat[, 2] ), ]    }

    ######
    ###     k equals 2
    ######

    if ( length(moda) > 0 ) {

      k <- 2
      pn <- p - k + 1
      mod <- list()

      if ( ncores <= 1 ) {
         if ( robust == FALSE ) {
          for (i in 1:pn) {
            ma = lm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ) )
            mat[i, 2] <- BIC( ma )
          }
         } else {
           for (i in 1:pn) {
             ma = MASS::rlm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), maxit = 2000 )
             mat[i, 2] <- BIC( ma )
           }
         }

      } else {

       if ( robust == FALSE ) {  ## Non robust
         cl <- makePSOCKcluster(ncores)
         registerDoParallel(cl)
         bico <- numeric(pn)
         mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
           ww <- lm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ) )
           bico[i] <- BIC( ww )
         }
         stopCluster(cl)

       } else {  ## Robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
            ww <- rlm( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ), maxit = 2000 )
            bico[i] <- BIC( ww )
          }
          stopCluster(cl)
       }

      mat[, 2] <- mod
    }

      ina <- which.min( mat[, 2] )
      sel <- mat[ina, 1]

      if ( mat[ina, 2] >= tool[1] ) {
        info <- rbind( info,  c( -10, 1e300 ) )

      } else {
         if ( robust == FALSE ) {
          ma <- lm( target ~ dataset[, info[1, 1] ] + dataset[, sel ] )
          tool[2] <- BIC( ma )
         } else {
           ma <- MASS::rlm( target ~ dataset[, info[1, 1] ] + dataset[, sel ], maxit = 2000 )
           tool[2] <- BIC( ma )
         }

        info <- rbind(info, mat[ina, ] )
        sela <- info[, 1]
        mat <- mat[-ina , ]
        if ( !is.matrix(mat) ) {
          mat <- matrix(mat, ncol = 2) 
        }
        mat <- mat[ order( mat[, 2] ), ]

        moda[[ 2 ]] <- ma
      }

      if ( sum( info[2, ] -  c( -10, 1e300 ) ) == 0 ) {
        info <- matrix( info[1, ], ncol = 2)
      }

    }

      #########
      ####      k is greater than 2
      #########

    if ( nrow(info) > 1 ) {
      while ( ( k < n - 10 ) & ( tool[ k - 1 ] - tool[ k ] > tol ) ) {

        k <- k + 1
        pn <- p - k + 1

        if ( ncores <= 1 ) {
           if ( robust == FALSE ) {
            for ( i in 1:pn ) {
              ma <- lm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ) )
              mat[i, 2] = BIC( ma )
            }

           } else {
             for ( i in 1:pn ) {
               ma <- MASS::rlm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), maxit = 2000 )
               mat[i, 2] = BIC( ma )
             }
           }

        } else {
           if ( robust == FALSE ) {  ## Non robust
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            bico <- numeric(pn)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- lm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ) )
              bico[i] = BIC( ww )
            }
          stopCluster(cl)

           } else {  ## Robust
             cl <- makePSOCKcluster(ncores)
             registerDoParallel(cl)
             bico <- numeric(pn)
             mod <- foreach( i = 1:p, .combine = rbind, export = "rlm", .packages = "MASS" ) %dopar% {
               ww <- rlm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), maxit = 2000 )
               bico[i] = BIC( ww )
             }
             stopCluster(cl)
           }
          
          mat[, 2] <- mod

        }

         ina <- which.min( mat[, 2] )
         sel <- mat[ina, 1]

         if ( mat[ina, 2] >= tool[k - 1] ) {
           info <- rbind( info,  c( -10, 1e300 ) )
           tool[k] <- 1e+300

         } else {
            if ( robust == FALSE ) {
             ma <- lm( target ~., data = as.data.frame( dataset[, c(sela, sel) ] )  )
             tool[k] <- BIC( ma )
            } else {
              ma <- MASS::rlm( target ~., data = as.data.frame( dataset[, c(sela, sel) ] ), maxit = 2000 )
              tool[k] <- BIC( ma )
            }

           info <- rbind(info, mat[ina, ] )
           sela <- info[, 1]
           mat <- mat[-ina , ]
           if ( !is.matrix(mat) ) {
             mat <- matrix(mat, ncol = 2) 
           }
           mat <- mat[ order( mat[, 2] ), ]

           moda[[ k ]] <- ma
        }

      }
    }

    runtime <- proc.time() - runtime
  }

    d <- length(sela)
    final <- NULL
    models <- NULL

    if ( d >= 1 ) {
      models <- NULL
      xx <- as.data.frame( dataset[, sela] )
      colnames(xx) <- paste("V", sela, sep = "")
      if ( d == 1 ) {

       models <- NULL
        xx <- as.data.frame( dataset[, sela] )
        colnames(xx) <- paste("V", sela, sep = "")

        if ( is.null(oiko) ) {
           if ( robust == FALSE ) {
            models <- final <- lm( target ~., data = as.data.frame( xx ) )
           } else {
             models <- final <- MASS::rlm( target ~., data = as.data.frame( xx ), maxit = 2000 )
           }

        } else {
          if ( robust == FALSE ) {
            models <- final <- glm( target ~., data = as.data.frame( xx ), family = oiko )
          } else {
            models <- final <- robust::glmRob( target ~., data = as.data.frame( xx ), family = oiko, maxit = maxit )
          }
        }

      } else {
        if ( is.null(oiko) ) {
          for (i in 1:d) {
             if ( robust == FALSE ) {
              models[[ i ]] <- lm( target ~., data = as.data.frame( xx[, 1:i] ) )
             } else {
               models[[ i ]] <- MASS::rlm( target ~., data = as.data.frame( xx[, 1:i] ), maxit = 2000 )
             }
          }

        } else {
          for (i in 1:d) {
            if ( robust == FALSE ) {
              models[[ i ]] <- glm( target ~., data = as.data.frame( xx[, 1:i] ), family = oiko )
            } else {
              models[[ i ]] <- robust::glmRob( target ~., data = as.data.frame( xx[, 1:i] ), family = oiko, maxit = maxit )
            }
          }
        }

        final <- summary( models[[ d ]] )
      }

      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      colnames(info) <- c( "variables", "BIC" )
      rownames(info) <- info[, 1]
    }

    list(mat = t(mat), info = info, models = models, final = final, runtime = runtime )

}

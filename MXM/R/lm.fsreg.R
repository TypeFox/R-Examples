lm.fsreg <- function(target, dataset, threshold = 0.05, stopping = "BIC", tol = 2, robust = FALSE, ncores = 1 ) {
  
  threshold = log(threshold)

  p <- ncol(dataset)  ## number of variables
  pval <- stat <- dof <- numeric( p )  
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  tool <- numeric( length( min(n, p) ) )
  
  # cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )

  runtime <- proc.time()
    
    if (ncores <= 1) {
     if ( robust == FALSE ) {  ## Non robust
       for (i in 1:p) {
         ww = lm( target ~ dataset[, i] )
         tab = anova(ww)
         stat[i] = tab[1, 4] 
         df1 = tab[1, 1]   ;  df2 = tab[2, 1]
         pval[i] = pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
       }
      mat <- cbind(1:p, pval, stat)
       
     } else {  ## Robust
       for (i in 1:p) {
        # ww = robust::lmRob( target ~ dataset[, i], control = cont )
        # tab = aovlmrob( ww )
        ww = MASS::rlm(target ~ dataset[, i ], maxit = 2000) 
        stat[i] = 2 * as.numeric( logLik(ww) )
        dof[i] = length( coef(ww) )
       }
       fit0 = MASS::rlm(target ~ 1, maxit = 2000) 
       stat0 = 2 * as.numeric( logLik(fit0) )
       difa = abs( stat - stat0 )
       pval = pchisq(difa, dof - 1, lower.tail = FALSE, log.p = TRUE)
       mat <- cbind(1:p, pval, difa)
     }

    } else {
      if ( robust == FALSE ) {  ## Non robust
         cl <- makePSOCKcluster(ncores)
         registerDoParallel(cl)
         mat <- matrix(0, p, 2)
         mod <- foreach( i = 1:p, .combine = rbind ) %dopar% {
           ww <- lm( target ~ dataset[, i] )
           tab <- anova( ww )
           stat[i] <- tab[1, 4] 
           df1 <- tab[1, 1]   ;  df2 = tab[2, 1]
           pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
           mat[i, ] <- c(pval[i], stat[i]) 
         }
         stopCluster(cl)

      } else {  ## Robust
        fit0 = MASS::rlm(target ~ 1, maxit = 2000) 
        stat0 = 2 * logLik(fit0)
        
         cl <- makePSOCKcluster(ncores)
         registerDoParallel(cl)
         mat <- matrix(0, p, 2)
         mod <- foreach( i = 1:p, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
           # ww <- lmRob( target ~ dataset[, i], control = cont )
           # tab <- aovlmrob( ww )
           # mat[i, ] <- c(tab[2, 3], tab[2, 2] ) 
           ww = rlm( target ~ dataset[, i], maxit = 2000 ) 
           mat[i, ] = c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) )
         }
         stopCluster(cl)
         
         difa = abs( mod[, 1] - stat0 )
         pval = pchisq(difa, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE)
         mod = cbind( pval, difa)
      }

      mat <- cbind(1:p, mod)      
  
    }
    
    colnames(mat) <- c( "variables", "p-value", "stat" )
    rownames(mat) <- 1:p
    sel <- which.min(mat[, 2])
    info <- matrix( numeric(3), ncol = 3 )
    sela <- sel
    
    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ] 
      if ( !is.matrix(mat) ) {
        mat <- matrix(mat, ncol = 3) 
      }
      mat <- mat[ order( mat[, 2] ), ] 

      if ( stopping == "adjrsq" ) {

        if ( robust == FALSE ) {
          mi = lm( target ~ dataset[, sel] )
          tool[1] <- as.numeric( summary( mi )[[ 9 ]] )
        } else {
          # mi = robust::lmRob( target ~ dataset[, sel], control = cont )
          # r2 = as.numeric( summary(mi)[[ 8 ]] )
          mi = MASS::rlm(target ~ dataset[, sel], maxit = 2000)
          r2 = cor( target, fitted(mi) )^2
          tool[1] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(mi) ) - 1 )
        }         
      } else if ( stopping  == "BIC" ) { 
        if ( robust == FALSE ) {
          mi = lm( target ~ dataset[, sel] )
          tool[1] <- BIC( mi )
        } else {
          # mi = robust::lmRob( target ~ dataset[, sel], control = cont )
          mi = MASS::rlm(target ~ dataset[, sel], maxit = 2000)
          tool[1] <-  BIC(mi)
        }
      }
      moda[[ 1 ]] <- mi
      
    }
    
 
     ######
     ###### k equal to 2
     ######

    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- 2
      pn <- p - k + 1   
      
      if ( ncores <= 1 ) {
        if ( robust == FALSE ) {      
          for (i in 1:pn) {
            ww = lm( target ~ dataset[, sel] + dataset[, mat[i, 1] ] )
            tab = anova( ww )
            mat[i, 3] = tab[2, 4] 
            df1 = tab[2, 1]   ;  df2 = tab[3, 1]
            mat[i, 2] = pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
          }
          
        } else {
          do = 2 * as.numeric( logLik( moda[[ 1 ]] ) )
          fr = length( coef( moda[[ 1 ]] ) )
          sta = dof = numeric(pn)
          
          for (i in 1:pn) {
            # ww = robust::lmRob( target ~ dataset[, sel] + dataset[, mat[i, 1] ], control = cont )
            # tab = aovlmrob( ww )
            # mat[i, 3] = tab[3, 2] 
            # mat[i, 2] = tab[3, 3]
            ww = MASS::rlm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], maxit = 2000 )
            sta[i] = 2 * as.numeric( logLik(ww) )
            dof[i] = length( coef(ww) )
          } 
          mat[, 3] = abs( sta - do )
          mat[, 2] = pchisq(mat[, 3], dof - fr, lower.tail = FALSE, log.p = TRUE)
        }

      } else {

        if ( robust == FALSE ) {  ## Non robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          stat <- pval <- numeric(pn) 
          mata <- matrix(0, pn, 2)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- lm( target ~ dataset[, sel] + dataset[, mat[i, 1] ] )
            tab <- anova( ww )
            stat[i] <- tab[2, 4] 
            df1 <- tab[2, 1]   ;  df2 = tab[3, 1]
            pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
            mata[i, ] <- c(pval[i], stat[i]) 
          }
          stopCluster(cl)

        } else {  ## Robust
          do = 2 * as.numeric( logLik( moda[[ 1 ]] ) )
          fr = length( coef( moda[[ 1 ]] ) )
          
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata <- matrix(0, pn, 2)
          mod <- foreach( i = 1:pn, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
            # ww <- lmRob( target ~ dataset[, sel] + dataset[, mat[i, 1] ], control = cont )
            # tab <- aovlmrob( ww )
            # stat[i] <- tab[3, 2] 
            # pval[i] <- tab[3, 3]
            # mata[i, ] <- c(pval[i], stat[i]) 
            ww <- rlm( target ~ dataset[, sel] + dataset[, mat[i, 1] ], maxit = 2000 )
            mata[i, ] = c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) )
          }
          stopCluster(cl)
          
          difa = abs( mod[, 1] - do )
          pval = pchisq(difa, mod[, 2] - fr, lower.tail = FALSE, log.p = TRUE)
          mod = cbind( pval, difa)
        } 
        mat <- cbind(mat[, 1], mod)   

      }

    }

      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]     
      
        if ( stopping == "adjrsq" ) {
          if ( mat[ina, 2] < threshold ) {
            if ( robust == FALSE ) {
              mi = lm( target ~ dataset[, sela] + dataset[, sel] )
              tool[k] <- as.numeric( summary( mi )[[ 9 ]] )
            } else {
              # mi = robust::lmRob( target ~ dataset[, sela] + dataset[, sel], control = cont )
              # r2 <- as.numeric( summary( mi )[[ 8 ]] )               
              mi = MASS::rlm(target ~ dataset[, sela] + dataset[, sel], maxit = 2000)
              r2 = cor( target, fitted(mi) )^2
              tool[k] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(mi) ) - 1)
            }

            if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
              info <- rbind(info, c( 1e300, 0, 0 ) )
            
            } else {  
              info <- rbind(info, mat[ina, ] )
              sela <- info[, 1]
              mat <- mat[-ina , ] 
              if ( is.matrix(mat) ) {
                mat <- mat[ order( mat[, 2] ), ]
              }

              moda[[ k ]] <- mi
            }

          } else {
           info <- rbind(info, c( 1e300, 0, 0 ) )
          }

        } else if ( stopping == "BIC" ) {
          if ( mat[ina, 2] < threshold ) {
            if ( robust == FALSE ) {            
              mi = lm( target ~ dataset[, sela] + dataset[, sel] )
              tool[2] <- BIC( mi )
            } else {
              # mi = robust::lmRob( target ~ dataset[, sela] + dataset[, sel], control = cont )
              mi = MASS::rlm(target ~ dataset[, sela] + dataset[, sel], maxit = 2000)
              tool[2] = BIC(mi)
            }
            if ( tool[ k - 1] - tool[ k ] <= tol ) {
              info <- rbind(info, c( 1e300, 0, 0 ) )
            
            } else {  
              info <- rbind(info, mat[ina, ] )
              sela <- info[, 1]
              mat <- mat[-ina , ]
              if ( !is.matrix(mat) ) {
                mat <- matrix(mat, ncol = 3) 
              }
              mat <- mat[ order( mat[, 2] ), ]

              moda[[ k ]] <- mi
            }
          } else {
            info <- rbind(info, c( 1e300, 0, 0 ) )
          }
        }

      if ( sum( info[2, ] - c(1e300, 0, 0) ) == 0 ) {
        info <- matrix( info[1, ], ncol = 3)    
      }
 
     ######
     ###### k greater than 2
     ######

    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( abs( tool[ k ] - tool[ k - 1 ] ) > tol ) & ( nrow(mat) > 0 ) )  {
         
        k <- k + 1   
        pn <- p - k + 1 
        
        if ( ncores <= 1 ) {
          if ( robust == FALSE ) {
            for ( i in 1:pn ) {
              ww <- lm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ) )
              tab <- anova( ww )
              mat[i, 3] <- tab[k, 4] 
              df1 <- tab[k, 1]   ;  df2 = tab[k + 1, 1]
              mat[i, 2] <- pf( mat[i, 3], df1, df2, lower.tail = FALSE, log.p = TRUE )
            }
          } else {
            
            do = 2 * as.numeric( logLik( moda[[ k - 1 ]] ) )
            fr = length( coef( moda[[ k - 1 ]] ) )
            sta = dof = numeric(pn)
            
              for (i in 1:pn) {
                # ww = robust::lmRob( target ~ dataset[, sel] + dataset[, mat[i, 1] ], control = cont )
                # tab = aovlmrob( ww )
                # mat[i, 3] = tab[3, 2] 
                # mat[i, 2] = tab[3, 3]
                ww = MASS::rlm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), maxit = 2000 )
                sta[i] = 2 * as.numeric( logLik(ww) )
                dof[i] = length( coef(ww) )
              } 
              mat[, 3] = abs( sta - do )
              mat[, 2] = pchisq(mat[, 3], dof - fr, lower.tail = FALSE, log.p = TRUE)
            }

        } else {
          if ( robust == FALSE ) {  ## Non robust
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            stat <- pval <- numeric(pn)
            mata <- matrix(0, pn, 2)
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- lm( target ~., data = as.data.frame( dataset[, c(sela, mat[ i, 1]) ] ) )
              tab <- anova(ww)
              stat[i] <- tab[k, 4] 
              df1 <- tab[k, 1]   ;  df2 = tab[k + 1, 1]
              pval[i] <- pf( stat[i], df1, df2, lower.tail = FALSE, log.p = TRUE )
              mata[i, ] <- c( pval[i], stat[i] ) 
            }
            stopCluster(cl)

          } else {  ## Robust
            do = 2 * as.numeric( logLik( moda[[ k - 1 ]] ) )
            fr = length( coef( moda[[ k - 1 ]] ) )
            
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mata <- matrix(0, pn, 2)
            mod <- foreach( i = 1:pn, .combine = rbind, .export = c("rlm"), .packages = "MASS" ) %dopar% {
              # ww <- robust::lmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), control = cont )
              # tab <- aovlmrob(ww)
              # stat[i] <- tab[k + 1, 2] 
              # pval[i] <- tab[k + 1, 3]
              # mata[i, ] <- c( pval[i], stat[i] ) 
              ww <- rlm( target ~., data = as.data.frame( dataset[, c(sela, mat[ i, 1]) ] ), maxit = 2000 )
              mata[i, ] = c( 2 * as.numeric( logLik(ww) ), length( coef(ww) ) )
            }
            stopCluster(cl)  
            difa = abs( mod[, 1] - do )
            pval = pchisq(difa, mod[, 2] - fr, lower.tail = FALSE, log.p = TRUE)
            mod = cbind( pval, difa)
          }

          mat <- cbind( mat[, 1], mod )   

        }
       
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]   
        
          if ( stopping == "BIC" ) {
            if ( mat[ina, 2] < threshold ) {
              if ( robust == FALSE ) {
                ma = lm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ) )
                tool[k] <- BIC( ma )
              } else {
                # ma = robust::lmRob( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), control = cont )
                ma = MASS::rlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), maxit = 2000 )
                tool[k] =  BIC(ma)
              }

              if ( tool[ k - 1] - tool[ k ] <= tol ) {
                info <- rbind(info, c( 1e300, 0, 0 ) )
              
              } else { 
                info <- rbind( info, mat[ina, ] )
                sela <- info[, 1]
                mat <- mat[-ina , ] 
                if ( !is.matrix(mat) ) {
                  mat <- matrix(mat, ncol = 3) 
                }
                mat <- mat[ order( mat[, 2] ), ]

                moda[[ k ]] <- ma
              }
             
            } else {
              info <- rbind(info, c( 1e300, 0, 0 ) )
            }

          } else if ( stopping == "adjrsq" ) {
            if ( mat[ina, 2] < threshold ) {

              if ( robust == FALSE ) {
                ma = lm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ) )
                tool[k] <- as.numeric( summary(ma)[[ 9 ]] )
              } else {
                # ma = robust::lmRob( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), control = cont )
                # r2 <- as.numeric( summary( ma )[[ 8 ]] ) 
                ma = MASS::rlm( target ~., data = as.data.frame( dataset[, c(sela, sel)] ), maxit = 2000 )
                r2 = cor(target, fitted(ma) )^2
                tool[k] = 1 - (1 - r2) * (n - 1) / ( n - length( coef(ma) ) - 1)
              }               
              if ( tool[ k ] - tool[ k - 1 ] <= tol ) {
                info <- rbind(info, c( 1e300, 0, 0 ) )
              
              } else { 
                info <- rbind( info, mat[ina, ] )
                sela <- info[, 1]
                mat <- mat[-ina , ]
                if ( is.matrix(mat) ) {
                  mat <- mat[ order( mat[, 2] ), ]
                }

                moda[[ k ]] <- ma
              }

            } else {
              info <- rbind(info, c( 1e300, 0, 0 ) )
            }
        
          }
      }

    }

    runtime <- proc.time() - runtime

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

          if ( robust == FALSE ) {
            models[[ 1 ]] <- final <- lm( target ~., data = as.data.frame( xx ) )
          } else {
            # models[[ 1]] <- final <- robust::lmRob( target ~., data = as.data.frame( xx ), control = cont )
            models[[ 1 ]] <- final <- MASS::rlm( target ~., data = as.data.frame( xx ), maxit = 2000 )
          }
        
      } else {
        for (i in 1:d) {
          if ( robust == FALSE ) {
            models[[ i ]] <- lm( target ~., data = as.data.frame( xx[, 1:i] ) )
          } else { 
            # models[[ i ]] <- lmRob( target ~., data = as.data.frame( xx[, 1:i]), control = cont )
            models[[ i ]] <- MASS::rlm( target ~., data = data.frame( xx[, 1:i]), maxit = 2000 )
            
          }
        }
      }

      final <- summary( models[[ d ]] )

      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "p-value", "stat", stopping )
      rownames(info) <- info[, 1]
    }

    list(mat = t(mat), info = info, models = models, final = final, runtime = runtime ) 
}
  
glm.fsreg <- function(target, dataset, threshold = 0.05, tol = 2, robust = FALSE, ncores = 1, maxit = 100 ) {
  
  ## target can be Real valued (normal), binary (binomial) or counts (poisson)
  ## dataset is a matrix or a data.frame with the predictor variables
  ## test is the test used, but not really required because depending on the data this will be decided.
  ## there is no hrm is psecifying this though
  ## threshold is the level of significance
  ## method can be either BIC or adjrsq (for non-robust linear models only). The BIC is a consistent method for selecting
  ## models, thus we also use it to avoid overfitting
  ## stopping is based on "BIC"
  ## tol is the tolerance value for the method. If BIC is used as the stopping rule, the default is 2, but usually can be 2 or 4.
  ## If BIC is used as a way to proceed, the tol is 0.
  ## robust is for robust modelling. TRUE or FALSE
  ## ncores is for parallel processing 
  
  threshold <- log(threshold)  ## log of the significance level
  
  p <- ncol(dataset)  ## number of variables
  devi <- dof <- numeric( p )  
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  tool <- numeric( length( min(n, p) ) )

  #########
  ## if it is binomial or poisson regression
  #########
 
  
    if ( length( unique(target) ) == 2 ) {
      oiko <- "binomial"  ## binomial regression
    } else {
      oiko <- "poisson"  ## poisson regression
    }

    runtime <- proc.time()
    
    devi = dof = numeric(p)
    if ( robust == FALSE ) {
      ini = glm( target ~ 1, family = oiko )$deviance  ## residual deviance
    } else {
      ini = robust::glmRob( target ~ 1, family = oiko, maxit = maxit )$deviance  ## residual deviance
    }
 
    if (ncores <= 1) {
      if ( robust == FALSE ) {  ## Non robust
        for (i in 1:p) {
          mi <- glm( target ~ dataset[, i], family = oiko )
          devi[i] <- mi$deviance
          dof[i] = length( coef( mi ) ) 
        }
        
      } else {  ## Robust
        for (i in 1:p) { 
          mi <- robust::glmRob( target ~ dataset[, i], family = oiko, maxit = maxit )
          devi[i] <- mi$deviance
          dof[i] = length( coef( mi ) )          
        }
      }

      stat = ini - devi
      pval = pchisq( stat, dof - 1, lower.tail = FALSE, log.p = TRUE )

    } else {
      if ( robust == FALSE ) {  ## Non robust
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mata <- matrix(0, p, 2)
        mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
          ww <- glm( target ~ dataset[, i], family = oiko )
          mata[i, ] <- c( ww$deviance, length( coef( ww ) )  )
        }

        stopCluster(cl)

      } else {  ## Robust
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mata <- matrix(0, p, 2)
        mod <- foreach( i = 1:p, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
          ww <- glmRob( target ~ dataset[, i], family = oiko, maxit = maxit  )
          mata[i, ] <- c( ww$deviance, length( coef( ww ) )  )
        }

        stopCluster(cl)
      }

      stat = ini - mod[, 1]
      pval = pchisq( stat, mod[, 2] - 1, lower.tail = FALSE, log.p = TRUE )
    }

    mat <- cbind(1:p, pval, stat) 
   
    colnames(mat)[1] <- "variables"
    rownames(mat) <- 1:p

    sel <- which.min(pval)
    info <- matrix( numeric(3), ncol = 3 )
    sela <- sel
    
    if ( mat[sel, 2] < threshold ) {
      info[1, ] <- mat[sel, ]
      mat <- mat[-sel, ] 
      if ( !is.matrix(mat) ) {
        mat <- matrix(mat, ncol = 3) 
      }
      mat <- mat[ order( mat[, 2] ), ]

      if ( robust == FALSE ) {
        mi <- glm( target ~ dataset[, sel], family = oiko )
        tool[1] <- BIC( mi )
      } else {
        mi <- robust::glmRob( target ~ dataset[, sel], family = oiko, maxit = maxit )
        tool[1] <- mi$deviance + length( coef( mi ) ) * log(n)
      }
      moda[[ 1 ]] <- mi
    }

    ############
    ###       k equals 2
    ############ 
    
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- 2
      pn <- p - k + 1   

      ini = moda[[ 1 ]]$deviance  ## residual deviance
      do = length( coef( moda[[ 1 ]]  ) ) 

      if ( ncores <= 1 ) {
        devi = dof = numeric(pn)
        if ( robust == FALSE ) {  ## Non robust
          for ( i in 1:pn ) {
            ww <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko )
            devi[i] <- ww$deviance
            dof[i] = length( coef( ww ) )          
          }

        } else {  ## Robust
          for ( i in 1:pn ) {
            ww <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
            devi[i] <- ww$deviance
            dof[i] = length( coef( ww ) ) 
          }   
        }
        stat = ini - devi
        pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )

      } else {
        
        if ( robust == FALSE ) {  ## Non robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata = matrix(0, pn, 2)  
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko )
            mata[i, ] <- c( ww$deviance, length( coef( ww ) ) )
          }

          stopCluster(cl)

        } else {  ## Robust
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata = matrix(0, pn, 2)  
          mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
            ww <- glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
            mata[i, ] <- c( ww$deviance, length( coef( ww ) ) )
          }

          stopCluster(cl)
        }

        stat = ini - mod[, 1]
        pval = pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )

      }

      mat[, 2:3] <- cbind(pval, stat)
  
      ina <- which.min(mat[, 2])
      sel <- mat[ina, 1]    
      
      if ( mat[ina, 2] < threshold ) {
        if ( robust == FALSE ) {
          ma <- glm( target ~ dataset[, sela] + dataset[, sel], family = oiko )
          tool[k] <- BIC( ma )
        } else {
          ma <- robust::glmRob( target ~  dataset[, sela] + dataset[, sel], family = oiko, maxit = maxit )
          tool[k] <- ma$deviance + length( coef( ma ) ) * log(n)
        }

        if ( tool[ k - 1 ] - tool[ k ] <= tol ) {
          info <- rbind(info, c( 1e300, 0, 0 ) )
          
        } else { 
          info <- rbind(info, c( mat[ina, ] ) )
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
    }


    ############
    ###       k greater than 2
    ############ 


    if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
      while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) )  {
        
        ini = moda[[ k ]]$deviance  ## residual deviance
        do = length( coef( moda[[ k ]]  ) ) 

        k <- k + 1   
        pn <- p - k  + 1
        
        if (ncores <= 1) {  
          devi = dof = numeric(pn) 
          if ( robust == FALSE ) {  ## Non robust
            for ( i in 1:pn ) {
              ma <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1] ) ] ), family = oiko )
              devi[i] <- ma$deviance
              dof[i] = length( coef( ma ) ) 
            }

          } else {  ## Robust
            for ( i in 1:pn ) {
              ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
              devi[i] <- ma$deviance
              dof[i] = length( coef( ma ) ) 
            }
          }

          stat = ini - devi
          pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )

        } else {
          if ( robust == FALSE ) {  ## Non robust
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            devi = dof = numeric(pn)
            mata = matrix(0, pn, 2)  
            mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
              ww <- glm( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko )
              mata[i, ] <- c( ww$deviance, length( coef( ww ) ) )
            }

           stopCluster(cl)

          } else {  ## Robust
            cl <- makePSOCKcluster(ncores)
            registerDoParallel(cl)
            mata = matrix(0, pn, 2)  
            mod <- foreach( i = 1:pn, .combine = rbind, .export = "glmRob", .packages = "robust" ) %dopar% {
              ww <- glmRob( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ), family = oiko, maxit = maxit )
              mata[i, ] <- c( ww$deviance, length( coef( ww ) ) )
            }

            stopCluster(cl)
          }

          stat = ini - mod[, 1]
          pval = pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )

        }

       mat[, 2:3] <- cbind(pval, stat)

       ina <- which.min(mat[, 2])
       sel <- mat[ina, 1]    
        
          if ( mat[ina, 2] < threshold ) {
            if ( robust == FALSE ) {
              ma <- glm( target ~., data = as.data.frame( dataset[, c(sela, sel) ] ), family = oiko )
              tool[k] <- BIC( ma )
            } else {
              ma <- robust::glmRob( target ~., data = as.data.frame( dataset[, c(sela, sel) ] ), family = oiko, maxit = maxit )
              tool[k] <- ma$deviance + length( coef( ma ) ) * log(n)
            } 
            if ( tool[ k - 1 ] - tool[ k  ] < tol ) {
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
            models[[ 1 ]] <- final <- glm( target ~., data = as.data.frame( xx ), family = oiko )
          } else {
            models[[ 1]] <- final <- robust::lmRob( target ~., data = as.data.frame( xx ), family = oiko, maxit = maxit )
          }
        
      } else {
        for (i in 1:d) {
          if ( robust == FALSE ) {
            models[[ i ]] <- glm( target ~., data = as.data.frame( xx[, 1:i] ), family = oiko )
          } else { 
            models[[ i ]] <- glmRob( target ~., data = as.data.frame( xx[, 1:i]), family = oiko, maxit = maxit )
          }
        }
      }

      final <- summary( models[[ d ]] )

      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "p-value", "stat", "BIC" )
      rownames(info) <- info[, 1]
    }

    

    list(mat = t(mat), info = info, models = models, final = final, runtime = runtime ) 
    
}    
  
  
  
  
  
  
  
  
  
  
  
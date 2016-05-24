fs.reg <- function(target, dataset, threshold = 0.05, test = NULL, stopping = "BIC", tol = 2, robust = FALSE, ncores = 1, maxit = 100 ) {
  
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
  info <- matrix( 0, ncol = 2 )
  result = NULL
  con = log(n)
  
  #check for NA values in the dataset and replace them with the variable median or the mode
  if(any(is.na(dataset)) == TRUE)
  {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    
    if(class(dataset) == "matrix")
    {
      dataset = apply( dataset, 2, function(x){ x[ which(is.na(x)) ] = median(x,na.rm = TRUE) } );
    }else{
      for(i in 1:ncol(dataset))
      {
        if ( any( is.na(dataset[, i]) ) )
        {
          xi = dataset[,i]
          if(class(xi) == "numeric")
          {                    
            xi[ which(is.na(xi)) ] = median(xi,na.rm = TRUE) 
          } else if( class(xi) == "factor" ){
            xi[ which(is.na(xi)) ] = levels(xi)[ which.max(xi) ]
          }
          dataset[, i] = xi
        }
      }
    }
  }
  
  ##################################
  # target checking and initialize #
  ##################################
  
  
  
  dataset <- as.data.frame(dataset)  ## just in case
  if ( is.null( colnames(dataset) ) )  {  ## checks for column names
    colnames(x) <- paste("X", 1:p, sep = "")
  }	
  
  ## dependent (target) variable checking if no test was given, 
  ## but other arguments are given. For some cases, these are default cases
  
  if ( is.null(test) ) {
    
    ## multivariate data
    if ( sum( class(target) == "matrix" ) > 0 ) {
      test <- "gaussian"
      a <- rowSums(target)
      if ( min( target ) > 0 & round( sd(a), 16 ) == 0 ) { ## are they compositional data?
        target <- log( target[, -1] / target[, 1] )
      }
    }
    
    ## percentages
    if ( min( target ) > 0 & max( target ) < 1 )  {  ## are they percentages?
      target <- log( target / (1 - target) ) 
      test <- "gaussian"
    }
    
    ## surival data
    if ( sum( class(target) == "Surv" ) > 0 ) {
      test <- "Cox"
    }
    
    ## binary data
    if ( length( unique(target) ) == 2 ) {
      test <- "binary"   
    }
    
    ## ordinal, multinomial or perhaps binary data
    if ( is.factor(target) ) {
      if ( !is.ordered(target) ) {
        if ( length(unique(target) ) == 2 ) {
          target <- as.vector(target)
          test <- "binary"
        } else {
          test <- "multinomial"
        }  
        
      } else {
        if ( length(unique(target) ) == 2 ) {
          target <- as.vector(target)
          test <- "binary"
        } else {
          test <- "ordinal"    
        }
      }
    }
    
    ## count data
    if ( is.vector(target) ) {
      
      if ( length( unique(target) ) > 2  &  sum( round(target) - target ) == 0 ) {
        test <- "poisson"
    
        ## binomial regression 
      
      } else if ( length( unique(target) ) == 2 ) {
        test <- "binary"

        ## linear regression 
      } else if ( sum( class(target) == "numeric" ) > 0 ) {
        test <- "gaussian"  
      }
    }  
    
  }
  
  #available conditional independence tests
  av_models = c("gaussian", "median", "beta", "Cox", "Weibull", "binary", "multinomial", "ordinal", "poisson", "nb", "zip", "speedglm");
  
  #cat(test)
  
  if ( test == "binary" || test == "poisson" ) {
    
    result <- glm.fsreg( target, dataset, threshold = exp(threshold), tol = tol, robust = robust, ncores = ncores, maxit = maxit ) 
    
  } else if ( test == "gaussian"  &  !is.matrix(target) ) {
    
    result <- lm.fsreg( target, dataset, threshold = exp(threshold), stopping = stopping, tol = tol, robust = robust, ncores = ncores ) 
    
    
  } else {
    
    test = match.arg(test , av_models ,TRUE);
    #convert to closure type
    
    if ( test == "median" ) {
      test = rq
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "beta" ) {
      test = betareg 
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "Cox" ) {
      test = coxph 
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "multinomial" ) {
      test = multinom 
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "ordinal" ) {
      test = clm 
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "neg.bin" ) {
      test = glm.nb
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "zip" ) {
      test = zeroinfl 
      robust = FALSE
      stopping = "BIC"
      
    } else if ( test == "gaussian"  &  is.matrix(target) ) {
      test = lm 
      robust = FALSE
    }
    
    
    runtime <- proc.time()
    
    devi = dof = numeric(p)
    ini = test( target ~ 1 ) 
    ini =  2 * as.numeric( logLik(ini) )  ## initial 
    
    if (ncores <= 1) {
      for (i in 1:p) {
        mi <- test( target ~ dataset[, i] )
        devi[i] <-  2 * as.numeric( logLik(mi) )
        dof[i] = length( coef( mi ) ) 
      }
      
      stat = abs( devi - ini )
      pval = pchisq( stat, dof - 1, lower.tail = FALSE, log.p = TRUE )
      
    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      mata <- matrix(0, p, 2)
      mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
        ww <- test( target ~ dataset[, i] )
        mata[i, ] <- c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) )  )
      }
      
      stopCluster(cl)
      
      stat =  abs( mod[, 1] - ini )
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
      
      mi <- test( target ~ dataset[, sel] )
      tool[1] <-  - 2 * as.numeric( logLik(mi) ) + length( coef(mi) ) * con
      
      moda[[ 1 ]] <- mi
    }
    
    ############
    ###       k equals 2
    ############ 
    
    if ( info[k, 2] < threshold  &  nrow(mat) > 0 ) {
      
      k <- 2
      pn <- p - k + 1   

      ini =  2 * as.numeric( logLik( moda[[ 1 ]] ) ) 
      do = length( coef( moda[[ 1 ]]  ) ) 
      
      if ( ncores <= 1 ) {
        devi = dof = numeric(pn)
        for ( i in 1:pn ) {
          ww <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ) )
          devi[i] <-  2 * as.numeric( logLik(ww) )
          dof[i] = length( coef( ww ) )          
        }
        
        stat = abs( devi - ini )
        pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
        
        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        mata = matrix(0, pn, 2)  
        mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
          ww <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ) )
          mata[i, ] <- c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) )  )
        }
        
        stopCluster(cl)
        
        stat = abs( mod[, 1] - ini )
        pval = pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
        
      }

    }
    
    mat[, 2:3] <- cbind(pval, stat)
    
    ina <- which.min(mat[, 2])
    sel <- mat[ina, 1]    
    
    if ( mat[ina, 2] < threshold ) {
      ma <- test( target ~ dataset[, sela] + dataset[, sel] )
      tool[2] <-  - 2 * as.numeric( logLik(ma) ) + length( coef(ma) ) * con
      
      
      if ( tool[ 1 ] - tool[ 2 ] <= tol ) {
        info <- rbind(info, c( 1e300, 0, 0 ) )
        
      } else { 
        info <- rbind(info, c( mat[ina, ] ) )
        sela <- info[, 1]
        mat <- mat[-ina , ] 
        if ( !is.matrix(mat) ) {
          mat <- matrix(mat, ncol = 3) 
        }
        mat <- mat[ order( mat[, 2] ), ]
        
        moda[[ 2 ]] <- ma
      }
      
    } else {
      info <- rbind(info, c( 1e300, 0, 0 ) )
    }
    
  
  ############
  ###       k greater than 2
  ############ 
  
  
  if ( nrow(info) > 1  &  nrow(mat) > 0 ) {
    while ( ( info[k, 2] < threshold ) &  ( k < n ) & ( tool[ k - 1 ] - tool[ k ] > tol ) & ( nrow(mat) > 0 ) )  {
      
      ini =  2 * as.numeric( logLik( moda[[ k ]] ) ) 
      do = length( coef( moda[[ k ]]  ) ) 
      
      k <- k + 1   
      pn <- p - k  + 1
      
      if (ncores <= 1) {  
        devi = dof = numeric(pn) 
        for ( i in 1:pn ) {
          ma <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1] ) ] ) )
          devi[i] <-  2 * as.numeric( logLik(ma) ) 
          dof[i] = length( coef( ma ) ) 
        }
        
        stat = abs( devi - ini )
        pval = pchisq( stat, dof - do, lower.tail = FALSE, log.p = TRUE )
        
      } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          mata = matrix(0, pn, 2)  
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ) )
            mata[i, ] <- c( 2 * as.numeric( logLik(ww) ), length( coef( ww ) )  )
          }
          
          stopCluster(cl)
          
          stat = abs( mod[, 1] - ini )
          pval = pchisq( stat, mod[, 2] - do, lower.tail = FALSE, log.p = TRUE )
          
        }
        
        mat[, 2:3] <- cbind(pval, stat)
        
        ina <- which.min(mat[, 2])
        sel <- mat[ina, 1]    
        
        if ( mat[ina, 2] < threshold ) {
            ma <- test( target ~., data = as.data.frame( dataset[, c(sela, sel) ] ) )
            tool[k] <-  - 2 * as.numeric( logLik(ma) ) + length( coef(ma) ) * con
 
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
        
        models[[ 1 ]] <- test( target ~., data = data.frame( xx ) )

        
      } else {
        for (i in 1:d) {
          models[[ i ]] <- test( target ~., data = data.frame( xx[, 1:i] ) )
        }
      }
      
      final <- summary( models[[ d ]] )
      
      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      info <- cbind( info, tool[ 1:d ] ) 
      colnames(info) <- c( "variables", "p-value", "stat", "BIC" )
      rownames(info) <- info[, 1]
    }
    
    result = list(mat = t(mat), info = info, models = models, final = final, runtime = runtime ) 
    
  }         
  
  result
  
}    











bic.fsreg <- function( target, dataset, test = NULL, robust = FALSE, tol = 0, ncores = 1, maxit = 100 ) {

  p <- ncol(dataset)  ## number of variables
  bico <- numeric( p )
  moda <- list()
  k <- 1   ## counter
  n <- length(target)  ## sample size
  con = log(n)
  tool <- NULL
  info <- matrix( 0, ncol = 2 )
  result = NULL

  #check for NA values in the dataset and replace them with the variable median or the mode
  if(any(is.na(dataset)) == TRUE)
  {
    #dataset = as.matrix(dataset);
    warning("The dataset contains missing values (NA) and they were replaced automatically by the variable (column) median (for numeric) or by the most frequent level (mode) if the variable is factor")
    
    if(class(dataset) == "matrix")
    {
      dataset = apply(dataset, 2, function(x){x[which(is.na(x))] = median(x,na.rm = TRUE)});
    }else{
      for(i in 1:ncol(dataset))
      {
        if(any(is.na(dataset[,i])))
        {
          xi = dataset[,i]
          if(class(xi) == "numeric")
          {                    
            xi[which(is.na(xi))] = median(xi,na.rm = TRUE) 
          }else if(class(xi) == "factor"){
            xi[which(is.na(xi))] = levels(xi)[which.max(xi)]
          }
          dataset[,i] = xi
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

    ## linear regression 
    if ( class(target) == "numeric" ) {
      test <- "gaussian"  
    }
    
    ## multivariate data
    if ( class(target) == "matrix" ) {
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
    if ( class(target) == "Surv" ) {
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
    
    ci_model = test
    #cat(test)
    
  if ( test == "binary" || test == "poisson" ||  ( test == "gaussian"  &  !is.matrix(target) ) ) {
   
    result <- bic.glm.fsreg( target, dataset, robust = robust, tol = tol, ncores = ncores, maxit = maxit ) 

  } else {
 
    test = match.arg(test , av_models ,TRUE);
    #convert to closure type
    
     if ( test == "median" ) {
      test = rq
      robust = FALSE
    
    } else if ( test == "beta" ) {
      test = betareg 
      robust = FALSE

    } else if ( test == "Cox" ) {
      test = coxph 
      robust = FALSE

    } else if ( test == "multinomial" ) {
      test = multinom 
      robust = FALSE

    } else if ( test == "ordinal" ) {
      test = clm 
      robust = FALSE

    } else if ( test == "neg.bin" ) {
      test = glm.nb
      robust = FALSE

    } else if ( test == "zip" ) {
      test = zeroinfl 
      robust = FALSE

    } else if ( test == "gaussian"  &  is.matrix(target) ) {
      test = lm 
      robust = FALSE
    }

    runtime <- proc.time()

      ini = test( target ~ 1 )
      ini =  - 2 * as.numeric( logLik(ini) ) + con  ## initial BIC
    
    if (ncores <= 1) {
        for (i in 1:p) {
          mi <- test( target ~ dataset[, i] )
          bico[i] <-  - 2 * as.numeric( logLik(mi) ) + length( coef( mi ) ) * con
        }

      mat <- cbind(1:p, bico)

      if( any( is.na(mat) ) ) {
        mat[ which( is.na(mat) ) ] = ini
      }

    } else {
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      bico <- numeric(p)
      mod <- foreach( i = 1:p, .combine = rbind) %dopar% {
        ww <- test( target ~ dataset[, i] )
        bico[i] <-  - 2 * as.numeric( logLik(ww) ) + length( coef( ww ) ) * con
      }
      stopCluster(cl)

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
     
      mi = test( target ~ dataset[, sel] )
      tool[1] <-  - 2 * as.numeric( logLik(mi) ) + length( coef( mi ) ) * con
      if ( is.na(tool[1]) )  tool[1] <- ini

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
        for ( i in 1:pn ) {
          ma <- test( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ) )
          bico[i] <-  - 2 * as.numeric( logLik(ma) ) + length( coef( ma ) ) * con
        }

        mat[, 2] <- bico

      } else {

        cl <- makePSOCKcluster(ncores)
        registerDoParallel(cl)
        bico <- numeric(pn)
        mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
          ww <- test( target ~., data = as.data.frame( dataset[, c(sel, mat[i, 1]) ] ) )
          bico[i] <-  - 2 * as.numeric( logLik(ma) ) + length( coef( ma ) ) * con
        }

        stopCluster(cl)

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

        mi = test( target ~., data = as.data.frame( dataset[, sela] ) )
        tool[2] <-  - 2 * as.numeric( logLik(mi) ) + length( coef( mi ) ) * con
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
          for ( i in 1:pn ) {
            ma <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ) )
            mat[i, 2] <-  - 2 * as.numeric( logLik(ma) ) + length( coef( ma ) ) * con
          }

        } else {
          cl <- makePSOCKcluster(ncores)
          registerDoParallel(cl)
          bico <- numeric(pn)
          mod <- foreach( i = 1:pn, .combine = rbind) %dopar% {
            ww <- test( target ~., data = as.data.frame( dataset[, c(sela, mat[i, 1]) ] ) )
            bico[i] <- - 2 * as.numeric( logLik(ww) ) + length( coef( ww ) ) * con
          }

          stopCluster(cl)

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

          ma = test( target ~., data =as.data.frame( dataset[, sela] ) )
          tool[k] <-  - 2 * as.numeric( logLik(ma) ) + length( coef( ma ) ) * con

          moda[[ k ]] <- ma

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

        xx <- as.data.frame( dataset[, sela] )
        colnames(xx) <- paste("V", sela, sep = "")

        models <- final <- test( target ~., data = as.data.frame( xx ) )

      } else {
        for (i in 1:d) {
          models[[ i ]] <- test( target ~., data = as.data.frame( xx[, 1:i] ) )
        }

        final <- summary( models[[ d ]] )

      }

      info <- info[1:d, ]
      if ( d == 1 )  info <- matrix(info, nrow = 1)
      colnames(info) <- c( "variables", "BIC" )
      rownames(info) <- info[, 1]

    }

    result = list(mat = t(mat), info = info, models = models, final = final, runtime = runtime )

  } 

  result

}

##  ############################################################################

cv <- function( object, ... ) UseMethod( "cv" )

##  ############################################################################

# cv.default <- function( 
#   object, 
#   formula = NULL, subset = NULL,
#   nset = 10, seed = NULL, sets = NULL,
#   
#   ... )
# {
#   
#   
#   
# }

##  ############################################################################

cv.georob <-
  function(
    object,
    formula = NULL, subset = NULL,
    method = c( "block", "random" ), nset = 10, seed = NULL, sets = NULL,
    duplicates.in.same.set = TRUE,
    re.estimate = TRUE, param = object[["param"]], 
    fit.param = object[["initial.objects"]][["fit.param"]],
    aniso = object[["aniso"]][["aniso"]], 
    fit.aniso = object[["initial.objects"]][["fit.aniso"]],
    return.fit = FALSE, reduced.output = TRUE,
    lgn = FALSE,
    mfl.action = c( "offset", "stop" ),
    ncores = min( nset, detectCores() ),
    verbose = 0,
    ...
  )
{
  
  ## Function computes nset-fold cross-validation predictions from a
  ## fitted georob object
  
  
  ## History:
  
  ## 2011-10-24 Korrektur Ausschluss von nichtbenoetigten Variablen fuer lognormal kriging
  ## 2011-12-23 AP modified for replicated observations and for parallel computing
  ## 2012-03-02 AP eliminated possibility for logging to file in parallel processing
  ## 2012-03-19 AP correction of error in parallel processing on Windows
  ## 2012-05-01 AP correct handling of NAs
  ## 2012-05-04 AP modifications for lognormal block kriging
  ## 2012-05-09 AP correction of error if a new formula is passed via update to georob
  ## 2012-05-22 AP correction of error in passing param and fit.param to georob
  ## 2012-06-05 AP correction of error in handling optional subset argument
  ## 2012-11-04 AP handling compressed cov.betahat
  ## 2012-12-04 AP modifiction for changes in predict.georob
  ## 2013-04-24 AP changes for parallelization on windows os
  ## 2013-05-23 AP correct handling of missing observations
  ## 2013-05-24 AP separate initial variogram parameters for each cross-validation set
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2013-07-02 AP passing initial values of aniso and fit.aniso to georob via update
  ## 2013-07-05 AP return "variogram.model" as part of fit componnent
  ## 2014-02-21 AP catching problem when factor are very unbalanced
  ## 2014-05-15 AP changes for version 3 of RandomFields
  ## 2014-05-28 AP catching error when all variogram parameters are fixed
  ## 2014-08-18 AP changes for parallelized computations
  ## 2015-03-13 AP minor changes if models are not re-estimated
  ## 2015-03-16 AP items to compute huberized prediction error
  ## 2015-06-23 AP modifications for robust prediction of response
  ## 2015-07-17 AP excluding blocks of observations
  ## 2015-07-20 AP inactivation of modifications for robust prediction of response 
  ##               (variables: robust, se.signal, scld.res, resscl)
  ## 2015-08-28 AP computation of hessian suppressed
  ## 2015-11-26 AP catching errors occurring during parallel model fits
  
  
  ## auxiliary function that fits the model and computes the predictions of
  ## a cross-validation set
  
  f.aux <- function( 
    ..i.., object, formula, data, sets,  
    param, fit.param, aniso, fit.aniso, lgn, 
    verbose, ...
  ){  ## cv function
    
    if (verbose) cat( "\n\n  processing cross-validation set", ..i.., "\n" ) 
    
    ## change environment of terms and formula so that subset selection works for update            
    
    environment( formula ) <- environment()
    environment( object[["terms"]] ) <- environment()
    
    ## read-off initial values of variogram parameters
    
    if( ( is.matrix( param ) || is.data.frame( param ) ) )  param <- param[..i..,]
    if( ( is.matrix( fit.param ) || is.data.frame( fit.param ) ) )  fit.param <- fit.param[..i..,]
    if( ( is.matrix( aniso ) || is.data.frame( param ) ) )  aniso <- aniso[..i..,]
    if( ( is.matrix( fit.aniso ) || is.data.frame( fit.aniso ) ) )  fit.aniso <- fit.aniso[..i..,]

    t.georob <- update( 
      object, 
      formula = formula,
      data = data,
      subset = -sets[[..i..]] ,
      param = param, fit.param = fit.param,
      aniso = aniso, fit.aniso = fit.aniso,
      verbose = verbose,
      object. = NULL,
      ...
    )
    
    if( verbose > 0 ){
      cat( "\n\n" )
      print( summary( t.georob ) )
    }
    
    ## compute predictions for current set
    
    ## note that the signal is predicted here; the predictive distribution
    ## of the response is then modelled from the pooled cross-validation
    ## predictions of th signal and the residuals of the object fitted to
    ## the full data set (cf.  main part of the function)
    
    t.predict <- predict( 
      t.georob, newdata = data[sets[[..i..]], ], type = "signal", signif = NULL,
      mmax = length( sets[[..i..]] ),
      control = control.predict.georob( ncores = 1, extended.output = lgn ) 
    )
        
    ## backtransformation for log-normal kriging
    
    if( lgn ){
      t.predict <- lgnpp( t.predict )
      t.predict <- t.predict[, -match( 
        c( "trend", "var.pred", "cov.pred.target", "var.target" ), names( t.predict ) 
      )]
    }
    
    t.predict <- data.frame( i = sets[[..i..]], t.predict )
        
    if( reduced.output ){
      
      if( !is.null( t.georob[["cov"]][["cov.betahat"]] ) ){
        t.se.coef <- sqrt( diag( expand( t.georob[["cov"]][["cov.betahat"]] ) ) )
      } else {
        t.se.coef <- NULL
      }
      
      t.georob <- t.georob[c(  
        "tuning.psi", "converged", "convergence.code",
        "gradient", "variogram.model", "param", "aniso",
        "coefficients"
      )]
      
      t.georob[["aniso"]] <- t.georob[["aniso"]][["aniso"]]
      
      if( !is.null( t.se.coef ) ) t.georob[["se.coefficients"]] <- t.se.coef
      
    }
    
    return( list( pred = t.predict, fit = t.georob ) )
    ## end cv function
  }
  
  ## begin of main body of function
  
  method <- match.arg( method )
      
  mfl.action <- match.arg( mfl.action )
  
  ## update terms of object is formula is provided
  
  if( !is.null( formula ) ){
    formula <- update( formula( object ), formula )
    object[["terms"]] <- terms( formula )
  } else {
    formula <- formula( object )
  }
    
  ## get data.frame with required variables (note that the data.frame passed
  ## as data argument to georob must exist in GlobalEnv)

  data <- cbind(
    get_all_vars( 
      formula( object ), data = eval( getCall(object)[["data"]] )
    ),
    get_all_vars( 
      object[["locations.objects"]][["locations"]], eval( getCall(object)[["data"]] )
    )
  )
  
  if( identical( class( object[["na.action"]] ), "omit" ) ) data <- na.omit(data)
  
  ## select subset if appropriate
  
  if( !is.null( subset ) ){
    data <- data[subset, ]
    object[["Tmat"]] <- object[["Tmat"]][subset]
  } else if( !is.null( getCall(object)[["subset"]] ) ){
   data <- data[eval( getCall(object)[["subset"]] ), ]
  }
  
#   if( !is.null( getCall(object)[["subset"]] ) )
   
  ## define cross-validation sets
  
  if( is.null( sets ) ){
    
    if( !is.null( seed ) ) set.seed( seed )
    
    sets <- switch(
      method,
      random = {
        sets <- runif( NROW( data ) )
        sets <- cut( 
          sets, 
          breaks = c( -0.1, quantile( sets, probs = ( 1:(nset-1)/nset ) ), 1.1 )
        )
        factor( as.integer( sets ) )
      },
      block = {
        sets <- kmeans( object[["locations.objects"]][["coordinates"]], centers = nset )[["cluster"]]
        if( !is.null( subset ) ) sets <- sets[subset]
        sets
      }
    )
    
  } else {
    
    nset <- length( unique( sets ) )
    if( length( sets ) != NROW( data ) ) stop(
      "sets must be an integer vector with length equal to the number of observations"      
    )
    
  }
  
  if( duplicates.in.same.set ){
    dups <- duplicated( object[["Tmat"]] )
    idups <- match( object[["Tmat"]][dups], object[["Tmat"]][!dups] )
    sets[dups] <- (sets[!dups])[idups]
  }
  
  sets <- tapply(
    1:NROW( data ),
    sets,
    function( x ) x
  )
  
  ## check whether all levels of factors are present in all calibration sets
  
  isf <- sapply( data, is.factor )
  
  mfl <- sapply(
    names( data )[isf],
    function( v, data, sets ){
      nolevel <- sapply( 
        sets, 
        function(s, x) any( table( x[-s] ) < 1 ),
        x = data[, v]
      )
      if( any( nolevel ) ){
        switch(
          mfl.action,
          "stop" = stop(
            "for factor '", v, "' some levels are missing in some calibration set(s),\n",
            "  try defining other cross-validation sets to avoid this problem"
          ),
          "offset" = {
            warning(
              "for factor '", v, "' some levels are missing in some calibration set(s),\n",
              "  respective factors are treated as offset terms in model"
            )
            cat(
              "  for factor '", v, "' some levels are missing in some calibration set(s),\n",
              "  respective factors are treated as offset terms in model"
            )
          }
        )
        TRUE
      } else {
        FALSE
      }
    },
    data = data, sets = sets
  )
  
  if( any( mfl ) ){
    v <- names( mfl )[mfl]
    
    ## construct offset and append to data
    
    offset <- apply( predict( object, type = "terms", terms = v )$fit, 1, sum )
    
    if( length( offset ) != nrow( data ) ) stop( "offset and data with incompatible dimensions" )
    data[, "..offset.."] <- offset
    
    ## update formula
    formula <- update( 
      formula, 
      as.formula(
        paste( 
          ". ~ . -", paste( v, collapse = " - " ) , "+ offset(..offset..)" )
      ) 
    )
    
  }
  
  
  ## redefine na.action component of object
  
  if( identical( class( object[["na.action"]] ), "exclude" ) ){
    class( object[["na.action"]] ) <- "omit"
  }
  
  ## check whether all variogram parameters are fixed
  
  if( !any( c( fit.param, fit.aniso ) ) && re.estimate ){
    re.estimate <- FALSE
    cat(
      "re.estimate set equal to FALSE because all variogram parameters are fixed\n\n"    
    )
    warnings(
      "re.estimate set equal to FALSE because all variogram parameters are fixed"    
    )
  }
  
  ## check dimension of param, fit.param, aniso, fit.aniso
  
  if( ( is.matrix( param ) || is.data.frame( param ) ) && nrow( param )!= nset ) stop(
    "'param' must have 'nset' rows if it is a matrix or data frame"  
  )
    
  if( ( is.matrix( fit.param ) || is.data.frame( fit.param ) ) && nrow( fit.param )!= nset ) stop(
    "'fit.param' must have 'nset' rows if it is a matrix or data frame"  
  )
    
  if( ( is.matrix( aniso ) || is.data.frame( aniso ) ) && nrow( aniso )!= nset ) stop(
    "'aniso' must have 'nset' rows if it is a matrix or data frame"  
  )
    
  if( ( is.matrix( fit.aniso ) || is.data.frame( fit.aniso ) ) && nrow( fit.aniso )!= nset ) stop(
    "'fit.aniso' must have 'nset' rows if it is a matrix or data frame"  
  )
  
  ## keep all variogram parameters fixed for re.estimate == FALSE
  
  if( !re.estimate ){
	fit.param <- default.fit.param( variance = FALSE, nugget = FALSE, scale = FALSE )[names( param )]
	fit.aniso <- default.fit.aniso()
  }
  
  
  ## set hessian equal to FALSE in control argument of georob call and update
  ## call
  
  if( reduced.output || !return.fit ){
    
    cl <- object[["call"]]
    
    if( "control" %in% names(cl) ){
      
      ## georob called with control argument
      
      cl.control <- as.list( cl[["control"]] )
      cl <- cl[ -match( "control", names(cl) ) ]
      if( "hessian" %in% names(cl.control) ){
        cl.control["hessian"] <- list( hessian = FALSE )
      } else {
        cl.control <- c( cl.control, hessian = FALSE )
      }
      
    } else {
      
      ## georob called without control argument
      
      cl.control <- list( as.symbol("control.georob"), hessian = FALSE )
      
    }
    
    object[["call"]] <- as.call(  c( as.list(cl), control = as.call(cl.control) ) )
  }
  
  ## loop over all cross-validation sets
  
  if( ncores > 1 ){
    
    ## parallel processing    

    if( .Platform[["OS.type"]] == "windows" ){
      
      ## create a SNOW cluster on windows OS
      
      cl <- makePSOCKcluster( ncores, outfile = "")
      
      ## export required items to workers
      
      junk <- clusterEvalQ( cl, require( georob, quietly = TRUE ) )
      
      t.result <- try(
        parLapply(
          cl, 
          1:length( sets ),
          f.aux, 
          object = object,
          formula = formula,
          data = data,
          sets = sets,
          param = param, fit.param = fit.param,
          aniso = aniso, fit.aniso = fit.aniso,
          lgn = lgn, 
          verbose = verbose,
          ...
        ), silent = TRUE
      )
      
      stopCluster(cl)
      
    } else {
      
      ## fork child processes on non-windows OS
      
      t.result <- try(
        mclapply(
          1:length( sets ),
          f.aux, 
          object = object,
          formula = formula,
          data = data,
          sets = sets,
          param = param, fit.param = fit.param,
          aniso = aniso, fit.aniso = fit.aniso,
          lgn = lgn, 
          verbose = verbose,
          mc.cores = ncores,
          mc.allow.recursive = FALSE,
          ...
        ) 
        , silent = TRUE
      )
            
    }
    
    has.error <- sapply(
      t.result, function( x ) identical( class(x), "try-error" ) 
    )
    
    if( any( has.error ) ){
      cat( "\nerror(s) occurred when fitting model in parallel to cross-validation sets:\n\n" )
      sapply( t.result[has.error], cat)
      cat( "\nrun cross-validation with arguments 'ncores = 1' and 'verbose = 2'\n\n" ) 
      stop()
    
    }
    
  } else {
    
    ## sequential processing
    
    t.result <- lapply(
      1:length( sets ),
      f.aux, 
      object = object,
      formula = formula,
      data = data,
      sets = sets,
      param = param, fit.param = fit.param,
      aniso = aniso, fit.aniso = fit.aniso,
      lgn = lgn, 
      verbose = verbose,
      ...
    )
    
  }
  
  ## create single data frame with cross-validation results 
  
  result <- t.result[[1]][["pred"]]
  result[["subset"]] <- rep( 1, nrow( t.result[[1]][["pred"]] ) )
  
  for( t.i in 2:length( t.result ) ) {
    result <- rbind( 
      result, 
      data.frame( 
        t.result[[t.i]][["pred"]],
        subset = rep( t.i, nrow( t.result[[t.i]][["pred"]] ) )
      )
    )
  }
  t.ix <- sort( result[["i"]], index.return = T )[["ix"]]
  result <- result[t.ix, ]
  result[["data"]] <- model.response( 
    model.frame( formula( object), data, na.action = na.pass ) 
  )
  
  if( lgn ) result[["lgn.data"]] <- exp( result[["data"]] )
  
  result <- result[, -match("i", colnames( result) )]
  
  isubset <- match( "subset", colnames( result ) )
  idata <- grep( "data", colnames( result ), fixed = TRUE )
  ipred<-  grep( "pred", colnames( result ), fixed = TRUE )
  ise <-   grep( "se", colnames( result ), fixed = TRUE )
  ise <- ise[ise != isubset]
  
  result <- cbind(
    result[, -c(isubset, idata, ipred, ise)],
    result[,  c(isubset, idata, ipred, ise)]
  )
  
  ## compute prediction standard errors of response 
  
  #   if( object[["tuning.psi"]] < object[["control"]][["tuning.psi.nr"]] ){
  #     
  #     ## robust
  #     
  #     resscl <- 1.
  #     warning( "scale factor for computing empirical distribution of residuals equals 1" )
  #     scld.res <- object[["residuals"]] / resscl
  #     se.signal <- result[, "se"]
  #     result[, "se"] <- sqrt( 
  #       result[, "se"]^2 + var(scld.res) * (length(scld.res) - 1) / length(scld.res) 
  #     )
  #   } else {
  
    ## Gaussian 
    
  #     scld.res <- NULL
  #     se.signal <- NULL
  result[, "se"] <- sqrt( result[, "se"]^2 + object[["param"]]["nugget"] )
    
  #   }
  
  ## prepare model fits for return
  
  t.fit <- lapply( t.result, function( x ) return( x[["fit"]] ) )
  
  if( re.estimate && !all( sapply( t.fit, function(x) x[["converged"]] ) ) ) warning(
    "lack of covergence for  ", 
    sum( !sapply( t.fit, function(x) x[["converged"]] ) ), " cross-validation sets"
  )
  
  result <- list( 
    pred = result, 
    fit = if( return.fit ) t.fit else NULL
  )
  
  #   attr( result[["pred"]], "nugget" )     <- sapply( t.fit, function(x) x[["param"]]["nugget"] )
  #   attr( result[["pred"]], "psi.func" )   <- object[["control"]][["psi.func"]]
  attr( result[["pred"]], "tuning.psi" ) <- object[["tuning.psi"]]
  #   attr( result[["pred"]], "exp.gauss.dpsi" )       <- object[["expectations"]][["exp.gauss.dpsi"]]
  #   attr( result[["pred"]], "se.signal" )  <- se.signal
  #   attr( result[["pred"]], "scaled.residuals" )   <- scld.res
  
  class( result ) <- "cv.georob"
  
  invisible( result )
  
}

##  ###########################################################################

plot.cv.georob <-
  function( 
    x, type = c( "sc", "lgn.sc", "ta", "qq", "hist.pit", "ecdf.pit", "mc", "bs" ), 
    smooth = TRUE, span = 2/3,
    ncutoff = NULL, 
    add = FALSE, 
    col, pch, lty,
    main, xlab, ylab, 
    ... 
  )
{
  
  ## plot method for class "cv.georob"  
  
  ## 2011-12-21 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-03-12 AP adding smooth curve of types sc and lgn.sc
  ## 2015-06-25 AP new method to compute pit, mc, bs and crps (Gaussian and robust)
  ## 2016-02-29 AP minor changes for adding plots to existing graphics
  
  x <- x[["pred"]]
  
  type = match.arg( type )
  
  if( type == "sc.lgn" && !"lgn.pred" %in% names( x ) ) stop(
    "lognormal kriging results missing, use 'lgn = TRUE' for cross-validation"
  )
  
    
  if( add && type %in% c("hist.pit", "mc") ) warning( 
    "plot will be replaced (adding elements not possible" 
  )
  
  ## extract scaled residuals and predictions standard error of signal for
  ## a robust fit
  
  #   robust    <- attr( x, "tuning.psi" ) < control.georob()[["tuning.psi.nr"]]  
  #   scld.res  <- attr( x, "scaled.residuals" )
  #   se.signal <- attr( x, "se.signal" )
  
  ## compute validation statistics

  if( type %in% c( "hist.pit", "ecdf.pit", "mc", "bs" ) ){
    
    result <- validate.predictions( 
      data = x[["data"]],
      pred = x[["pred"]],
      se.pred = x[["se"]],
      statistic = gsub( "ecdf.", "", gsub( "hist.", "", type ) ), 
    #       robust = robust,
      ncutoff = ncutoff#,
    #       se.signal = se.signal,
    #       scld.res <- scld.res
    )
  }

  if( missing( col ) ) col <- 1
  if( missing( pch ) ) pch <- 1
  if( missing( lty ) ) lty <- 1

  
    
  switch(
    type,
    sc = {
      
      ##  scatterplot of (transformed) measurements vs. predictions
      
      if( missing( main ) ) main <- "data vs. predictions"
      if( missing( xlab ) ) xlab <- "predictions"
      if( missing( ylab ) ) ylab <- "data"
      
      if( add ){
        points( data ~ pred, x, col = col, pch = pch, ... )
      } else {        
        plot( 
          data ~ pred, x, col = col, pch = pch, 
          main = main, xlab = xlab, ylab = ylab, ... 
        )
      }
      
      if( smooth ){
        lines( loess.smooth( x[["pred"]], x[["data"]], span = span ) )
      }
      
      
    },        
    lgn.sc = {
      
      ##  scatterplot of original measurements vs.  back-transformded
      ##  lognormal predictions
      
      if( missing( main ) ) main <- "data vs. back-transformed predictions"
      if( missing( xlab ) ) xlab <- "back-transformed predictions"
      if( missing( ylab ) ) ylab <- "data"
      
      if( add ){
        points( lgn.data ~ lgn.pred, x, col = col, pch = pch, ... )
      } else {        
        plot( 
          lgn.data ~ lgn.pred, x, col = col, pch = pch, 
          main = main, xlab = xlab, ylab = ylab, ... 
        )
      }
      
      if( smooth ){
        lines( loess.smooth( x[["lgn.pred"]], x[["lgn.data"]], span = span ) )
      }
      
      
    }, 
    ta = {
      
      ##  Tukey-Anscombe plot
      
      if( missing( main ) ) main <- "Tukey-Anscombe plot"
      if( missing( xlab ) ) xlab <- "predictions"
      if( missing( ylab ) ) ylab <- "standardized prediction errors"
      
      if( add ){
        points( I((data-pred)/se) ~ pred, x, col = col, pch = pch, ... )
      } else {        
        plot( 
          I((data-pred)/se) ~ pred, x, col = col, pch = pch, 
          main = main, xlab = xlab, ylab = ylab, ... 
        )
      }
      
    },
    qq = {
      
      ##  normal QQ-Plot of standardized prediction errors
      
      if( missing( main ) ) main <- "normal-QQ-plot of standardized prediction errors"
      if( missing( xlab ) ) xlab <- "quantile N(0,1)"
      if( missing( ylab ) ) ylab <- "quantiles of standardized prediction errors"
      
      r.qq <- with( x, qqnorm( ( data - pred ) / se, plot.it = FALSE ) )
      
      if( add ){
        points( r.qq, col = col, pch = pch, ... )
      } else {
        plot( r.qq, col = col, pch = pch, main = main, xlab = xlab, ylab = ylab, ... )
      }
    },
    hist.pit = {
      
      ##  histogramm of probability-integral-transformation
      
      if( missing( main ) ) main <- "histogramm PIT-values"
      if( missing( xlab ) ) xlab <- "PIT"
      if( missing( ylab ) ) ylab <- "density"
      
      r.hist <- hist( 
        result,
        col = col, lty = lty, 
        main = main, xlab = xlab, ylab = ylab, freq = FALSE, ... )
    },
    ecdf.pit = {
      
      ##  ecdf of probability-integral-transformation
      
      if( missing( main ) ) main <- "ecdf PIT-values"
      if( missing( xlab ) ) xlab <- "PIT"
      if( missing( ylab ) ) ylab <- "probability"
      
      r.hist <- plot(
        ecdf(result),  add = add,
        col = col, lty = lty, 
        main = main, xlab = xlab, ylab = ylab, ... 
      )
    },
    mc = {
      
      ##  narginal calibration plots: ecdf of measurements and mean
      ##  predictive cdf
      
      if( missing( main ) ) main <- "empirical cdf of data and mean predictive cdfs"
      if( missing( xlab ) ) xlab <- "data or predicitons"
      if( missing( ylab ) ) ylab <- "probability"
      
      matplot( 
        result[["y"]], 
        result[, c( "ghat", "fbar" )], type = "l",
        col = c( "black", "red" ),
        lty = c( "solid", "dashed" ),
        main = main, xlab = xlab, ylab = ylab,
        ...
      )
      
      t.usr <- par( "usr" )
      t.usr[3:4] <- with( result, range( fbar - ghat ) ) *1.04
      par( usr = t.usr )
      with( result, lines( y, fbar-ghat, col= "blue", lty = "dotted" ) )
      axis(2, pos = t.usr[2], col.axis = "blue", col.ticks = "blue" )
      legend( 
        "topleft",
        lty = c("solid", "dashed", "dotted" ), 
        col = c( "black", "red", "blue" ), bty = "n", cex = 1,
        legend = c(
          expression( paste( "empirical cdf ", hat(G) ) ),
          expression( paste( "mean predictive cdf ", bar(F) ) ),
          expression( bar(F)-hat(G) )
        )
      )
    },
    bs  ={
      
      # plot of brier score vs. cutoff
      
      if( missing( main ) ) main <- "Brier score vs. cutoff"
      if( missing( xlab ) ) xlab <- "cutoff"
      if( missing( ylab ) ) ylab <- "Brier score"
      
      if( add ){
        lines( result[["y"]], result[["bs"]], col = col, lty = lty, ... )
      } else {
        plot( result[["y"]], result[["bs"]], type = "l", col = col, lty = lty, 
          main = main, xlab = xlab, ylab = ylab, ...
        )
      }
    }
  )
  invisible( NULL )
}

  

##  ###########################################################################

print.cv.georob <-
  function( 
    x, digits = max(3, getOption("digits") - 3), ...
  )
{   ## print method for class "cv.georob"  
  
  ## 2011-10-13 A. Papritz
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-06-25 AP new method to compute pit, mc, bs and crps (Gaussian and robust)
  
  x <- x[["pred"]]
  
  st <- validate.predictions( 
    data = x[["data"]],
    pred = x[["pred"]],
    se.pred = x[["se"]],
    statistic = "st",
    ...
  )
  
  print(
    format( st, digits = digits ), print.gap = 2, 
    quote = FALSE
  )
  
  invisible( x )
}


##  ###########################################################################

summary.cv.georob <-
  function( object, se = FALSE, ... )
{
  
  ## summary method for class "cv.georob" 
  
  ## function computes statistics of the cross-validation errors
  
  ## 2011-10-13 A. Papritz
  ## 2012-05-21 ap
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-06-25 AP new method to compute pit, mc, bs and crps (Gaussian and robust)
  
  object <- object[["pred"]]
  
  ## extract scaled residuals and predictions standard error of signal for
  ## a robust fit
  
  #   robust    <- attr( object, "tuning.psi" ) < control.georob()[["tuning.psi.nr"]]  
  #   scld.res  <- attr( object, "scaled.residuals" )
  #   se.signal <- attr( object, "se.signal" )
  
  crps <- validate.predictions( 
    data = object[["data"]],
    pred = object[["pred"]],
    se.pred = object[["se"]],
    statistic = "crps",
    #     robust = robust,
    ncutoff = length( object[["data"]] )#,
    #     se.signal = se.signal,
    #     scld.res = scld.res
  )
  
  st <- validate.predictions( 
    data = object[["data"]],
    pred = object[["pred"]],
    se.pred = object[["se"]],
    statistic = "st"
  )
  
  if( !is.null( object[["lgn.pred"]] ) ){
    st.lgn <- validate.predictions( 
      data = object[["lgn.data"]],
      pred = object[["lgn.pred"]],
      se.pred = object[["lgn.se"]],
      statistic = "st"
    )
  } else {
    st.lgn <- NULL
  }
  
  ## collect results
  
  result <- list( st = st, crps = crps, st.lgn = st.lgn )

  ## compute standard errors of criteria across cross-validation sets
  
  criteria <- NULL
  
  if( se && !is.null( object[["subset"]] ) ){
    
    criteria <- t( sapply(
        tapply(
          1:nrow( object ),
          factor( object[["subset"]] ),
          function( 
            i, data, pred, se.pred,
            #              robust, se.signal, scld.res, 
            lgn.data, lgn.pred, lgn.se.pred 
          ){
            
            crps <- validate.predictions( 
              data = data[i],
              pred = pred[i],
              se.pred = se.pred[i],
              statistic = "crps",
              #               robust = robust,
              ncutoff = length( data[i] )#,
              #               se.signal = se.signal[i],
              #               scld.res = scld.res
            )
            
            st <- validate.predictions( 
              data = data[i],
              pred = pred[i],
              se.pred = se.pred[i],
              statistic = "st"
            )
            
            if( !is.null( lgn.pred ) ){
              st.lgn <- validate.predictions( 
                data = lgn.data[i],
                pred = lgn.pred[i],
                se.pred = lgn.se.pred[i],
                statistic = "st"
              )
              names( st.lgn ) <- paste( names( st.lgn ), "lgn", sep = "." )
            } else {
              st.lgn <- NULL
            }
            
            return( c( st, st.lgn, crps = crps ) )
            
          },
          data = object[["data"]],
          pred = object[["pred"]],
          se.pred = object[["se"]],
          #           robust = robust,
          #           se.signal = se.signal,
          #           scld.res = scld.res,
          lgn.data = object[["lgn.data"]],
          lgn.pred = object[["lgn.pred"]],
          lgn.se.pred = object[["lgn.se"]]
        ),
        function( x ) x
      ))
    
    se.criteria <- apply(
      criteria, 2,
      function( x ) sd( x ) / sqrt( length( x ) )
    )
    
    result[["se.st"]] <- se.criteria[c( "me", "mede", "rmse", "made", "qne", "msse", "medsse")]
    result[["se.crps"]] <- se.criteria["crps"]
    if( !is.null( st.lgn ) ){
      result[["se.st.lgn"]] <- se.criteria[
        c( "me.lgn", "mede.lgn", "rmse.lgn", "made.lgn", "qne.lgn", "msse.lgn", "medsse.lgn")
      ]
      names( result[["se.st.lgn"]] ) <- gsub( ".lgn", "", names( result[["se.st.lgn"]] ) )
    }
    
  }
  
  class( result ) <- "summary.cv.georob"
  
  if( !is.null( criteria ) ) attr( result, "statistics" ) <- criteria
 
  
  return( result )
  
}

##  ###########################################################################

print.summary.cv.georob <-
  function( 
    x, digits = max(3, getOption("digits") - 3), ...
  )
{
  
  ## print method for class "summary.cv.georob"  
  
  ## 2011-12-20 A. Papritz
  ## 2012-05-21 ap
  ## 2012-12-18 AP invisible(x)
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  

  result <- c( x[["st"]], crps = x[["crps"]] )
  if( !is.null( x[["se.st"]] ) ){
    result <- rbind( result, c( x[["se.st"]], crps = x[["se.crps"]] ) )
    rownames( result ) <- c( "", "se" )
  }
  
  cat( "\nStatistics of cross-validation prediction errors\n" )
  print(
    format( result, digits = digits ), print.gap = 2, 
    quote = FALSE
  )
  
  if( !is.null( x[["st.lgn"]] ) ){
    result <- x[["st.lgn"]]
    if( !is.null( x[["se.st.lgn"]] ) ){
      result <- rbind( x[["st.lgn"]], x[["se.st.lgn"]] )
      rownames( result ) <- c( "", "se" )
    }
    
    cat( "\nStatistics of back-transformed cross-validation prediction errors\n" )
    print(
      format( result, digits = digits ), print.gap = 2, 
      quote = FALSE
    )
  }
  
  invisible( x )
  
}

##  ###########################################################################

rstudent.cv.georob <-
  function( model, ... )
{
  
  ## Function extracts studentized residuals from cv.georob object
  
  ## Arguments:
  
  ## model     cv.georob object
  ## ...       further arguments (currently not used)
  
  ## 2011-10-13 A. Papritz
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  
  if( !identical( class( model )[1], "cv.georob" ) ) stop(
    "model is not of class 'cv.georob'" 
  )
  
  model <- model[["pred"]]
  
  ( model[["data"]] - model[["pred"]] ) / model[["se"]]
  
}

# ##  ###########################################################################
# 
# cv.variomodel <-
#   function( object, geodata, ... )
# {
#   
#   ## Wrapper function for cross-validation of object of class variomodel{geoR}
#   ## by function xvalid{geoR}
#   
#   ## Arguments:
#   
#   ## model     an object of class "variomodel{geoR}
#   ## ...       further arguements passed to xvalid{geoR), cf. respective help page
#   
#   ## 2012-11-22 A. Papritz
#   
#   call.fc <- match.call()
# 
#   res <- geoR::xvalid( model = object, ... )
#   
#   if( !is.null( attr( res, "geodata.xvalid" ) ) ){
#     attr( res, "geodata.xvalid" ) <- call.fc[["geodata"]]
#   }
#   if( !is.null( attr( res, "locations.xvalid" ) ) ){
#     attr( res, "locations.xvalid" ) <- call.fc[["locations.xvalid"]]
#   }
#   
#   return(res)
#   
# }
# 
# cv.likGRF <-
#   function( object, geodata, ... )
# {
#   
#   ## Wrapper function for cross-validation of object of class variomodel{geoR}
#   ## by function xvalid{geoR}
#   
#   ## Arguments:
#   
#   ## model     an object of class "likGRF{geoR}
#   ## ...       further arguements passed to xvalid{geoR), cf. respective help page
#   
#   ## 2012-11-22 A. Papritz
#   
#   call.fc <- match.call()
# 
#   res <- geoR::xvalid( model = object, geodata = geodata, ... )
#   
#   if( !is.null( attr( res, "geodata.xvalid" ) ) ){
#     attr( res, "geodata.xvalid" ) <- call.fc[["geodata"]]
#   }
#   if( !is.null( attr( res, "locations.xvalid" ) ) ){
#     attr( res, "locations.xvalid" ) <- call.fc[["locations.xvalid"]]
#   }
#   
#   return(res)
#   
# }

## ======================================================================

validate.predictions <- 
  function( 
    data,
    pred,
    se.pred,
    statistic = c( "crps", "pit", "mc", "bs", "st" ),	
    #     robust = FALSE,
    ncutoff = NULL#,
    #     se.signal = NULL,
    #     scld.res = NULL
  )
{
  
  ## function computes several statistics to validate probabilistic
  ## predictions, cf.  Gneiting et al., 2007, JRSSB
  
  ## 2011-20-21 A. Papritz
  ## 2012-05-04 AP coping with NAs
  ## 2013-06-12 AP substituting [["x"]] for $x in all lists
  ## 2015-06-25 AP new method to compute pit, mc, bs and crps (Gaussian and robust)
  ## 2015-11-27 AP checking whether mandatory arguments were provided
  
  ## checking whether mandatory arguments were provided
  
  if( missing( data ) || missing( pred ) || missing( se.pred ) ) stop(
	"some mandatory arguments are missing" 
  )
  
  statistic = match.arg( statistic )
  
  ## exclude item with NAs
  
  param <- na.exclude( cbind( data, pred, se.pred ) )
  #   param <- na.exclude( cbind( data, pred, se.pred, se.signal ) )
  #   if( robust ) scld.res <- scld.res[!is.na(scld.res)]
  
  if( missing( ncutoff ) || is.null( ncutoff ) ) ncutoff <- min( 500, NROW( param ) )
  
  result <- switch(
    statistic,
    pit = {
      
      ## probability integral transformation
      
      #       if( !robust ){
      pnorm( param[, "data"], mean = param[, "pred"], sd = param[, "se.pred"] )
      #       } else {
      #         apply( 
      #           param, 1,
      #           function( x, r ){
      #             pnorMix( x["data"], norMix( x["pred"] + r, sigma = x["se.signal"] ) )
      #           }, r = scld.res
      #         )
      #       }
      
    },
    mc = ,
    bs = {
      
      ## marginal calibration and brier score
      
      margin.calib <- data.frame(
        y = t.x <- unique( t.y <- sort( c( param[, "data"] ) ) ),    # y: cutoff
        ghat = cumsum( tabulate( match( t.y, t.x ) ) ) / length(t.y) # Ghat(y)
      )
      
      ## 'dilute' margin.calib
      
      t.sel <- trunc( 
        seq( 
          from = as.integer(1), 
          to = nrow( margin.calib ), 
          length.out = min( nrow( margin.calib ), ncutoff ) 
        ) 
      )
      margin.calib <- margin.calib[t.sel,]
      
      ## compute mean of predictive distriutions Fhat_i(y) and Brier score 
      
      t.bla <- t(
        sapply(
          margin.calib[, "y"],
          function( cutoff, param ){
            #           function( cutoff, param, scld.res ){
            
            ## compute Fhat_i(y)
            
            #             if( !robust ){
            t.p <- pnorm( cutoff, mean = param[, "pred"], sd = param[, "se.pred"] )
            #             } else {
            #               t.p <- ppd.resp.rob( cutoff, m = param[, "pred"], s = param[, "se.signal"], r = scld.res )
            #             }
            
            ## compute barFhat(y) and BS(y)
            
            c( 
              fbar = mean( t.p ), 
              bs = mean( ( t.p - as.numeric( param[, "data"] <= cutoff ) )^2 )
            )
          },
          param = param#,
          #           scld.res = scld.res
        )
      )
      cbind(
        margin.calib, as.data.frame( t.bla ) 
      )
    },
    crps = {
      
      ## continuous ranked probability score
      
      #       if( !robust ){
      mean( crpsnorm(
          y = param[, "data"], m = param[, "pred"], s = param[, "se.pred"] 
        ))  
      #       } else {
      #         mean( crpspd.resp.rob( 
      #           y = param[, "data"], m = param[, "pred"], s = param[, "se.signal"], r = scld.res
      #         ))
      #       }
      
    },
    st = {
      
      ## statistics of (standardized) prediction errors
      
      error <- param[, "data"] - param[, "pred"]
      std.error <- error / param[, "se.pred"]
      
      statistics <- c( 
        me = mean( error ),
        mede = median( error ),
        rmse = sqrt( mean( error^2 ) ),
        made = mad( error, center = 0 ),
        qne = Qn( error, finite.corr = TRUE ),
        msse = mean( std.error^2 ),
        medsse = median( std.error^2 )
      )
      
    }
  )
  
  return( result )   
  
}


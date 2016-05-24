#######################################################################
##   NORM Version 2 client functions for R                         ####
#######################################################################
## Functions to export:
##     emNorm
##     impNorm
##     loglikNorm
##     logpostNorm
##     mcmcNorm
##     miInference (source code in miInference.R)
#######################################################################
emNorm <- function(obj, ...){
   # S3 generic function
   UseMethod("emNorm")}
#######################################################################
emNorm.default <- function(obj, x=NULL, intercept=TRUE,
   iter.max=1000, criterion=NULL, estimate.worst=TRUE, 
   prior="uniform", prior.df=NULL, 
   prior.sscp=NULL, starting.values=NULL, ...){
   #########################################
   y <- obj
   if( class(y) == "data.frame" ){
      status <- logical( length(y) )
      for( j in 1:length(y) ) status[j] <- is.factor( y[[j]] )
      if( any(status) ){ 
         msg <- paste( "Factors in argument \"y\"",
               "converted to mode \"numeric\"." )
         warning( msg )
      }
   }  
   y <- data.matrix(y)
   storage.mode(y) <- "double"
   if( is.null( colnames(y) ) ){
      colnames(y) <- paste("Y", as.character( 1:ncol(y) ), sep="" )
   }
   #########################################
   # get matrix of predictors
   if( is.null(x) ){
      # default to a column of ones
      x <- matrix( 1., nrow(y), 1)
      colnames(x) <- "CONST"
      }
   else{
      # coerce x to a numeric matrix and check dimensions 
      if( class(x) == "data.frame" ){
         status <- logical( length(x) )
         for( j in 1:length(x) ) status[j] <- is.factor( x[[j]] )
         if( any(status) ){
            msg <- paste( "Factors in argument \"x\"",
               "converted to mode \"numeric\"." )
            warning( msg )
         }
      }  
      x <- data.matrix(x)
      if( nrow(x) != nrow(y) ){
         stop("Arguments x and y do not have the same number of rows.")
      }
      if( any( is.na(x) ) ){
         stop("Missing values in x are not allowed.")
      }
      if( is.null( colnames(x) )) {
      	    colnames(x) <- paste("X", as.character( 1:ncol(x) ),
	    sep="" )
      }
      if( intercept ){
         x <- cbind( CONST=1., x )
      }
   }
   rownames(x) <- rownames(y)
   storage.mode(x) <- "double"
   #########################################
   # missing-value code to pass to Fortran
   mvcode <- .Machine$double.xmax
   y[ is.na(y) ] <- mvcode
   #########################################
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 8L )
   #########################################
   # get starting values
   if( is.null(starting.values) ){
      # dimension arrays to pass to Fortran
      beta <- matrix( 0., ncol(x), ncol(y) )
      sigma <- matrix( 0., ncol(y), ncol(y) )
      startval.present <- FALSE
      }
   else{
      # get user-supplied starting values
      beta <- data.matrix( starting.values$beta )
      sigma <- data.matrix( starting.values$sigma )
      if( ( nrow(beta)!=ncol(x) ) | (ncol(beta)!=ncol(y)) )
         stop("Incorrect dimensions for argument beta.")
      if( ( nrow(sigma)!=ncol(y) ) | (ncol(sigma)!=ncol(y)) )
         stop("Incorrect dimensions for argument sigma.")
      startval.present <- TRUE
      }   
   storage.mode(beta) <- "double"
   storage.mode(sigma) <- "double"
   rownames( beta ) <- colnames( x )
   colnames( beta ) <- colnames( y )
   rownames( sigma ) <- colnames( y )
   colnames( sigma ) <- colnames( y )
   #########################################
   if( prior == "uniform" ){
      prior.type.int <- 1L
      prior.df <- -( ncol(y) + 1 ) 
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )
   }
   else if( prior == "jeffreys" ){
      prior.type.int <- 2L
      prior.df <- 0.
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )
   }
   else if( prior == "ridge" ){
      prior.type.int <- 3L
      if( is.null( prior.df ) )
         stop("Argument prior.df must be provided when prior = \"ridge\".")
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )
   }
   else if( prior == "invwish" ){
      prior.type.int <- 4L
      if( is.null( prior.df ) )
         stop("Argument prior.df must be provided when prior = \"invwish\".")
      length( prior.df ) <- 1
      if( is.null( prior.sscp ) )
         stop("Argument prior.sscp must be provided when prior = \"invwish\".")
      prior.sscp <- data.matrix( prior.sscp )
      if( ( nrow( prior.sscp )!=ncol(y) ) | ( ncol( prior.sscp )!=ncol(y) ) )
         stop("Argument prior.sscp has incorrect dimensions.")
   }
   else{
      msg <- paste( "Prior type \"", prior, "\" not recognized.", sep="")
      stop( msg )
   }
   storage.mode( prior.type.int ) <- "integer"
   length( prior.type.int ) <- 1L
   storage.mode( prior.df ) <- "double"
   length( prior.df ) <- 1L
   storage.mode( prior.sscp ) <- "double"
   #########################################
   # check remaining user-supplied arguments
   if( !missing(iter.max) ){
      length(iter.max) <- 1L
      }
   storage.mode(iter.max) <- "integer"
   if( iter.max < 0 )
      stop("Argument iter.max cannot be negative.")
   #
   if( !is.null( criterion ) ){
      length(criterion) <- 1L
      storage.mode( criterion ) <- "double"
      if( criterion < 0. )
         stop("Argument criterion cannot be negative.")
   }
   else{
      criterion <- as.double( 1.e-5 )
   }
   if( is.null( estimate.worst ) ) estimate.worst <- FALSE
   storage.mode( estimate.worst ) <- "logical"
   length(estimate.worst) <- 1L
   #########################################
   nparam <- as.integer( ncol(x) * ncol(y) + 
       ( ncol(y) * (ncol(y) + 1L) ) / 2L )
   #########################################
   tmp <- .Fortran("norm_em",
      n = nrow(y),
      r = ncol(y),
      p = ncol(x),
      x = x,
      y = y,
      mvcode = mvcode,
      prior.type.int = prior.type.int,
      prior.df = prior.df,
      prior.sscp = prior.sscp,
      iter.max = iter.max,
      criterion = criterion,
      estimate.worst = estimate.worst,
      startval.present = startval.present,
      iter = integer(1),
      converged = logical(1),
      rel.diff = numeric(1),
      loglik = numeric(iter.max),
      logpost = numeric(iter.max),
      beta = beta,
      sigma = sigma,
      y.imp = y,
      npatt = integer(1),
      mis = matrix( logical(1), nrow(y), ncol(y) ),
      n.in.patt = integer( nrow(y) ),
      n.obs = integer( ncol(y) ),
      which.patt = integer( nrow(y) ),
      ybar = numeric( ncol(y) ), 
      ysdv = numeric( ncol(y) ), 
      rate.beta = beta,
      rate.sigma = sigma,
      em.worst.ok = logical(1),
      worst.frac = numeric(1),
      nparam = nparam,
      worst.linear.coef = numeric(nparam),
      status = integer(1),
      msg.len.max = msg.len.max,
      msg.codes = msg.codes,
      msg.len.actual = integer(1),
      PACKAGE="norm2"
      )
   #########################################
   # display message from Fortran
   msg.lines <- msgNorm( tmp$msg.codes, tmp$msg.len.actual )
   if( is.null( msg.lines ) ){
      msg <- "OK"
      }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
      }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #########################################
   if( tmp$status != 0 ){
      result <- invisible( NULL )
   }
   else{
      if( (!tmp$converged) & (iter.max>0) ){
         warning( paste( "Algorithm did not converge by iteration", 
            format(tmp$iter) ) ) 
      }
      #########################################
      miss.patt <- tmp$mis[ 1:tmp$npatt, , drop=FALSE]
      colnames(miss.patt) <- colnames(y)
      miss.patt.freq <- tmp$n.in.patt[ 1:tmp$npatt ]
      n.obs <- tmp$n.obs
      names( n.obs ) <- colnames( y )
      which.patt <- tmp$which.patt
      names( which.patt ) <- rownames( y )
      ybar <- tmp$ybar
      names( ybar ) <- colnames( y )
      ysdv <- tmp$ysdv
      names( ysdv ) <- colnames( y )
      rel.diff <- tmp$rel.diff
      if( tmp$iter == 0 ) rel.diff <- NA
      #########################################
      y <- tmp$y
      y[ y == tmp$mvcode ] <- NA
      #########################################
      if( startval.present ){
      	 starting.values <- list( beta=beta, sigma=sigma )
      }
      else{
         starting.values <- NULL
      }
      #########################################
      if( tmp$iter >= 1 ){
         loglik <- tmp$loglik[1:tmp$iter]
         logpost <- tmp$logpost[1:tmp$iter]
      }
      else{
         loglik <- NULL
	 logpost <- NULL
      }
      #########################################
      if( prior == "ridge" ){
         prior.sscp <- NULL
      }
      else if( prior == "uniform" ){
         prior.df <- NULL
         prior.sscp <- NULL
      }
      else if( prior == "jeffreys" ){
         prior.df <- NULL
         prior.sscp <- NULL
      }
      #########################################
      em.worst.ok <- tmp$em.worst.ok
      if( em.worst.ok ){
         worst.frac <- tmp$worst.frac
         if( worst.frac == 0. ){
            worst.linear.coef <- NULL
            }
         else{
            worst.linear.coef <- tmp$worst.linear.coef
            }
         }
      else{
         worst.frac <- NULL
         worst.linear.coef <- NULL
         }
      #########################################
      result <- list(
         y = y,
         x = x,
         method = "EM",
         prior = prior,
         prior.df = prior.df,
         prior.sscp = prior.sscp, 
	 starting.values = starting.values,
         iter = tmp$iter,
         converged = tmp$converged,
         criterion = criterion,
	 estimate.worst = estimate.worst,
         loglik  = loglik,
         logpost = logpost,
         param = list(beta=tmp$beta, sigma=tmp$sigma),
	 param.rate = list(beta=tmp$rate.beta, sigma=tmp$rate.sigma),
         y.mean.imp = tmp$y.imp,
         miss.patt = miss.patt,
         miss.patt.freq = miss.patt.freq,
         n.obs = n.obs,
         which.patt = which.patt,
         rel.diff = rel.diff,
         ybar = ybar,
         ysdv = ysdv,
         em.worst.ok = em.worst.ok,
	 worst.frac = worst.frac,
	 worst.linear.coef = worst.linear.coef,
         msg = msg)
      # set class and return
      class(result) <- "norm"
   }
   return(result)}
#######################################################################
emNorm.norm <- function( obj, iter.max = 1000, 
   criterion = obj$criterion, estimate.worst = obj$estimate.worst,
   prior = obj$prior, prior.df = obj$prior.df,
   prior.sscp = obj$prior.sscp,
   starting.values = obj$param, ...){
   ###############################################
   # S3 method for class "norm"
   ###############################################
   result <- emNorm.default( obj$y, x = obj$x,
      intercept = FALSE, iter.max = iter.max, criterion = criterion, 
      estimate.worst = estimate.worst, starting.values = starting.values, 
      prior = prior, prior.df = prior.df, prior.sscp = prior.sscp,
      ... )
   return(result)}
#######################################################################
emNorm.formula <- function( formula, data, iter.max=1000,
   criterion=NULL, estimate.worst=TRUE, 
   prior="uniform", prior.df=NULL, prior.sscp=NULL,
   starting.values=NULL, ...){
   ###############################################
   # S3 method for formula object
   ###############################################
   # Prepare the model frame
   mf <- match.call(expand.dots=FALSE)
   mf[[1]] <- as.name("model.frame")
   mf$iter.max <- mf$criterion <- mf$prior <- mf$prior.df <-
       mf$prior.sscp <- mf$starting.values <- NULL
   m <- match( c("formula","data"), names(mf), 0L)
   mf <- mf[c(1L,m)]
   mf$na.action <- as.name("na.pass")
   mf <- eval( mf, parent.frame() )
   ###############################################
   # Extract the matrix of response variables
   y <- model.response( mf )
   y <- data.matrix( y )
   storage.mode(y) <- "double"
   if( ncol(y) == 1 ) colnames(y) <- all.vars(formula)[1]
   ###############################################
   # Extract the matrix of predictors
   x <- model.matrix( formula, mf )
   if( any( is.na(x) ) ){
      stop( "Missing values in predictors not allowed." )
      }
   ###############################################
   result <- emNorm.default( y, x = x, intercept = FALSE,
      iter.max = iter.max, criterion = criterion,
      estimate.worst = estimate.worst,
      starting.values = starting.values, 
      prior = prior, prior.df = prior.df, prior.sscp = prior.sscp,
      ... )
   return(result)}
#######################################################################
#######################################################################
#######################################################################
#######################################################################
mcmcNorm <- function(obj, ...){
   # S3 generic function
   UseMethod("mcmcNorm")}
#######################################################################
mcmcNorm.default <- function(obj, x=NULL, intercept=TRUE, 
   starting.values,
   iter=1000, multicycle=NULL, seeds=NULL, prior="uniform", 
   prior.df=NULL, prior.sscp=NULL,
   save.all.series=TRUE,
   save.worst.series=FALSE, worst.linear.coef=NULL,
   impute.every=NULL, ...){ 
   #########################################
   y <- obj
   if( class(y) == "data.frame" ){
      status <- logical( length(y) )
      for( j in 1:length(y) ) status[j] <- is.factor( y[[j]] )
      if( any(status) ){ 
         msg <- paste( "Factors in argument \"y\"",
               "converted to mode \"numeric\"." )
         warning( msg )
      }
   }  
   y <- data.matrix(y)
   storage.mode(y) <- "double"
   if( is.null( colnames(y) ) ){
      colnames(y) <- paste("Y", as.character( 1:ncol(y) ), sep="" )
   }
   #########################################
   # get matrix of predictors
   if( is.null(x) ){
      # default to a column of ones
      x <- matrix( 1., nrow(y), 1)
      colnames(x) <- "CONST"
      }
   else{
      # coerce x to a numeric matrix and check dimensions 
      if( class(x) == "data.frame" ){
         status <- logical( length(x) )
         for( j in 1:length(x) ) status[j] <- is.factor( x[[j]] )
         if( any(status) ){ 
         msg <- paste( "Factors in argument \"x\"",
               "converted to mode \"numeric\"." )
         warning( msg )
         }
      }  
      x <- data.matrix(x)
      if( nrow(x) != nrow(y) ){
         stop("Arguments x and y do not have the same number of rows.")
      }
      if( any( is.na(x) ) ){
         stop("Missing values in x are not allowed.")
      }
      if( is.null( colnames(x) )) {
      	    colnames(x) <- paste("X", as.character( 1:ncol(x) ),
	    sep="" )
      }
      if( intercept ){
         x <- cbind( CONST=1., x )
      }
   }
   rownames(x) <- rownames(y)
   storage.mode(x) <- "double"
   #########################################
   # missing-value code to pass to Fortran
   mvcode <- .Machine$double.xmax
   y[ is.na(y) ] <- mvcode
   #########################################
   # create a blank string of length 255 to 
   # pass to Fortran for messaging
   msg.len <- as.integer(255)
   msg <- format("", width=msg.len)
   #########################################
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 8L )
   #########################################
   # check user-supplied starting values
   beta <- data.matrix( starting.values$beta )
   sigma <- data.matrix( starting.values$sigma )
   if( ( nrow(beta)!=ncol(x) ) | (ncol(beta)!=ncol(y)) )
      stop("Incorrect dimensions for argument beta.")
   if( ( nrow(sigma)!=ncol(y) ) | (ncol(sigma)!=ncol(y)) )
      stop("Incorrect dimensions for argument sigma.")
   storage.mode(beta) <- "double"
   storage.mode(sigma) <- "double"
   rownames( beta ) <- colnames( x )
   colnames( beta ) <- colnames( y )
   rownames( sigma ) <- colnames( y )
   colnames( sigma ) <- colnames( y )
   #########################################
   if( prior == "uniform" ){
      prior.type.int <- 1L
      prior.df <- -( ncol(y) + 1. ) 
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )

   }
   else if( prior == "jeffreys" ){
      prior.type.int <- 2L
      prior.df <- 0.
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )
   }
   else if( prior == "ridge" ){
      prior.type.int <- 3L
      if( is.null( prior.df ) )
         stop("Argument prior.df must be provided when prior = \"ridge\".")
      length( prior.df ) <- 1L
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )
   }
   else if( prior == "invwish" ){
      prior.type.int <- 4L
      if( is.null( prior.df ) )
         stop("Argument prior.df must be provided when prior = \"invwish\".")
      length( prior.df ) <- 1L
      if( is.null( prior.sscp ) )
         stop("Argument prior.sscp must be provided when prior = \"invwish\".")
      prior.sscp <- data.matrix( prior.sscp )
      if( ( nrow( prior.sscp )!=ncol(y) ) | ( ncol( prior.sscp )!=ncol(y) ) )
         stop("Argument prior.sscp has incorrect dimensions.")
   }
   else{
      msg <- paste( "Prior type \"", prior, "\" not recognized.", sep="")
      stop( msg )
   }
   storage.mode( prior.type.int ) <- "integer"
   length( prior.type.int ) <- 1L
   storage.mode( prior.df ) <- "double"
   length( prior.df ) <- 1L
   storage.mode( prior.sscp ) <- "double"
   #########################################
   # check remaining user-supplied arguments
   if( !missing(iter) ){
      length(iter) <- 1L
      }
   storage.mode(iter) <- "integer"
   if( iter < 1L )
      stop("Argument iter must be positive.")
   #
   if( !is.null(multicycle) ){
      length(multicycle) <- 1L
      storage.mode(multicycle) <- "integer"
   }
   else{
      multicycle <- 1L
   }
   if( multicycle < 1 )
      stop("Argument multicycle must be positive.")
   #
   if( !is.null(seeds) ){ 
      if( length(seeds) != 2L ) 
         stop("Two integer seeds were expected.")}
   else{
      seeds=ceiling( runif(2)*.Machine$integer.max ) }
   storage.mode(seeds) <- "integer"
   #
   if( is.null(impute.every) ){
      impute.every <- as.integer(0)
      }
   else{
      impute.every <- as.integer(impute.every)
      length( impute.every ) <- 1L
      if( impute.every < 0L ) 
         stop("Argument impute.every cannot be negative.")
      if( impute.every > iter ) 
         stop("Argument impute.every cannot exceed iter.")
      }
   if( impute.every == 0 ){
       nimps <- as.integer(0)
   }
   else{
      nimps <- as.integer( floor( iter / impute.every ) )
   }
   # 
   if( mode(save.all.series) != "logical" ) 
      stop("Argument save.all.series should be of mode \"logical\".")
   length( save.all.series ) <- 1L
   #########################################
   if( mode(save.worst.series) != "logical" ) 
      stop("Argument save.worst.series should be of mode \"logical\".")
   length( save.worst.series ) <- 1L
   if( save.worst.series ){
      if( is.null( worst.linear.coef ) ){
         stop("If save.worst.series, worst.linear.coef must be provided.")
         }
      }
   if( is.null( worst.linear.coef ) ) worst.linear.coef <- 
      rep(0.,  ncol(x)*ncol(y) + ( ncol(y)*(ncol(y)+1L) ) / 2L )
   storage.mode(worst.linear.coef) <- "double"
   #########################################
   # dimension other arrays
   y.imp <- y
   if( save.all.series ){
      series.length <- as.integer( iter )
      series.beta <- array(0., c( ncol(x), ncol(y), iter) )
      dimnames(series.beta) <- 
         list( rownames(beta), colnames(beta), NULL )
      series.sigma <- array(0., c( ncol(y), ncol(y), iter) )
      dimnames(series.sigma) <-
         list( rownames(sigma), colnames(sigma), NULL )
   }
   else{
      series.length <- as.integer(0)
      series.beta <- array(0., c( ncol(x), ncol(y),  0) )
      series.sigma <- array(0., c( ncol(y), ncol(y), 0) )
   }
   storage.mode(series.beta) <- "double"
   storage.mode(series.sigma) <- "double"
   imp.list <- array(0., c( nrow(y), ncol(y), nimps ) )
   storage.mode(imp.list) <- "double"
   series.worst <- numeric(iter)
   storage.mode( series.worst ) <- "double"
   #########################################
   tmp <- .Fortran("norm_mcmc",
      n = nrow(y),
      r = ncol(y),
      p = ncol(x),
      x = x,
      y = y,
      mvcode = mvcode,
      prior.type.int = prior.type.int,
      prior.df = prior.df,
      prior.sscp = prior.sscp,
      iter = iter,
      multicycle = multicycle,
      seeds = seeds,
      impute.every = impute.every,
      nimps = nimps,
      save.all.series = save.all.series,
      save.worst.series = save.worst.series,
      worst.linear.coef = worst.linear.coef,
      series.length = series.length,
      beta = beta,
      sigma = sigma,
      y.imp = y.imp,
      iter.actual = integer(1), 
      series.beta = series.beta,
      series.sigma = series.sigma,
      loglik = numeric(iter),     
      logpost = numeric(iter),
      series.worst = series.worst,
      npatt = integer(1),
      mis = matrix( logical(1), nrow(y), ncol(y) ),
      n.in.patt = integer( nrow(y) ),
      n.obs = integer( ncol(y) ),
      which.patt = integer( nrow(y) ),
      ybar = numeric( ncol(y) ), 
      ysdv = numeric( ncol(y) ), 
      imp.list = imp.list,
      nimps.actual = integer(1), 
      status = integer(1),
      msg.len.max = msg.len.max,
      msg.codes = msg.codes,
      msg.len.actual = integer(1),
      PACKAGE="norm2"
      )
   #########################################
   # display message from Fortran, if present
   msg.lines <- msgNorm( tmp$msg.codes, tmp$msg.len.actual )
   if( is.null( msg.lines ) ){
      msg <- "OK"
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( paste("Note: ", msg, sep="") )
   #########################################
   if( tmp$status != 0 ){
      result <- invisible(NULL)
   }
   else{
      ###########################################
      iter <- tmp$iter.actual
      ###########################################
      if( iter >= 1 ){
         loglik <- tmp$loglik[1:iter]
         logpost <- tmp$logpost[1:iter]
      }
      else{
         loglik <- NULL
	 logpost <- NULL
      }
      ###########################################
      nimps <- tmp$nimps.actual
      if( ( impute.every == 0 ) | ( nimps == 0 ) ){
         imp.list <- NULL
         impute.every <- NULL
      } 
      else{
         imp.list <- as.list(1:nimps)
         for( i in 1:nimps ){
            imp.list[[i]] <- tmp$imp.list[,,i]
            dimnames( imp.list[[i]] ) <- list( rownames(y),colnames(y) )
         }
      }
      #########################################
      y <- tmp$y
      y[ y == tmp$mvcode ] <- NA
      #########################################
      if( save.all.series & (iter>0) ){
          # extract beta series and reshape into matrix
          beta.array  <- tmp$series.beta[,,1:iter, drop=FALSE]
          series.beta <- matrix( numeric(1), iter, ncol(x)*ncol(y) )
          beta.names <- character( ncol(x)*ncol(y) )
	  posn <- 0
	  for( j in 1:ncol(x) ){
             for( k in 1:ncol(y) ){
	     posn <- posn + 1
	     series.beta[,posn] <- beta.array[j,k,]
             beta.names[posn] <- paste( 
                colnames(x)[j], ",", colnames(y)[k], sep="")
             }
          }
	  colnames(series.beta) <- beta.names
	  series.beta <- as.ts( series.beta )
	  # do same for sigma series
          sigma.array  <- tmp$series.sigma[,,1:iter, drop=FALSE]
          series.sigma <- matrix( numeric(1), iter, ncol(y)*(ncol(y)+1)/2)
          sigma.names  <- character( ncol(y)*(ncol(y)+1)/2 )
	  posn <- 0
	  for( k in 1:ncol(y) ){
             for( j in k:ncol(y) ){
	     posn <- posn + 1
	     series.sigma[,posn] <- sigma.array[j,k,]
             sigma.names[posn] <- paste( 
                colnames(y)[j], ",", colnames(y)[k], sep="")
             }
          }
	  colnames(series.sigma) <- sigma.names
	  series.sigma <- as.ts( series.sigma )
      }
      else{ 
         series.beta <- NULL
         series.sigma <- NULL
      }
      #########################################
      if( save.worst.series & (iter > 0) ){
         series.worst <- as.ts( tmp$series.worst[1:iter] )
         }
      else{
         series.worst <- NULL
         worst.linear.coef <- NULL
         }
      #########################################
      starting.values <- list( beta=beta, sigma=sigma )
      #########################################
      if( prior == "ridge" ){
         prior.sscp <- NULL
      }
      else if( prior == "uniform" ){
         prior.df <- NULL
         prior.sscp <- NULL
      }
      else if( prior == "jeffreys" ){
         prior.df <- NULL
         prior.sscp <- NULL
      }
      #########################################
      miss.patt <- tmp$mis[ 1:tmp$npatt, , drop=FALSE]
      colnames(miss.patt) <- colnames(y)
      miss.patt.freq <- tmp$n.in.patt[ 1:tmp$npatt ]
      n.obs <- tmp$n.obs
      names( n.obs ) <- colnames( y )
      which.patt <- tmp$which.patt
      names( which.patt ) <- rownames( y )
      ybar <- tmp$ybar
      names( ybar ) <- colnames( y )
      ysdv <- tmp$ysdv
      names( ysdv ) <- colnames( y )
      #########################################
      result <- list(
         y = y,
         x = x,
         method = "MCMC",
         prior = prior,
         prior.df = prior.df,
         prior.sscp = prior.sscp, 
         iter = iter, 
	 multicycle = multicycle, 
         seeds = seeds, 
	 starting.values = starting.values,
         param = list(beta=tmp$beta, sigma=tmp$sigma),
	 loglik = loglik,
	 logpost = logpost,
         series.worst = series.worst,
	 series.beta = series.beta,
	 series.sigma = series.sigma,
         y.imp = tmp$y.imp,
         impute.every = impute.every,
         imp.list = imp.list,
         miss.patt = miss.patt,
         miss.patt.freq = miss.patt.freq,
         n.obs = n.obs,
         which.patt = which.patt,
         worst.linear.coef = worst.linear.coef,
         ybar = ybar,
         ysdv = ysdv,
         msg = msg)
      # set class and return
      class(result) <- "norm"
   }
   return(result)}
#######################################################################
mcmcNorm.norm <- function( obj, starting.values = obj$param,
   iter = 1000, multicycle = obj$multicycle, 
   seeds = NULL, prior = obj$prior, prior.df = obj$prior.df, 
   prior.sscp = obj$prior.sscp,
   save.all.series = !(obj$method=="MCMC" & is.null( obj$series.beta )), 
   save.worst.series = !is.null( obj$worst.linear.coef ),
   worst.linear.coef = obj$worst.linear.coef,
   impute.every = obj$impute.every, ... ){
   #########################################
   # S3 method for object of class "norm"
   #########################################
   result <- mcmcNorm.default( obj$y, x = obj$x,
      intercept = FALSE, 
      starting.values = starting.values,
      iter = iter, multicycle = multicycle,
      seeds = seeds, prior = prior, prior.df = prior.df, 
      prior.sscp = prior.sscp,
      save.all.series = save.all.series, 
      save.worst.series = save.worst.series,
      worst.linear.coef = worst.linear.coef,
      impute.every = impute.every, ... )
   return(result)}
#######################################################################
mcmcNorm.formula <- function( formula, data, starting.values, 
   iter=1000, multicycle=NULL, seeds=NULL, prior="uniform", 
   prior.df=NULL, prior.sscp=NULL, save.all.series=TRUE, 
   save.worst.series=FALSE, worst.linear.coef=NULL,
   impute.every=NULL, ...){ 
   ###############################################
   # S3 method for formula object
   ###############################################
   # Prepare the model frame
   mf <- match.call(expand.dots=FALSE)
   mf[[1]] <- as.name("model.frame")
   mf$starting.values <- mf$iter <-
       mf$multicycle <- mf$seeds <- mf$prior <- 
       mf$prior.df <- mf$prior.sscp <- 
       mf$save.all.series <- mf$impute.every <- NULL
   m <- match( c("formula","data"), names(mf), 0L)
   mf <- mf[c(1L,m)]
   mf$na.action <- as.name("na.pass")
   mf <- eval( mf, parent.frame() )
   ###############################################
   # Extract the matrix of response variables
   y <- model.response( mf )
   y <- data.matrix( y )
   storage.mode(y) <- "double"
   if( ncol(y) == 1 ) colnames(y) <- all.vars(formula)[1]
   ###############################################
   # Extract the matrix of predictors
   x <- model.matrix( formula, mf )
   if( any( is.na(x) ) ){
      stop( "Missing values in predictors not allowed." )
      }
   ###############################################
   result <- mcmcNorm.default( y, x = x, intercept = FALSE,
      starting.values = starting.values, 
      iter = iter, multicycle = multicycle, seeds = seeds,
      prior = prior, prior.df = prior.df, prior.sscp = prior.sscp,
      save.all.series = save.all.series, 
      save.worst.series = save.worst.series,
      worst.linear.coef = worst.linear.coef,
      impute.every = impute.every, ... )
   return(result)}
#######################################################################
#######################################################################
#######################################################################
#######################################################################
impNorm <- function(obj, ...){
   # S3 generic function
   UseMethod("impNorm")}
#######################################################################
impNorm.default <- function(obj, x=NULL, intercept=TRUE, param,
   seeds=NULL, method="random", ... ){
   #########################################
   y <- obj
   if( class(y) == "data.frame" ){
      status <- logical( length(y) )
      for( j in 1:length(y) ) status[j] <- is.factor( y[[j]] )
      if( any(status) ){ 
         msg <- paste( "Factors in argument \"y\"",
               "converted to mode \"numeric\"." )
         warning( msg )
      }
   }  
   y <- data.matrix(y)
   storage.mode(y) <- "double"
   if( is.null( colnames(y) ) ){
      colnames(y) <- paste("Y", as.character( 1:ncol(y) ), sep="" )
   }
   #########################################
   # get matrix of predictors
   if( is.null(x) ){
      # default to a column of ones
      x <- matrix( 1., nrow(y), 1)
      colnames(x) <- "CONST"
      }
   else{
      # check dimensions of x and coerce to a numeric matrix
      if( class(x) == "data.frame" ){
         status <- logical( length(x) )
         for( j in 1:length(x) ) status[j] <- is.factor( x[[j]] )
         if( any(status) ){ 
         msg <- paste( "Factors in argument \"x\"",
               "converted to mode \"numeric\"." )
         warning( msg )
         }
      }  
      x <- data.matrix(x)
      if( nrow(x) != nrow(y) ){
         stop("Arguments x and y do not have the same number of rows.")
      }
      if( any( is.na(x) ) ){
         stop("Missing values in x are not allowed.")
      }
      if( is.null( colnames(x) )) {
      	    colnames(x) <- paste("X", as.character( 1:ncol(x) ),
	    sep="" )
      }
      if( intercept ){
         x <- cbind( CONST=1., x )
      }
   }
   rownames(x) <- rownames(y)
   storage.mode(x) <- "double"
   #########################################
   # missing-value code to pass to Fortran
   mvcode <- .Machine$double.xmax
   y[ is.na(y) ] <- mvcode
   #########################################
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 8L )
   #########################################
   # check user-supplied parameter values
   beta <- data.matrix( param$beta )
   sigma <- data.matrix( param$sigma )
   if( ( nrow(beta)!=ncol(x) ) | (ncol(beta)!=ncol(y)) )
      stop("Incorrect dimensions for argument beta.")
   if( ( nrow(sigma)!=ncol(y) ) | (ncol(sigma)!=ncol(y)) )
      stop("Incorrect dimensions for argument sigma.")
   storage.mode(beta) <- "double"
   storage.mode(sigma) <- "double"
   rownames( beta ) <- colnames( x )
   colnames( beta ) <- colnames( y )
   rownames( sigma ) <- colnames( y )
   colnames( sigma ) <- colnames( y )
   #########################################
   if( method == "random" ){
      # check seeds
      if( !is.null(seeds) ){ 
      if( length(seeds) != 2 ) 
         stop("Two integer seeds were expected.")
      }
      else{
         seeds=ceiling( runif(2)*.Machine$integer.max )
      }
      storage.mode(seeds) <- "integer"
      tmp <- .Fortran("norm_imp_rand",
         n = nrow(y),
         r = ncol(y),
         p = ncol(x),
         x = x,
         y = y,
         mvcode = mvcode,
         seeds = seeds,
         beta = beta,
         sigma = sigma,
         y.imp = y,
         status = integer(1),
         msg.len.max = msg.len.max,
         msg.codes = msg.codes,
         msg.len.actual = integer(1),
         PACKAGE="norm2"
         )
      }
   else if( method == "predict" ){
      storage.mode(seeds) <- "integer"
      tmp <- .Fortran("norm_imp_mean",
         n = nrow(y),
         r = ncol(y),
         p = ncol(x),
         x = x,
         y = y,
         mvcode = mvcode,
         beta = beta,
         sigma = sigma,
         y.imp = y,
         status = integer(1),
         msg.len.max = msg.len.max,
         msg.codes = msg.codes,
         msg.len.actual = integer(1),
         PACKAGE="norm2"
         )
      }
   else{
      stop( paste( "Method \"", method, "\" not recognized.", sep="") )
   }
   #########################################
   # display message from Fortran
   msg.lines <- msgNorm( tmp$msg.codes, tmp$msg.len.actual )
   if( is.null( msg.lines ) ){
      msg <- "OK"
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat( msg )
   #########################################
   if( tmp$status != 0 ){
      result <- invisible(NULL)
   }
   else{
      result <- tmp$y.imp
   }
   return(result)}
#######################################################################
impNorm.norm <- function( obj, param = obj$param,
   seeds = NULL, method="random", ...){
   ##################################
   # S3 method for class "norm"
   ##################################
   result <- impNorm.default( obj$y, x = obj$x,
      intercept = FALSE, param = param, seeds = seeds,
      method = method, ... )
   return(result)}
#######################################################################
impNorm.formula <- function( formula, data, param, 
   seeds=NULL, method="random", ... ){
   ###############################################
   # S3 method for formula object
   ###############################################
   # Prepare the model frame
   mf <- match.call(expand.dots=FALSE)
   mf[[1]] <- as.name("model.frame")
   mf$param <- mf$seeds <- mf$method <- NULL
   m <- match( c("formula","data"), names(mf), 0L)
   mf <- mf[c(1L,m)]
   mf$na.action <- as.name("na.pass")
   mf <- eval( mf, parent.frame() )
   ###############################################
   # Extract the matrix of response variables
   y <- model.response( mf )
   y <- data.matrix( y )
   storage.mode(y) <- "double"
   if( ncol(y) == 1 ) colnames(y) <- all.vars(formula)[1]
   ###############################################
   # Extract the matrix of predictors
   x <- model.matrix( formula, mf )
   if( any( is.na(x) ) ){
      stop( "Missing values in predictors not allowed." )
      }
   ###############################################
   result <- impNorm.default( y, x = x,
      intercept = FALSE, param = param, 
      seeds = seeds, method = method, ... )
   return(result)}
#######################################################################
#######################################################################
#######################################################################
#######################################################################
logpostNorm <- function( obj, ...){
   # S3 generic function
   UseMethod("logpostNorm")}
#######################################################################
logpostNorm.default <- function(obj, x=NULL, intercept=TRUE, param,
   prior="uniform", prior.df=NULL, prior.sscp=NULL, ...){
   ##################################
   y <- obj
   if( class(y) == "data.frame" ){
      status <- logical( length(y) )
      for( j in 1:length(y) ) status[j] <- is.factor( y[[j]] )
      if( any(status) ){ 
         msg <- paste( "Factors in argument \"y\"",
               "converted to mode \"numeric\"." )
         warning( msg )
      }
   }  
   y <- data.matrix(y)
   storage.mode(y) <- "double"
   if( is.null( colnames(y) ) ){
      colnames(y) <- paste("Y", as.character( 1:ncol(y) ), sep="" )
   }
   ##################################
   if( mode(y) != "numeric" )
      stop("First argument should be of mode \"numeric\".")
   storage.mode(y) <- "double"
   if( is.null( colnames(y) ) ){
      colnames(y) <- paste("Y", as.character( 1:ncol(y) ), sep="" )
   }
   #########################################
   # get matrix of predictors
   if( is.null(x) ){
      # default to a column of ones
      x <- matrix( 1., nrow(y), 1)
      colnames(x) <- "CONST"
      }
   else{
      # check dimensions of x and coerce to a numeric matrix
      if( class(x) == "data.frame" ){
         status <- logical( length(x) )
         for( j in 1:length(x) ) status[j] <- is.factor( x[[j]] )
         if( any(status) ){ 
         msg <- paste( "Factors in argument \"x\"",
               "converted to mode \"numeric\"." )
         warning( msg )
         }
      }  
      x <- data.matrix(x)
      if( nrow(x) != nrow(y) ){
         stop("Arguments x and y do not have the same number of rows.")
      }
      if( any( is.na(x) ) ){
         stop("Missing values in x are not allowed.")
      }
      if( is.null( colnames(x) )) {
      	    colnames(x) <- paste("X", as.character( 1:ncol(x) ),
	    sep="" )
      }
      if( intercept ){
         x <- cbind( CONST=1., x )
      }
   }
   rownames(x) <- rownames(y)
   storage.mode(x) <- "double"
   #########################################
   # missing-value code to pass to Fortran
   mvcode <- .Machine$double.xmax
   y[ is.na(y) ] <- mvcode
   #########################################
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 8L )
   #########################################
   # check user-supplied parameter values
   beta <- data.matrix( param$beta )
   sigma <- data.matrix( param$sigma )
   if( ( nrow(beta)!=ncol(x) ) | (ncol(beta)!=ncol(y)) )
      stop("Incorrect dimensions for argument beta.")
   if( ( nrow(sigma)!=ncol(y) ) | (ncol(sigma)!=ncol(y)) )
      stop("Incorrect dimensions for argument sigma.")
   storage.mode(beta) <- "double"
   storage.mode(sigma) <- "double"
   rownames( beta ) <- colnames( x )
   colnames( beta ) <- colnames( y )
   rownames( sigma ) <- colnames( y )
   colnames( sigma ) <- colnames( y )
   #########################################
   if( prior == "uniform" ){
      prior.type.int <- 1L
      prior.df <- -( ncol(y) + 1 ) 
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )

   }
   else if( prior == "jeffreys" ){
      prior.type.int <- 2L
      prior.df <- 0.
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )
   }
   else if( prior == "ridge" ){
      prior.type.int <- 3L
      if( is.null( prior.df  ) )
         stop("Argument prior.df must be provided when prior = \"ridge\".")
      length( prior.df ) <- 1
      prior.sscp <- matrix( 0., ncol(y), ncol(y) )
   }
   else if( prior == "invwish" ){
      prior.type.int <- 4L
      if( is.null( prior.df ) )
         stop("Argument prior.df must be provided when prior = \"invwish\".")
      length( prior.df ) <- 1
      if( is.null( prior.sscp ) )
         stop("Argument prior.sscp must be provided when prior = \"invwish\".")
      prior.sscp <- data.matrix( prior.sscp )
      if( ( nrow( prior.sscp )!=ncol(y) ) | ( ncol( prior.sscp )!=ncol(y) ) )
         stop("Argument prior.sscp has incorrect dimensions.")
   }
   else{
      msg <- paste( "Prior type \"", prior, "\" not recognized.", sep="")
      stop( msg )
   }
   storage.mode( prior.type.int ) <- "integer"
   length( prior.type.int ) <- 1
   storage.mode( prior.df ) <- "double"
   storage.mode( prior.sscp ) <- "double"
   #########################################
   tmp <- .Fortran("norm_logpost",
         n = nrow(y),
         r = ncol(y),
         p = ncol(x),
         x = x,
         y = y,
         mvcode = mvcode,
         beta = beta,
         sigma = sigma,
         prior.type.int = prior.type.int,
         prior.df = prior.df,
         prior.sscp = prior.sscp,
	 logpost = numeric(1), 
         status = integer(1),
	 msg.len.max = msg.len.max,
         msg.codes = msg.codes,
         msg.len.actual = integer(1),
         PACKAGE="norm2"
         )
   #########################################
   # display message from Fortran
   msg.lines <- msgNorm( tmp$msg.codes, tmp$msg.len.actual )
   if( is.null( msg.lines ) ){
      msg <- "OK"
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat(msg)
   #########################################
   if( tmp$status != 0 ){
      result <- invisible(NULL)
   }
   else{
      result <- tmp$logpost
   }
   return(result)}
#######################################################################
logpostNorm.norm <- function( obj, param = obj$param,
   prior = obj$prior, prior.df = obj$prior.df,
   prior.sscp = obj$prior.sscp, ...){
   ##################################
   # S3 method for class "norm"
   ##################################
   result <- logpostNorm.default( obj$y, x = obj$x,
      intercept = FALSE, param = param,
      prior = prior, prior.df = prior.df, prior.sscp = prior.sscp)
   return(result)}
#######################################################################
logpostNorm.formula <- function( formula, data, param, 
   prior="uniform", prior.df=NULL, prior.sscp=NULL, ...){
   ###############################################
   # S3 method for formula object
   ###############################################
   # Prepare the model frame
   mf <- match.call(expand.dots=FALSE)
   mf[[1]] <- as.name("model.frame")
   mf$param <- mf$prior <- mf$prior.df <-
      mf$prior.sscp <- NULL
   m <- match( c("formula","data"), names(mf), 0L)
   mf <- mf[c(1L,m)]
   mf$na.action <- as.name("na.pass")
   mf <- eval( mf, parent.frame() )
   ###############################################
   # Extract the matrix of response variables
   y <- model.response( mf )
   y <- data.matrix( y )
   storage.mode(y) <- "double"
   if( ncol(y) == 1 ) colnames(y) <- all.vars(formula)[1]
   ###############################################
   # Extract the matrix of predictors
   x <- model.matrix( formula, mf )
   if( any( is.na(x) ) ){
      stop( "Missing values in predictors not allowed." )
      }
   ###############################################
   result <- logpostNorm.default( y, x = x,
      intercept = FALSE, param = param, 
      prior = prior, prior.df = prior.df, prior.sscp = prior.sscp)
   return(result)}
#######################################################################
#######################################################################
#######################################################################
#######################################################################
loglikNorm <- function(obj, ...){
   # S3 generic function
   UseMethod("logpostNorm")}
#######################################################################
loglikNorm.default <- function(obj, x=NULL, intercept=TRUE,
   param, ...){
   ##################################
   y <- obj
   if( class(y) == "data.frame" ){
      status <- logical( length(y) )
      for( j in 1:length(y) ) status[j] <- is.factor( y[[j]] )
      if( any(status) ){ 
         msg <- paste( "Factors in argument \"y\"",
               "converted to mode \"numeric\"." )
         warning( msg )
      }
   }  
   y <- data.matrix(y)
   storage.mode(y) <- "double"
   if( is.null( colnames(y) ) ){
      colnames(y) <- paste("Y", as.character( 1:ncol(y) ), sep="" )
   }
   ##################################
   if( mode(y) != "numeric" )
      stop("First argument should be of mode \"numeric\".")
   storage.mode(y) <- "double"
   if( is.null( colnames(y) ) ){
      colnames(y) <- paste("Y", as.character( 1:ncol(y) ), sep="" )
   }
   #########################################
   # get matrix of predictors
   if( is.null(x) ){
      # default to a column of ones
      x <- matrix( 1., nrow(y), 1)
      colnames(x) <- "CONST"
      }
   else{
      # check dimensions of x and coerce to a numeric matrix
      if( class(x) == "data.frame" ){
         status <- logical( length(x) )
         for( j in 1:length(x) ) status[j] <- is.factor( x[[j]] )
         if( any(status) ){ 
         msg <- paste( "Factors in argument \"x\"",
               "converted to mode \"numeric\"." )
         warning( msg )
         }
      }  
      x <- data.matrix(x)
      if( nrow(x) != nrow(y) ){
         stop("Arguments x and y do not have the same number of rows.")
      }
      if( any( is.na(x) ) ){
         stop("Missing values in x are not allowed.")
      }
      if( is.null( colnames(x) )) {
      	    colnames(x) <- paste("X", as.character( 1:ncol(x) ),
	    sep="" )
      }
      if( intercept ){
         x <- cbind( CONST=1., x )
      }
   }
   rownames(x) <- rownames(y)
   storage.mode(x) <- "double"
   #########################################
   # missing-value code to pass to Fortran
   mvcode <- .Machine$double.xmax
   y[ is.na(y) ] <- mvcode
   #########################################
   # create a matrix for holding message codes
   msg.len.max <- 40L
   msg.codes <- matrix( 0L, msg.len.max, 8L )
   #########################################
   # check user-supplied parameter values
   beta <- data.matrix( param$beta )
   sigma <- data.matrix( param$sigma )
   if( ( nrow(beta)!=ncol(x) ) | (ncol(beta)!=ncol(y)) )
      stop("Incorrect dimensions for argument beta.")
   if( ( nrow(sigma)!=ncol(y) ) | (ncol(sigma)!=ncol(y)) )
      stop("Incorrect dimensions for argument sigma.")
   storage.mode(beta) <- "double"
   storage.mode(sigma) <- "double"
   rownames( beta ) <- colnames( x )
   colnames( beta ) <- colnames( y )
   rownames( sigma ) <- colnames( y )
   colnames( sigma ) <- colnames( y )
   #########################################
   # uniform prior
   prior.type.int <- 1L
   storage.mode( prior.type.int ) <- "integer"
   prior.df <- -( ncol(y) + 1 ) 
   prior.sscp <- matrix( 0., ncol(y), ncol(y) )
   storage.mode( prior.df ) <- "double"
   storage.mode( prior.sscp ) <- "double"
   #########################################
   tmp <- .Fortran("norm_logpost",
         n = nrow(y),
         r = ncol(y),
         p = ncol(x),
         x = x,
         y = y,
         mvcode = mvcode,
         beta = beta,
         sigma = sigma,
         prior.type.int = prior.type.int,
         prior.df = prior.df,
         prior.sscp = prior.sscp,
	 logpost = numeric(1), 
         status = integer(1),
         msg.len.max = msg.len.max,
         msg.codes = msg.codes,
         msg.len.actual = integer(1),
         PACKAGE="norm2"
         )
   #########################################
   # display message from Fortran
   msg.lines <- msgNorm( tmp$msg.codes, tmp$msg.len.actual )
   if( is.null( msg.lines ) ){
      msg <- "OK"
   }
   else{
      msg <- paste0( msg.lines, collapse="\n" )
   }
   msg <- paste( msg, "\n", sep="")
   if( msg!= "OK\n" ) cat(msg)
   #########################################
   if( tmp$status != 0 ){
      result <- invisible(NULL)
   }
   else{
      result <- tmp$logpost
   }
   return(result)}
#######################################################################
loglikNorm.norm <- function( obj, param = obj$param, ...){
   ##################################
   # S3 method for class "norm"
   ##################################
   result <- loglikNorm.default( obj$y, x = obj$x,
      intercept = FALSE, param = param)
   return(result)}
#######################################################################
loglikNorm.formula <- function( formula, data, param, ...){
   ###############################################
   # S3 method for formula object
   ###############################################
   # Prepare the model frame
   mf <- match.call(expand.dots=FALSE)
   mf[[1]] <- as.name("model.frame")
   mf$param <- NULL
   m <- match( c("formula","data"), names(mf), 0L)
   mf <- mf[c(1L,m)]
   mf$na.action <- as.name("na.pass")
   mf <- eval( mf, parent.frame() )
   ###############################################
   # Extract the matrix of response variables
   y <- model.response( mf )
   y <- data.matrix( y )
   storage.mode(y) <- "double"
   if( ncol(y) == 1 ) colnames(y) <- all.vars(formula)[1]
   ###############################################
   # Extract the matrix of predictors
   x <- model.matrix( formula, mf )
   if( any( is.na(x) ) ){
      stop( "Missing values in predictors not allowed." )
      }
   ###############################################
   result <- loglikNorm.default( y, x = x,
      intercept = FALSE, param = param)
   return(result)}
#######################################################################
#######################################################################
#######################################################################
#######################################################################
print.norm <- function( x, ... ){
   ###############################
   # S3 method for class "norm"
   ###############################
   cat( "This is an object of class \"norm\";\n" )
   cat( "use summary() to view its contents.\n" )
   return( invisible(x) )}
#######################################################################
summary.norm <- function( object,
   show.variables = (object$method=="EM"), 
   show.patterns = (object$method=="EM"),
   show.params = (object$method=="EM"), ... ){
   ##########################################
   # S3 method for class "norm"
   ##########################################
   if( class(object) != "norm" )
      stop("Argument should be of class \"norm.\"")
   result <- object
   ###############################################
   result$show.variables <- show.variables
   result$show.patterns <- show.patterns
   result$show.params <- show.params
   ###############################################
   result$miss.patt.strings <- 
      turn.miss.patt.into.strings( object$miss.patt )
   ###############################################
   # summary table for X variables
   x.table <- matrix( numeric(1), ncol(object$x), 5 )
   rownames(x.table) <- colnames(object$x)
   colnames(x.table) <- c("Mean", "SD", "Observed", 
      "Missing", "Pct.Missing")
   x.table[,1] <- apply( object$x, 2, mean )
   x.table[,2] <- sqrt( diag( var( object$x ) ) )
   x.table[,3] <- nrow( object$x )
   x.table[,4] <- 0
   x.table[,5] <- 0.
   result$x.table <- x.table
   # summary table for Y variables
   y.table <- matrix( numeric(1), ncol(object$y), 5 )
   rownames(y.table) <- colnames( object$y )
   colnames(y.table) <- c("Mean", "SD", "Observed", 
      "Missing", "Pct.Missing")
   y.table[,1] <- object$ybar
   y.table[,2] <- object$ysdv
   y.table[,3] <- object$n.obs
   y.table[,4] <- nrow( object$y ) - object$n.obs
   y.table[,5] <- 100*( nrow( object$y ) - object$n.obs ) / nrow(object$y)
   result$y.table <- y.table
   ###############################################
   # summary of EM or MCMC run
   if( object$method == "EM" ){
      # estimate componentwise rate of convergence from last
      # two iterations
      if( object$iter < 2 ){
         conv.rate <- NULL
      }
      else{
         # estimate rate of convergence
         tmpA <- object$param.rate$beta
         tmpB <- object$param.rate$sigma
         tmpB <- tmpB[ row(tmpB) >= col(tmpB) ]
         tmp <- c(tmpA,tmpB)
         if( all( tmp == 0. ) ){
            conv.rate <- 0.
         }
         else{
            conv.rate <- round( median( tmp[(tmp>0)&(tmp<1)] ), 5)
            if( is.na(conv.rate) ) conv.rate <- NULL
         }
      }
      # get worst fraction of missing information from power method
      if( object$em.worst.ok ){
         worst.frac <- object$worst.frac
         }
      else{
         worst.frac <- NA
         }
      em.summary <- c(
         "Method:" = "EM",
         "Prior:" = paste("\"", object$prior, "\"", sep=""),
         "Prior df:" = object$prior.df,
         "Convergence criterion:" = format(object$criterion),
         "Iterations:" = format( object$iter ),
         "Converged:" = format( object$converged ),
         "Max. rel. difference:" = format(object$rel.diff, digits=5),
         "-2 Loglikelihood:" = format(-2.*object$loglik[object$iter]),
         "-2 Log-posterior density:" = format(-2.*object$logpost[object$iter]),
         "Worst fraction missing information:" = round(worst.frac,4) )
      result$em.summary <- em.summary
   }
   else if( object$method == "MCMC" ){
      mcmc.summary <- c(
         "Method:" = "MCMC",
         "Prior:" = paste("\"", object$prior, "\"", sep=""),
         "Prior df:" = object$prior.df,
         "Iterations:" = format( object$iter ),
         "Cycles per iteration:" = format(object$multicycle),
         "Impute every k iterations, k =" = format(object$impute.every),
         "No. of imputations created:" = format(length(object$imp.list)),
         "series.worst present:" = format( !is.null(object$series.worst) ),
         "series.beta  present:" = format( !is.null(object$series.beta ) ),
         "series.sigma present:" = format( !is.null(object$series.sigma) ) )
      result$mcmc.summary <- mcmc.summary
   }
   #
   # set class and return
   class(result) <- "summary.norm"
   return( result )}
#######################################################################
# private function used by print.summary.norm
turn.miss.patt.into.strings <- function(miss.patt){
   result <- ""
   for( j in 1:ncol(miss.patt) ){
      tmp <- rep(".", nrow(miss.patt) )
      tmp[ miss.patt[,j] ] <- "m"
      result <- paste(result, tmp, sep="")
      }
   return(result)}
#######################################################################
print.summary.norm <- function( x, ... ){ 
   #########################################
   # S3 method for class "summary.norm"
   #########################################
   if( class(x) != "summary.norm" )
      stop("Argument should be of class summary.norm")
   # 
   if( x$show.variables ){
      # print x.table
      cat("Predictor (X) variables:\n")
      print( x$x.table )
      # print y.table
      cat("\n")
      cat("Response (Y) variables:\n")
      print( x$y.table )
   }
   #
   if( x$show.patterns ){
      # print missing pattern strings
      cat("\n")
      cat("Missingness patterns for response (Y) variables\n")
      cat("   (. denotes observed value, m denotes missing value)\n")
      cat("   (variable names are displayed vertically)\n")
      cat("   (rightmost column is the frequency):\n")
      tmp <- format( colnames(x$y), justify="right")
      header <- ""
      for(i in 1:ncol( x$y )) header <- paste(header, " ", sep="")
      header <- rep( header, nchar( tmp[1] ) )
      for( i in 1: nchar( tmp[1] ) ){
         for( j in 1:ncol( x$y ) ){
         substring( header[i], j, j ) <- substring( tmp[j], i, i )
         }
      }
      cat( header, sep="\n")
      cat( paste(x$miss.patt.strings, 
          format(x$miss.patt.freq)), sep="\n")
   }
   #
   if( x$method == "EM" ){
      cat("\n")
      cat( paste( format(names(x$em.summary)), 
         x$em.summary ), sep="\n")
   }
   else if( x$method == "MCMC" ) {
      cat("\n")
      cat( paste( format(names(x$mcmc.summary)),
         x$mcmc.summary ), sep="\n")
   }
   #
   if( x$msg != "OK\n"){
      cat("\n")
      cat( x$msg )
   }
   #
   if( x$show.params ){
      # print beta
      cat("\n")
      cat("Estimated coefficients (beta):\n")
      print( x$param$beta )
      # print sigma
      cat("\n")
      cat("Estimated covariance matrix (sigma):\n")
      print( x$param$sigma )
   }
   return( invisible(x) )}
#######################################################################
msgNorm <- function( msg.codes, msg.len.actual ){
   ###########################################
   # private function 
   # converts matrix of integer message codes
   # to character text;
   # relies on internal package data object
   # icodesNorm
   ###########################################
   if( msg.len.actual == 0 ){
      msg.lines <- NULL
   }
   else{
      msg.lines <- character( msg.len.actual )
      for( i in 1:msg.len.actual ){
         if( msg.codes[i,1L] == 1L ){
            msg.lines[i] <- icodesNorm$comments[ msg.codes[i,4L] ]
         }
         else if( msg.codes[i,1L] == 2L ){
            msg.lines[i] <- paste( "OCCURRED IN:", 
               icodesNorm$subnames[ msg.codes[i,3L] ],
               "in MOD",
               icodesNorm$modnames[ msg.codes[i,2L] ],
               sep=" " )
         }
         else if( msg.codes[i,1L] == 3L ){
            msg.lines[i] <- paste( "Observation", 
   	    format( msg.codes[i,5L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 4L ){
            msg.lines[i] <- paste( "Variable", 
   	    format( msg.codes[i,6L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 5L ){
            msg.lines[i] <- paste( "Iteration", 
   	    format( msg.codes[i,7L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 6L ){
            msg.lines[i] <- paste( "Iteration ", 
   	    format( msg.codes[i,7L] ),
               ", Cycle ",
   	    format( msg.codes[i,8L] ),
               sep="" )
         }
         else{
            msg.lines[i] <- "???"
         }
      }
   }
   return( msg.lines ) }
#######################################################################
#######################################################################
#######################################################################

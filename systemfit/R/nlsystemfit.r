###   $Id: nlsystemfit.r 1132 2015-06-08 20:51:21Z arne $
###
###            Simultaneous Nonlinear Least Squares for R
###
### Copyright 2003-2004 Jeff D. Hamann <jeff.hamann@forestinformatics.com>
###
### This file is part of the nlsystemfit library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

### uses the Dennis + Schnabel Minimizer which is the one utilized by R's nlm()

### remove before release
# rm(list=ls(all=TRUE))

### this function is the "driver" function for the minimization...
knls <- function( theta, eqns, data, fitmethod="OLS", parmnames, instr=NULL, S=NULL ) {

  r  <- matrix()               # residuals equation wise
  r <- NULL

  gmm.resids <- matrix()
  gmm.resids <- NULL

  residi  <- list()               # residuals equation wise
  lhs <- list()
  rhs <- list()
  neqs <- length( eqns )
  nobs <- dim( data )[[1]]            # number of nonmissing observations

  ## GMM specific variables, in this case... g = 2, k = 3
#  V <- matrix( 0, g*k, g*k )          # V is a 6x6 matrix
  moments <- list()
  mn <- array()

  moments <- NULL
  mn <- NULL
  lhs <- NULL
  rhs <- NULL
  residi <- NULL
  # partial derivatives of the residuals with respect to the parameters
  dResidTheta <- NULL
  dResidThetai <- list()
  d2ResidTheta <- array( NA, c( 0, 0, 0 ) )
  d2ResidThetai <- list()

  ## get the values of the parameters
  for( i in 1:length( parmnames ) ) {
    name <- names( parmnames )[i]
    val <- theta[i]
    storage.mode( val ) <-  "double"
    assign( name, val )
  }

  ## build the residual vector...
  for( i in 1:length( eqns ) ) {
    lhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[2]], envir = data ) )
    rhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[3]], envir = data ) )
    residi[[i]] <- lhs[[i]] - rhs[[i]]
    r <- rbind( r, as.matrix( residi[[i]] ) )
    if( fitmethod == "GMM" ) {
      gmm.resids <- cbind( gmm.resids, as.matrix( residi[[i]] ) )
    }
    dResidThetai[[ i ]] <- - attributes( with( data, with( as.list( theta ),
      eval( deriv( eqns[[i]], names( parmnames )), envir = data ))))$gradient
    dResidTheta <- rbind( dResidTheta, dResidThetai[[ i ]] )
    d2ResidThetai[[ i ]] <- - attributes( with( data, with( as.list( theta ),
      eval( deriv3( eqns[[i]], names( parmnames ) ), envir = data ))))$hessian
    temp <- array( NA, c( dim( d2ResidTheta )[ 1 ] +
      dim( d2ResidThetai[[ i ]] )[ 1 ], dim( d2ResidThetai[[ i ]] )[ 2:3 ] ) )
    if( i > 1 ) {
      temp[ 1:dim( d2ResidTheta )[ 1 ], , ] <- d2ResidTheta
    }
    temp[ ( dim( d2ResidTheta )[ 1 ] + 1 ):( dim( temp )[ 1 ] ),
      , ] <- d2ResidThetai[[ i ]]
    d2ResidTheta <- temp
  }
  ## these are the objective functions for the various fitting methods
  if( fitmethod == "OLS" ) {
    obj <- crossprod( r )
    gradient <- 2 * t( r ) %*% dResidTheta
#print( d2ResidTheta )
#print( t( dResidTheta ) %*% dResidTheta )
    hessian <- matrix( NA, nrow = length( parmnames ), ncol = length( parmnames ) )
    for( i in 1:length( parmnames ) ) {
      hessian[ i, ] <- 2 * t( r ) %*% d2ResidTheta[ , , i ]
    }
    hessian <- hessian + 2 * t( dResidTheta ) %*% dResidTheta
    rownames( hessian ) <- colnames( dResidTheta )
    colnames( hessian ) <- colnames( dResidTheta )
#print( hessian - t( hessian ) )
    hessian <- NULL
    attributes( obj ) <- list( gradient = gradient, hessian = hessian )
  }
  if( fitmethod == "2SLS" ) {
    ## W is premultiplied == ( diag( neqs ) %x% W )
    ##obj <- ( t(r) %*% S %*% r )
    obj <- crossprod(t(crossprod(r,S)),r)
    attributes( obj ) <- list( gradient = 2 * t( r ) %*% S %*% dResidTheta )
  }
  if( fitmethod == "SUR" ) {
    ## S is premultiplied == ( qr.solve( S ) %x% diag( nobs ) )
    ##obj <- ( t(r) %*% S %*% r )
    obj <- crossprod(t(crossprod(r,S)),r)
    attributes( obj ) <- list( gradient = 2 * t( r ) %*% S %*% dResidTheta )
  }
  if( fitmethod == "3SLS" ) {
    ## S is premultiplied == ( qr.solve( S ) %x% W )
    ##obj <- ( t(r) %*% S %*% r )
    obj <- crossprod(t(crossprod(r,S)),r)
    attributes( obj ) <- list( gradient = 2 * t( r ) %*% S %*% dResidTheta )
  }
  if( fitmethod == "GMM" ) {
    ## this just can't be correct... or can it...
    ## S is a gx x gk matrix
    ## g = number of eqns, k = number of inst variables
    z <- as.matrix( model.frame( instr, data = data ) )
    for(t in 1:nobs) {
      moments <- rbind( moments, gmm.resids[t,] %x% z[t,] )
    }
    g <- length( eqns )                 # number of equations
    k <- dim( as.matrix( model.frame( instr, data = data ) ) )[[2]]
    gk <- g*k
    for( i in 1:gk ) {
      mn <- rbind( mn, mean( moments[,i] ) )
    }
    ##obj <- ( t(nobs*mn) %*% S %*% (nobs*mn) ) / nobs
    ##obj <- ( t(mn) %*% S %*% (mn) )
    obj <- crossprod(t(crossprod(mn,S)),mn)
  }

  ## it would be nice to place the gradient and/or hessian attributes...
  ## how can I make this work???
  ## attr( obj, "gradient" ) <- "hi mom"
  ## attr( obj, "hessian" ) <- hessian...

  return( obj )
}


nlsystemfit <- function( method="OLS",
                        eqns,
                        startvals,
                        eqnlabels=c(as.character(1:length(eqns))),
                        inst=NULL,
                        data=list(),
                        solvtol=.Machine$double.eps,
                        maxiter=1000, ... ) {

  ## some tests
  if(!(method=="OLS" | method=="SUR" | method=="2SLS" | method=="3SLS" | method=="GMM" )){
    stop("The method must be 'OLS', 'SUR', '2SLS', or '3SLS'")}
  if((method=="2SLS" | method=="3SLS" | method=="GMM") & is.null(inst)) {
    stop("The methods '2SLS', '3SLS' and 'GMM' need instruments!")}

  lhs <- list()
  rhs <- list()
  derivs <- list()

  results <- list()               # results to be returned
  results$eq <- list()            # results for the individual equations
  resulti <- list()               # results of the ith equation
  residi  <- list()               # residuals equation wise
  iter    <- NULL                 # number of iterations
  G       <- length( eqns )       # number of equations
  n       <- array( 0, c(G))      # number of observations in each equation
  k       <- array( 0, c(G) )     # number of (unrestricted) coefficients/regressors in each equation
  df       <- array( 0, c(G) )     # degrees of freedom in each equation
  instl   <- list()               # list of the instruments for each equation
  ssr     <- array( 0, c(G))      # sum of squared residuals of each equation
  mse     <- array( 0, c(G))      # mean square error (residuals) of each equation
  rmse    <- array( 0, c(G))      # root of mse
  r2      <- array( 0, c(G))      # R-squared value
  adjr2   <- array( 0, c(G))      # adjusted R-squared value
  nobs <- dim( data )[[1]]
  S <- matrix( 0, G, G )               # covariance matrix of the residuals
  X <- array()
  x <- list()

  resids <- array()
  resids <- NULL

  if( method == "OLS" ) {
    if( TRUE ) {
      est <- nlm( knls, startvals,
               typsize=abs(startvals),iterlim=maxiter,
               eqns=eqns,
               data=data,
               fitmethod=method,
               parmnames=startvals,
               ... )
    } else {
      est <- optim( fn = knls, par = startvals,
               eqns=eqns,
               data=data,
               fitmethod=method,
               parmnames=startvals )
    }
  }
  if( method == "2SLS" ) {
    ## just fit and part out the return structure
    z <- as.matrix( model.frame( inst, data = data ) )
    Wt <- z %*% qr.solve( crossprod( z ), tol=solvtol ) %*% t(z)
    W2 <- diag( length( eqns ) ) %x% Wt
    est <- nlm( knls, startvals,
               typsize=abs(startvals),iterlim=maxiter,
               eqns=eqns,
               data=data,
               fitmethod=method,
               parmnames=startvals,
               S=W2, ... )
  }
  if( method == "SUR" || method == "3SLS" || method == "GMM" ) {
    ## fit ols/2sls, build the cov matrix for estimation and refit
    if( method == "SUR" ) {
      estols <- nlm( knls, startvals,
                    typsize=abs(startvals),iterlim=maxiter,
                    eqns=eqns,
                    data=data,
                    fitmethod="OLS",
                    parmnames=startvals, ... )
    }
    if( method == "3SLS" || method == "GMM" ) {
      z <- as.matrix( model.frame( inst, data = data ) )
      W <- z %*% qr.solve( crossprod( z ), tol=solvtol ) %*% t(z)
      W2 <- ( diag( length( eqns ) ) %x% W )
      estols <- nlm( knls, startvals,
                    typsize=abs(startvals),iterlim=maxiter,
                    eqns=eqns,
                    data=data,
                    fitmethod="2SLS",
                    parmnames=startvals,
                    instr=inst,
                    S=W2, ... )
    }

    ## build the S matrix
    names( estols$estimate ) <- names( startvals )
    for( i in 1:length( estols$estimate ) ) {
      name <- names( estols$estimate )[i]
      val <- estols$estimate[i]
      storage.mode( val ) <-  "double"
      assign( name, val )
    }

    ## get the rank for the eqns, compute the first-stage
    ## cov matrix to finish the SUR and 3SLS methods
    for(i in 1:G) {
      lhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[2]], envir = data ) )
      rhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[3]], envir = data ) )
      residi[[i]] <- lhs[[i]] - rhs[[i]]
      derivs[[i]] <- deriv( as.formula( eqns[[i]] ), names( startvals ) )
      ## computing the jacobian to get the rank to get the number of variables...
      jacobian <- attr( eval( derivs[[i]], envir = data ), "gradient" )
      n[i]   <-  length( lhs[[i]] )
      k[i] <- qr( jacobian )$rank
      df[i] <- n[i] - k[i]
    }

    ## covariance matrix of the residuals from the first stage...
    Sols <- matrix( 0, G, G )
    rcovformula <- 1
    for(i in 1:G) {
      for(j in 1:G) {
        Sols[i,j] <- sum(residi[[i]]*residi[[j]])/(
                                                   sqrt((n[i]-rcovformula*k[i])*(n[j]-rcovformula*k[j])))
      }
    }

    if( method == "SUR" ) {
      Solsinv <- qr.solve( Sols, tol=solvtol ) %x% diag( nobs )
      est <- nlm( knls,estols$estimate,
                 typsize=abs(estols$estimate),iterlim=maxiter,
                 eqns=eqns, data=data, fitmethod=method, parmnames=startvals,
                 S=Solsinv, ... )
    }
    if( method == "3SLS" ) {
      z <- as.matrix( model.frame( inst, data = data ) )
      W <- z %*% qr.solve( crossprod( z ), tol=solvtol ) %*% t(z)
      Solsinv <- qr.solve( Sols, tol=solvtol ) %x% W
      est <- nlm( knls, estols$estimate,
                 typsize=abs(estols$estimate),iterlim=maxiter,
                 eqns=eqns, data=data, fitmethod=method, parmnames=startvals,
                 S=Solsinv, instr=z, ... )
    }
    if( method == "GMM" ) {
      resids <- NULL
      for(i in 1:G) {
        resids <- cbind( resids, residi[[i]] )
      }
      z <- as.matrix( model.frame( inst, data = data ) )
      moments <- list()
      moments <- NULL
      for(t in 1:nobs) {
        moments <- rbind( moments, resids[t,] %x% z[t,] )
      }
      v2sls <- qr.solve( var( moments ), tol=solvtol )
      est <- nlm( knls,estols$estimate,
                 typsize=abs(estols$estimate),iterlim=maxiter,
                 eqns=eqns, data=data, fitmethod="GMM", parmnames=startvals,
                 S=v2sls, instr=inst, ... )
    }
  }

  ## done with the fitting...
  ## now, part out the results from the nlm function
  ## to rebuild the equations and return object
  ## get the parameters for each of the equations and


  ## evaluate the residuals for eqn
  ## get the values of the final parameters
  if( TRUE ) {
    estimate <- est$estimate
  } else {
    estimate <- est$par
  }
  names( estimate ) <- names( startvals )
  for( i in 1:length( estimate ) ) {
    name <- names( estimate )[i]
    ### I wonder if I need to clear out the variables before assigning them for good measure...
    assign( name, NULL )
    val <- estimate[i]
    storage.mode( val ) <-  "double"
    assign( name, val )
  }

  ## get the rank for the eqns, compute the first-stage
  ## cov matrix to finish the SUR and 3SLS methods
  X <- NULL
  results$resids <- array()
  results$resids <- NULL

  ## you're working on parsing out the parameters and the estiamtes for the return structure...
  for(i in 1:G) {
    lhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[2]], envir = data ) )
    rhs[[i]] <- as.matrix( eval( as.formula( eqns[[i]] )[[3]], envir = data ) )
    residi[[i]] <- lhs[[i]] - rhs[[i]]
    derivs[[i]] <- deriv( as.formula( eqns[[i]] ), names( startvals ) )

    ## computing the jacobian to get the rank to get the number of variables...
    jacobian <- attr( eval( derivs[[i]], envir = data ), "gradient" )
    n[i]   <-  length( lhs[[i]] )
    k[i] <- qr( jacobian )$rank
    df[i] <- n[i] - k[i]
    ssr[i] <- crossprod( residi[[i]] )
    mse[i] <- ssr[i] / ( n[i] - k[i] )
    rmse[i] <- sqrt( mse[i] )

    X <- rbind( X, jacobian )
    results$resids <- cbind( results$resids, as.matrix( residi[[i]] ) )
  }

  ## compute the final covariance matrix
  ## you really should use the code below to handle weights...
  rcovformula <- 1
  for(i in 1:G) {
    for(j in 1:G) {
      S[i,j] <- sum(residi[[i]]*residi[[j]])/(
                                              sqrt((n[i]-rcovformula*k[i])*(n[j]-rcovformula*k[j])))
    }
  }

### for when you get the weights working...
#     vardef <- 1
#     if( vardef == 1 ) {
#       D <- diag( G ) * 1 / sqrt( nrow( data ) )
#     }
#     if( vardef == 2 ) {
#       D <- diag( G ) * 1 / sqrt( sum( weights ) )
#     }
#     if( vardef == 3 ) {
#       D <- diag( G ) * 1 / sqrt( sum( weights ) - ( sum( n ) - sum( k ) ) )
#     }
#     if( vardef == 4 ) {
#       for(i in 1:G) {
#         D <- diag( G )
#         D[i,i] <- D[i,i] * 1 / sqrt( nrow( data ) - n[i] )
#       }
#     }
#     ## the docs have this, but the table contains the above equations
#     R <- crossprod( results$resids )
#     S <- D %*% R %*% D
#     SI <- qr.solve( S, tol=solvtol ) %x% diag( nrow( data ) )
#     covb <- qr.solve(t(X) %*% SI %*% X, tol=solvtol )



  ## get the variance-covariance matrix
  if( method == "OLS" ) {
    SI <- diag( diag( qr.solve( S, tol=solvtol ) ) ) %x% diag( nrow( data ) )
    covb <- qr.solve(t(X) %*% SI %*% X, tol=solvtol )
  }
  if( method == "2SLS" ) {
    Z <- model.matrix(inst, data = data )
    W <- Z %*% qr.solve( crossprod( Z ), tol=solvtol ) %*% t(Z)
    SW <- diag( diag( qr.solve( S, tol=solvtol ) ) ) %x% W
    covb <- qr.solve(t(X) %*% SW %*% X, tol=solvtol )
  }
  if( method == "SUR" ) {
    SI <- qr.solve( S, tol=solvtol ) %x% diag( nrow( data ) )
    covb <- qr.solve(t(X) %*% SI %*% X, tol=solvtol )
  }
  if( method == "3SLS" ) {
    Z <- model.matrix(inst, data = data )
    W <- Z %*% qr.solve( crossprod( Z ), tol=solvtol ) %*% t(Z)
    SW <- qr.solve( S, tol=solvtol ) %x% W
    covb <- qr.solve(t(X) %*% SW %*% X, tol=solvtol )
  }
  if( method == "GMM" ) {
#    print( "obtaining GMM(SE)" )
    z <- as.matrix( model.frame( inst, data = data ) )
    moments <- list()
    moments <- NULL
    resids <- NULL
    for(i in 1:G) {
      resids <- cbind( resids, residi[[i]] )
    }
    for(t in 1:nobs) {
      moments <- rbind( moments, resids[t,] %x% z[t,] )
    }
#    print( var( moments ) )
    Vinv <- qr.solve( var( moments ), tol=solvtol )
#    print( Vinv )
    Y <- diag( length( eqns ) ) %x% t(z)
#    print( "covb now..." )
#    print( dim( Y ) )
#    print( dim( X ) )
    covb <- qr.solve( t( Y %*% X ) %*% Vinv %*% ( Y %*% X  ), tol=solvtol )
#    print( covb )
  }
  colnames( covb ) <- rownames( covb )


  ## bind the standard errors to the parameter estimate matrix
  se2 <- sqrt( diag( covb ) )
  t.val <- estimate / se2
  prob  <- 2*( 1 - pt( abs( t.val ), sum( n ) - sum( k ) ) ) ### you better check this...

  results$method       <- method
  results$n <- sum( n )
  results$k <- sum( k )
  results$b <- estimate
  results$se <- se2
  results$t <- t.val
  results$p <- prob

  ## build the results structure...
  for(i in 1:G) {
    ## you may be able to shrink this up a little and write the values directly to the return structure...
    eqn.terms <- vector()
    eqn.est <- vector()
    eqn.se <- vector()
    jacob <- attr( eval( deriv( as.formula( eqns[[i]] ), names( startvals ) ),
       envir = data ), "gradient" )
    for( v in 1:length( estimate ) ) {
      if( qr( jacob[,v] )$rank > 0 ) {
        eqn.terms <- rbind( eqn.terms, name <- names( estimate )[v] )
        eqn.est <- rbind( eqn.est, estimate[v] )
        eqn.se <- rbind( eqn.se, se2[v] )
      }
    }


    ## build the "return" structure for the equations
    resulti$method       <- method
    resulti$i            <- i               # equation number
    resulti$eqnlabel     <- eqnlabels[[i]]
    resulti$formula      <- eqns[[i]]
    resulti$b <- as.vector( eqn.est )
    names( resulti$b )   <- eqn.terms
    resulti$se           <- eqn.se
    resulti$t            <- resulti$b / resulti$se
    resulti$p            <- 2*( 1-pt(abs(resulti$t), n[i] - k[i] ))
    resulti$n            <- n[i]            # number of observations
    resulti$k            <- k[i]            # number of coefficients/regressors
    resulti$df           <- df[i]           # degrees of freedom of residuals
    resulti$predicted    <- rhs[[i]]           # predicted values
    resulti$residuals    <- residi[[i]]     # residuals
    resulti$ssr          <- ssr[i]             # sum of squared errors/residuals
    resulti$mse          <- mse[i]             # estimated variance of the residuals (mean squared error)
    resulti$s2           <- mse[i]             #        the same (sigma hat squared)
    resulti$rmse         <- rmse[i]            # estimated standard error of the residuals
    resulti$s            <- rmse[i]            #        the same (sigma hat)

#     ## you'll need these to compute the correlations...
#     print( paste( "eqn ", i ) )
    coefNames <- rownames( covb )[ rownames( covb ) %in%
      strsplit( as.character( eqns[[ i ]] )[ 3 ], "[^a-zA-Z0-9.]" )[[ 1 ]] ]
    resulti$covb <- covb[ coefNames, coefNames ]

#     resulti$x <- model.frame( as.formula( eqns[[i]] )[[3]], data = data )
#     print( resulti$x )
#    print( model.frame( eval( eqns[[i]], envir = data ), data = data ) )



    ## fix this to allow for multiple instruments?
    if(method=="2SLS" | method=="3SLS" | method=="GMM") {
      resulti$inst         <- inst
      ##resulti$inst         <- inst[[i]]
      ##resulti$inst         <- instl[[i]]
      ## resulti$h            <- h[[i]]          # matrix of instrumental variables
    }

    resulti$r2     <- 1 - ssr[i] / ( ( crossprod( lhs[[i]]) ) - mean( lhs[[i]] )^2 * nobs )
    resulti$adjr2  <- 1 - ((n[i]-1)/df[i])*(1-resulti$r2)

    class(resulti)        <- "nlsystemfit.equation"
    results$eq[[i]]      <- resulti
  }

  results$solvtol <- solvtol
  results$covb <- covb
  results$rcov <- S
  results$rcor <- cor( results$resids )
  results$drcov <- det(results$rcov)          # det(rcov, tol=solvetol)

  if(method=="2SLS" || method=="3SLS") {
    ##      results$h       <- H            # matrix of all (diagonally stacked) instrumental variables
  }
  if(method=="SUR" || method=="3SLS" || method=="GMM" ) {
    results$rcovest <- Sols      # residual covarance matrix used for estimation
    ##results$mcelr2  <- mcelr2       # McElroy's R-squared value for the equation system
  }

  ## build the "return" structure for the whole system
  results$method  <- method
  results$g       <- G              # number of equations
  results$nlmest  <- est

  class(results)  <- "nlsystemfit.system"

  if( results$nlmest$code >= 4 ) {
    warning( "Estimation did not converge!" )
  }

  return( results )
}


## print the (summary) results that belong to the whole system
summary.nlsystemfit.system <- function(object,...) {
  summary.nlsystemfit.system <- object
  summary.nlsystemfit.system
}


## print the results that belong to the whole system
print.nlsystemfit.system <- function( x, digits=6,... ) {
  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  table <- NULL
  labels <- NULL

  cat("\n")
  cat("nlsystemfit results \n")
  cat("method: ")

#  if(!is.null(object$iter)) if(object$iter>1) cat("iterated ")
  cat( paste( object$method, "\n\n"))
#   if(!is.null(object$iter)) {
#     if(object$iter>1) {
#       if(object$iter<object$maxiter) {
#         cat( paste( "convergence achieved after",object$iter,"iterations\n\n" ) )
#       } else {
#         cat( paste( "warning: convergence not achieved after",object$iter,"iterations\n\n" ) )
#       }
#     }
#   }

  cat( paste( "convergence achieved after",object$nlmest$iterations,"iterations\n" ) )
  cat( paste( "nlsystemfit objective function value:",object$nlmest$minimum,"\n\n" ) )


  for(i in 1:object$g) {
    row <- NULL
    row <- cbind( round( object$eq[[i]]$n,     digits ),
                  round( object$eq[[i]]$df,    digits ),
                  round( object$eq[[i]]$ssr,   digits ),
                  round( object$eq[[i]]$mse,   digits ),
                  round( object$eq[[i]]$rmse,  digits ),
                  round( object$eq[[i]]$r2,    digits ),
                  round( object$eq[[i]]$adjr2, digits ))
    table  <- rbind( table, row )
    labels <- rbind( labels, object$eq[[i]]$eqnlabel )
  }
  rownames(table) <- c( labels )
  colnames(table) <- c("N","DF", "SSR", "MSE", "RMSE", "R2", "Adj R2" )

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE )
  cat("\n")

  ## check this code before release...
  if(!is.null(object$rcovest)) {
    cat("The covariance matrix of the residuals used for estimation\n")
    rcov <- object$rcovest
    rownames(rcov) <- labels
    colnames(rcov) <- labels
    print( rcov )
    cat("\n")
    if( min(eigen( object$rcovest, only.values=TRUE)$values) < 0 ) {
      cat("warning: this covariance matrix is NOT positive semidefinit!\n")
      cat("\n")
    }
  }

  cat("The covariance matrix of the residuals\n")
  rcov <- object$rcov
  rownames(rcov) <- labels
  colnames(rcov) <- labels
  print( rcov )
  cat("\n")

  cat("The correlations of the residuals\n")
  rcor <- object$rcor
  rownames(rcor) <- labels
  colnames(rcor) <- labels
  print( rcor )
  cat("\n")

  cat("The determinant of the residual covariance matrix: ")
  cat(object$drcov)
  cat("\n")

### check this code before release
#   cat("OLS R-squared value of the system: ")
#   cat(object$olsr2)
#   cat("\n")

#   if(!is.null(object$mcelr2)) {
#     cat("McElroy's R-squared value for the system: ")
#     cat(object$mcelr2)
#     cat("\n")
#   }

  ## now print the individual equations
  for(i in 1:object$g) {
    print( object$eq[[i]], digits )
  }

}


## print the (summary) results for a single equation
summary.nlsystemfit.equation <- function(object,...) {
  summary.nlsystemfit.equation <- object
  summary.nlsystemfit.equation
}


## print the results for a single equation
print.nlsystemfit.equation <- function( x, digits=6, ... ) {
  object <- x

  save.digits <- unlist(options(digits=digits))
  on.exit(options(digits=save.digits))

  cat("\n")
  cat( paste( object$method, " estimates for ", object$eqnlabel, " (equation ", object$i, ")\n", sep="" ) )

  cat("Model Formula: ")
  print(object$formula)
  if(!is.null(object$inst)) {
    cat("Instruments: ")
    print(object$inst)
  }
  cat("\n")

  Signif <- symnum(object$p, corr = FALSE, na = FALSE,
                   cutpoints = c(0,  .001,.01,.05, .1, 1),
                   symbols   = c("***","**","*","."," "))

  table <- cbind(round( object$b,  digits ),
                 round( object$se, digits ),
                 round( object$t,  digits ),
                 round( object$p,  digits ),
                 Signif)

  rownames(table) <- names(object$b)
  colnames(table) <- c("Estimate","Std. Error","t value","Pr(>|t|)","")

  ##print.matrix(table, quote = FALSE, right = TRUE )
  ##prmatrix(table, quote = FALSE, right = TRUE )
  print(table, quote = FALSE, right = TRUE )
  cat("---\nSignif. codes: ",attr(Signif,"legend"),"\n")

  cat(paste("\nResidual standard error:", round(object$s, digits),  ## s ist the variance, isn't it???
            "on", object$df, "degrees of freedom\n"))

  cat( paste( "Number of observations:", round(object$n, digits),
              "Degrees of Freedom:", round(object$df, digits),"\n" ) )

  cat( paste( "SSR:",     round(object$ssr,    digits),
              "MSE:", round(object$mse, digits),
              "Root MSE:",   round(object$rmse,  digits), "\n" ) )

   cat( paste( "Multiple R-Squared:", round(object$r2,    digits),
               "Adjusted R-Squared:", round(object$adjr2, digits),
               "\n" ) )
  cat("\n")
}

# from Model Selection and Inference: A Practical Information-Theoretic Approach
# Kenneth P. Burnham and David R. Anderson, 1998. Springer-Verlag, New York, New York.

## Akaike's Information Criterion
## AIC = n * log( sigmahat^2 ) + 2K
## n = number of obs
## sigmahat^2 = sum( error^2 ) / n == residual sums of squares
## K is the total number if estimated parameters, including the intercept and sigma^2 (nparams + 1)
## second order AIC
## AICc = AIC + (2K*(K+1))/(n-K-1)
## unless the sample size is large with repsect to the number of estiamted parameters, use AICc.


### This is a stripped down version of rq that
### only estimates betas without bells and whistles uses "br" method
shortrq.fit.br <- function (x, y, tau = 0.5)
{
    tol <- .Machine$double.eps^(2/3)
    eps <- tol
    big <- .Machine$double.xmax
    x <- as.matrix(x)
    p <- ncol(x)
    n <- nrow(x)
    ny <- NCOL(y)
    nsol <- 2
    ndsol <- 2
    # Check for Singularity of X since br fortran isn't very reliable about this
    if (qr(x)$rank < p)
        stop("Singular design matrix")
    lci1 <- FALSE
    qn <- rep(0, p)
    cutoff <- 0
 
    z <- .Fortran("rqbr", as.integer(n), as.integer(p), as.integer(n +
        5), as.integer(p + 3), as.integer(p + 4), as.double(x),
        as.double(y), as.double(tau), as.double(tol), flag = as.integer(1),
        coef = double(p), resid = double(n), integer(n), double((n +
            5) * (p + 4)), double(n), as.integer(nsol), as.integer(ndsol),
        sol = double((p + 3) * nsol), dsol = double(n * ndsol),
        lsol = as.integer(0), h = integer(p * nsol), qn = as.double(qn),
        cutoff = as.double(cutoff), ci = double(4 * p), tnmat = double(4 *
            p), as.double(big), as.logical(lci1))
    if (z$flag != 0)
        warning(switch(z$flag, "Solution may be nonunique", "Premature end - possible conditioning problem in x"))
    coef <- z$coef
    coef
}

### This is a stripped down version of rq that
### only estimates betas without bells and whistles uses "fn" method
shortrq.fit.fnb <- function (x, y, tau = 0.5, beta = 0.99995, eps = 1e-06)
{
    n <- length(y)
    p <- ncol(x)
    if (n != nrow(x))
        stop("x and y don't match n")
    rhs <- (1 - tau) * apply(x, 2, sum)
    d   <- rep(1,n)
    u   <- rep(1,n)
    wn <- rep(0,10*n)
    wn[1:n] <- (1-tau) #initial value of dual solution
    z <- .Fortran("rqfnb", as.integer(n), as.integer(p), a = as.double(t(as.matrix(x))),
        c = as.double(-y), rhs = as.double(rhs), d = as.double(d),as.double(u),
        beta = as.double(beta), eps = as.double(eps),
        wn = as.double(wn), wp = double((p + 3) * p),
        it.count = integer(3), info = integer(1))
    if (z$info != 0)
        stop(paste("Error info = ", z$info, "in stepy: singular design"))
    coefficients <- -z$wp[1:p]
    names(coefficients) <- dimnames(x)[[2]]
    return ( coefficients )
}




QICD <- function(y, x, tau=.5,lambda=NULL, weights=NULL, beta_initial=NULL, intercept=TRUE,  penalty="SCAD",a=3.7, 
                converge_criteria=1e-06, maxin=100,maxout=20, method="br",...)
#x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
#y: response variable, length n vector
#lambda is the tuning parameter (numeric value > 0)
#tau is the quantile value
#beta_intial: initial value for x-covariates, the default value is NULL (lasso estimates will be used)
#intercept is a logical value,should intercept(s) be fitted (default=TRUE) or set to zero(FALSE)
#a is scale parameter, the default value is 3.7 for SCAD
#penalty is the name of nonconvex penalty function, could be "SCAD" and "MCP", the default
#value is "SCAD"
#converge_criteria is the convergence threshold for coordinate descent and majorization minimization step.
#Default value is 1E-6
#maxin: maximum number of iterations for inside coordinate descent,default value is 100
#maxout: maximum number of iterations for outside MM step,default value is 20
{
#######################################################

   if(is.null(dim(x))){                                                                                    
      stop('x needs to be a matrix with more than 1 column')
   }

  # Make sure all lambdas are positive
  if( lambda <= 0){
    stop("lambda must be positive.")
  }

  # Make sure we use br or fn method ("?rq.fit.br" or "?rq.fit.fnc")
  if( method != "br" & method != "fn"){
    stop("Incorrect method.  Choose br or fn")
  }

  # Get some basic info
  nyrow  <- as.integer( length(y) )
  nxcol  <- as.integer( ncol(x) )
  nxzcol <- as.integer( nxcol + intercept ) # Total number of betas to estimate
  p      <- as.integer( nxcol ) # Number of covariates (not intercept)
  nonpen <- as.integer( intercept ) # Number of nonpenalized betas

   if(dim(x)[1] != nyrow){
      stop('length of y and rows of x do not match')
    }

  if( is.null(beta_initial) ){
    # Put intercept as last value
    beta_initial <- rq.lasso.fit(x=x,y=y,tau=tau,lambda=lambda,weights=weights,intercept=intercept,...)$coefficients[c((1:p)+intercept, nonpen)]
  } else {
    if(intercept){ beta_initial <- c(beta_initial,0) }
    if( length(beta_initial) != nxzcol ){
      stop("beta_initial does not have length ncol(x)")
	}
  }



  # Check penalty functions
  if(penalty == "SCAD"){
      pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
      pentype <- as.integer(1)
  } else{
    stop("wrong penalty function")
  }


  # Name the coefficients
  rowname <- paste("x", 1:nxcol, sep="")
  if(intercept){
    rowname <- c("(Intercept)", rowname )
  }


  # Create design matrix (Penalized, Intercept, Nonpenalized)
  xdes <- x
  if( intercept ){
      xdes <- cbind(xdes,1)
  }

  if( !is.null(weights) ){
      xdes <- xdes * weights
      y <- y * weights
  }

  lambda.new <- lambda/nyrow
  beta_final <- matrix( NA, ncol=1, nrow=nxzcol )
  df <- rep( NA, length(lambda) ) # Number of nonzero coefficients

  # Choose function for finding new beta and do some housekeeping
  y         <- as.double(y)
  xdescol   <- as.integer( ncol(xdes) )
  tau       <- as.double(tau)
  int       <- as.integer(intercept)
  a         <- as.double(a)
  thresh    <- as.double(converge_criteria)
  maxin     <- as.integer(maxin)

    beta1 <- QICDx( y=y, nyrow=nyrow, x=xdes, nxcol=xdescol, beta=beta_initial,
              tau=tau, intercept=int, pentype=pentype, lambda=as.double(lambda.new), a=a,
              thresh=thresh, maxin=maxin, maxout=maxout, method=method )
    #final beta coefficients, Order of coefficients is linear,intercept

  residuals <-  y - xdes%*%beta1
  rho <- sum(check(residuals, tau=tau ))

  if(intercept){ 
    coefficients <- c(beta1[ c( nxzcol, 1:nxcol )])
    PenRho <- rho + sum(scad( abs(coefficients[-1]), lambda=lambda, a=a ))
  } else {
    coefficients <- c(beta1)
    PenRho <- rho + sum(scad( abs(coefficients), lambda=lambda, a=a ))
  }
  names(coefficients) <- rowname

  sub_fit <- list(coefficients = coefficients,
                 PenRho = PenRho,
                 residuals = residuals,
                 rho = rho,
                 tau = tau,
                 n = nyrow,
				 intercept=intercept,
                 penalty = penalty)

  class(sub_fit) <-  c("rq.pen", "rqNC")

  return( sub_fit )
}



# QICD.nonpen <- function(x, z, y, tau=.5,lambda=NULL, weights=NULL, beta_initial=NULL, intercept=TRUE,  penalty="SCAD",a=3.7, 
                # converge_criteria=1e-06, maxin=100,maxout=20, method="br",...)
# #x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
# #y: response variable, length n vector
# #z: spline basis matrix
# #beta: initial value for x-covariates, the default value is NULL
# #tau is the quantile value
# #lambda is the tuning parameter sequence
# #a is scale parameter, the default value is 3.7 for SCAD
# #penalty is the name of nonconvex penalty function, could be scad and mcp, the default
# #value is scad
# #intercept is a logical value,should intercept(s) be fitted (default=TRUE) or set to zero(FALSE)
# #thresh is the convergence threshold for coordinate descent and majorization minimization step.
# #Default value is 1E-6
# #maxin: maximum number of iterations for inside coordinate descent,default value is 100
# #maxout: maximum number of iterations for outside MM step,default value is 20
# {
# #######################################################
   # if(is.null(dim(x))){                                                                                    
      # stop('x needs to be a matrix with more than 1 column')
   # }
   # if(is.null(dim(z))){                                                                                    
      # stop('z needs to be a matrix with 1 or more columns')
   # }

  # # Make sure all lambdas are positive
  # if( lambda <= 0){
    # stop("lambda must be positive.")
  # }

  # # Make sure we use br or fn method ("?rq.fit.br" or "?rq.fit.fnc")
  # if( method != "br" & method != "fn"){
    # stop("Incorrect method.  Choose br or fn")
  # }

  # # Get some basic info
  # nyrow  <- as.integer( length(y) )
  # nzcol  <- as.integer( ncol(z) )
  # nxcol  <- as.integer( ncol(x) )
  # nxzcol <- as.integer( nxcol + nzcol + intercept ) # Total number of betas to estimate
  # p      <- as.integer( nxcol + nzcol ) # Number of covariates (not intercept)
  # nonpen <- as.integer( nzcol + intercept ) # Number of nonpenalized betas

   # if(dim(x)[1] != nyrow){
      # stop('length of y and rows of x do not match')
    # }
   # if(dim(z)[1] != nyrow){
      # stop('length of y and rows of z do not match')
    # }

  # if( is.null(beta_initial) ){
    # # Put intercept after linear coefficients, but before nonlinear coefficients
    # beta_initial <- rq.lasso.fit(x=x,y=y,tau=tau,lambda=lambda,weights=weights,intercept=intercept,...)$coefficients
    # if(intercept) { beta_initial <-  beta_initial[ c( (1:nxcol)+intercept, 1 ) ] }
    # beta_initial <- c( beta_initial,  rep(0, nzcol) )
  # } else {
    # if(intercept){ beta_initial <- c(beta_initial[1:nxcol], 0, beta_initial[-(1:nxcol)]) }
    # if( length(beta_initial) != nxzcol ){
      # stop("beta_initial does not have length ncol(x)+ncol(z)")
    # }
  # }


  # # Check penalty functions
  # if(penalty == "SCAD"){
      # pentype <- as.integer(0)  
  # } else if (penalty == "MCP"){
      # pentype <- as.integer(1)
  # } else{
    # stop("wrong penalty function")
  # }


  # # Name the coefficients
  # rowname <- c( paste("x", 1:nxcol, sep=""), paste("z", 1:nzcol, sep="") )
  # if(intercept){
    # rowname <- c("(Intercept)", rowname )
  # }


 # # Create design matrix (Penalized, Intercept, Nonpenalized)
  # xdes <- x
  # zdes <- z
  # if( intercept ){
    # zdes <- cbind(1,zdes)
  # }
  # if( !is.null(weights) ){
      # xdes <- xdes * weights
      # zdes <- zdes * weights
      # y <- y * weights
  # }  


  # lambda.new <- lambda/nyrow
  # beta_final <- matrix( NA, ncol=1, nrow=nxzcol )
  # df <- rep( NA, length(lambda) ) # Number of nonzero coefficients

  # # Choose function for finding new beta and do some housekeeping
  # y         <- as.double(y)
  # xdescol   <- as.integer( ncol(xdes) )
  # zdescol   <- as.integer( ncol(zdes) )
  # tau       <- as.double(tau)
  # int       <- as.integer(intercept)
  # a         <- as.double(a)
  # thresh    <- as.double(converge_criteria)
  # maxin     <- as.integer(maxin)

    # beta1 <- QICDz( y=y, nyrow=nyrow, x=xdes, nxcol=xdescol, z=zdes, nzcol=zdescol, beta=beta_initial,
              # tau=tau, pentype=pentype, lambda=as.double(lambda.new), a=a,
              # thresh=thresh, maxin=maxin, maxout=maxout, method=method )
    # #final beta coefficients, Order of coefficients is linear,intercept, nonlinear

  # residuals <-  y - cbind(xdes,zdes) %*% beta1
  # rho <- sum(check(residuals, tau=tau ))

  # if(intercept){ 
    # coefficients <- c(beta1[ c( nxcol+1, 1:nxcol, (nxcol+2):nxzcol )])
    # PenRho <- rho + sum(scad( abs(coefficients[ 1+(1:nxcol) ]), lambda=lambda, a=a ))
  # } else {
    # coefficients <- c(beta1)
    # PenRho <- rho + sum(scad( abs(coefficients[ 1:nxcol ]), lambda=lambda, a=a ))
  # }
  # names(coefficients) <- rowname

  # sub_fit <- list(coefficients = coefficients,
                 # PenRho = PenRho,
                 # residuals = residuals,
                 # rho = rho,
                 # tau = tau,
                 # n = nyrow,
				 # intercept = intercept,
                 # penalty = penalty)

  # class(sub_fit) <-  c("rq.pen", "rqNC")

  # return( sub_fit )
# }

QICD.nonpen <- function(x, z, y, tau=.5,lambda=NULL, weights=NULL, beta_initial=NULL, intercept=TRUE,  penalty="SCAD",a=3.7, 
                converge_criteria=1e-06, maxin=100,maxout=20, method="br",...)
#x: input nxp matrix, of dimension nobs x nvars; each row is an observation vector. 
#y: response variable, length n vector
#z: spline basis matrix
#beta: initial value for x-covariates, the default value is NULL
#tau is the quantile value
#lambda is the tuning parameter sequence
#a is scale parameter, the default value is 3.7 for SCAD
#penalty is the name of nonconvex penalty function, could be scad and mcp, the default
#value is scad
#intercept is a logical value,should intercept(s) be fitted (default=TRUE) or set to zero(FALSE)
#thresh is the convergence threshold for coordinate descent and majorization minimization step.
#Default value is 1E-6
#maxin: maximum number of iterations for inside coordinate descent,default value is 100
#maxout: maximum number of iterations for outside MM step,default value is 20
{
#######################################################
   if(is.null(dim(x))){                                                                                    
      stop('x needs to be a matrix with more than 1 column')
   }

  # Make sure all lambdas are positive
  if( lambda <= 0){
    stop("lambda must be positive.")
  }

  # Make sure we use br or fn method ("?rq.fit.br" or "?rq.fit.fnc")
  if( method != "br" & method != "fn"){
    stop("Incorrect method.  Choose br or fn")
  }

  # Get some basic info
  nyrow  <- as.integer( length(y) )
  if(is.null(dim(z))){
	nzcol <- 1
	nzrow <- length(z)
  } else{
	nzcol  <- as.integer( ncol(z) )
	nzrow <- nrow(z)
  }
  nxcol  <- as.integer( ncol(x) )
  nxzcol <- as.integer( nxcol + nzcol + intercept ) # Total number of betas to estimate
  p      <- as.integer( nxcol + nzcol ) # Number of covariates (not intercept)
  nonpen <- as.integer( nzcol + intercept ) # Number of nonpenalized betas

   if(dim(x)[1] != nyrow){
      stop('length of y and rows of x do not match')
    }
   if(nzrow != nyrow){
      stop('length of y and rows of z do not match')
    }

  if( is.null(beta_initial) ){
    # Put intercept after linear coefficients, but before nonlinear coefficients
    beta_initial <- rq.lasso.fit(x=x,y=y,tau=tau,lambda=lambda,weights=weights,intercept=intercept,...)$coefficients
    if(intercept) { beta_initial <-  beta_initial[ c( (1:nxcol)+intercept, 1 ) ] }
    beta_initial <- c( beta_initial,  rep(0, nzcol) )
  } else {
    if(intercept){ beta_initial <- c(beta_initial[1:nxcol], 0, beta_initial[-(1:nxcol)]) }
    if( length(beta_initial) != nxzcol ){
      stop("beta_initial does not have length ncol(x)+ncol(z)")
    }
  }


  # Check penalty functions
  if(penalty == "SCAD"){
      pentype <- as.integer(0)  
  } else if (penalty == "MCP"){
      pentype <- as.integer(1)
  } else{
    stop("wrong penalty function")
  }


  # Name the coefficients
  rowname <- c( paste("x", 1:nxcol, sep=""), paste("z", 1:nzcol, sep="") )
  if(intercept){
    rowname <- c("(Intercept)", rowname )
  }


 # Create design matrix (Penalized, Intercept, Nonpenalized)
  xdes <- x
  zdes <- z
  if( intercept ){
    zdes <- cbind(1,zdes)
  }
  if( !is.null(weights) ){
      xdes <- xdes * weights
      zdes <- zdes * weights
      y <- y * weights
  }  


  lambda.new <- lambda/nyrow
  beta_final <- matrix( NA, ncol=1, nrow=nxzcol )
  df <- rep( NA, length(lambda) ) # Number of nonzero coefficients

  # Choose function for finding new beta and do some housekeeping
  y         <- as.double(y)
  xdescol   <- as.integer( ncol(xdes) )
  zdescol   <- as.integer( ncol(zdes) )
  tau       <- as.double(tau)
  int       <- as.integer(intercept)
  a         <- as.double(a)
  thresh    <- as.double(converge_criteria)
  maxin     <- as.integer(maxin)

    beta1 <- QICDz( y=y, nyrow=nyrow, x=xdes, nxcol=xdescol, z=zdes, nzcol=zdescol, beta=beta_initial,
              tau=tau, pentype=pentype, lambda=as.double(lambda.new), a=a,
              thresh=thresh, maxin=maxin, maxout=maxout, method=method )
    #final beta coefficients, Order of coefficients is linear,intercept, nonlinear

  residuals <-  y - cbind(xdes,zdes) %*% beta1
  rho <- sum(check(residuals, tau=tau ))

  if(intercept){ 
    coefficients <- c(beta1[ c( nxcol+1, 1:nxcol, (nxcol+2):nxzcol )])
    PenRho <- rho + sum(scad( abs(coefficients[ 1+(1:nxcol) ]), lambda=lambda, a=a ))
  } else {
    coefficients <- c(beta1)
    PenRho <- rho + sum(scad( abs(coefficients[ 1:nxcol ]), lambda=lambda, a=a ))
  }
  names(coefficients) <- rowname

  sub_fit <- list(coefficients = coefficients,
                 PenRho = PenRho,
                 residuals = residuals,
                 rho = rho,
                 tau = tau,
                 n = nyrow,
				 intercept = intercept,
                 penalty = penalty)

  class(sub_fit) <-  c("rq.pen", "rqNC")

  return( sub_fit )
}


QICDx <- function( y, nyrow, x, nxcol, beta, tau, intercept, pentype, lambda, a,
                    thresh, maxin, maxout, ... ){

  beta1 <- beta
  i=0
  distance <- thresh+1

  while( (i < maxout) & (distance >= thresh) ){

    beta0 <- beta1
    xbeta1 <- x%*%beta1

    out <- .C("QCD", y, as.double(x), as.double(beta0), as.double(beta1), as.double(xbeta1),
                     nyrow, nxcol, intercept, tau, lambda, a, pentype, thresh, maxin)
    beta1 <- out[[4]]
    i <- i+1
    distance <- sqrt( sum((beta1 - beta0)^2) )
  }

  if(i == maxout & distance > thresh){
      warning(paste("did not converge after ", maxout, " iterations", sep=""))
    }

  return( beta1 )
}


QICDz <- function( y, nyrow, x, nxcol, z, nzcol, beta, tau, intercept, pentype, lambda, a,
                    thresh, maxin, maxout, method="br", ... ){

  if( method=="br" ) zcoef.func <- shortrq.fit.br
  if( method=="fn" ) zcoef.func <- shortrq.fit.fnb

  beta1 <- beta[1:nxcol]
  betaz <- beta[-(1:nxcol)]
  i=0
  distance <- thresh+1
  xbeta1 <- x%*%beta1

  while( (i < maxout) & (distance >= thresh) ){
    iter <- 0
    distance.inner <- thresh + 1
    beta0 <- beta1
    betaz0 <- betaz

    while( (iter < maxin) & (distance.inner >= thresh) ){

      beta0in <- beta1
      betazin <- betaz
      y.aug <- y - z%*%betazin

      out <- .C("QCD", as.double(y.aug), as.double(x), as.double(beta0), as.double(beta1), as.double(xbeta1),
                     nyrow, nxcol, as.integer(0), as.double(tau),
                     as.double(lambda), as.double(a), as.integer(pentype),
                     as.double(thresh), as.integer(1))
      beta1 <- out[[4]]

      xbeta1 <- x%*%beta1
      y.aug <- y - xbeta1

      betaz <- zcoef.func(x=z, y=y.aug, tau=tau)
      # coef(rq(y.aug ~ z+0, tau=tau, method=method))

      distance.inner <- sqrt(sum( (beta1 - beta0in)^2 )) +
                        mean( ( z%*%(betaz-betazin) )^2 )

      iter <- iter+1

    }

    distance <- sqrt(sum( (beta1 - beta0)^2 )) +
                mean( ( z%*%(betaz-betaz0) )^2 )

    i <- i+1
  }

  if(i == maxout & distance > thresh){
      warning(paste("did not converge after ", maxout, " iterations", sep=""))
    }

  return( c(beta1, betaz) )
}
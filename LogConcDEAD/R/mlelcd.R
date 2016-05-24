## Actually does the computation of the MLE
## and cleans up the output for feeding in to other functions


'mlelcd' <- function(x,
                     w=rep(1/nrow(x),nrow(x)),
                     y=initialy(x),
                     verbose=-1,
                     alpha=5,
                     c=1,
                     sigmatol=10^-8,
                     integraltol=10^-4,
                     ytol=10^-4,
                     Jtol=0.001,
                     chtol = 10^-6 ) {

  if (is.matrix(x)==FALSE) {
    if(is.numeric(x)) {
      x <- matrix(x, ncol=1)
    }
    else {
      stop("x must be a numeric vector or matrix")
    }
  }
  if(length(w)!=nrow(x)) stop("there must be one w_i for each observation")

  if(sum(w <= 0)) stop("all weights must be strictly positive")
  
  ## Correct the weights if necessary
  if( abs( sum( w ) - 1 ) > 10^-5 ) {
    w <- w/sum(w)
    warning("weights have been renormalized to sum to 1")
  }

  ## If we have a 1-d x, things are simpler
  if(ncol(x)==1) {
    outerpoints <- c(which.min(x),which.max(x))
  } else {
    chull <- convhullnew(x)
    outerpoints <- unique(c(chull))
  }
  
  innerpoints <- (1:nrow(x))[-outerpoints]
  nouter <- length(outerpoints)
  lcdsort <- c(outerpoints,innerpoints)
  
  ##Make the initial vector
  opts <- rep(0,11)
  opts[1] <- as.double(-c) #-for minimization; c for initial step length
  opts[2] <- ytol #distance
  opts[3] <- sigmatol #function
  opts[4] <- 15000 #maximum number of iterations
  opts[5] <- as.double(verbose) #display control
  opts[6] <- integraltol #integral
  opts[7] <- as.double(alpha) #dilation factor
  opts[8] <- 1*10^(-11) #for numerical gradient approx


  ## Set up the qhull options
  chopts <- paste( "Qt", sep="" )
  
  out <- .C( "logconestw",
            yvalue = as.double( y [ lcdsort ] ),
            as.double( x [ lcdsort, ] ),
            as.integer( ncol( x ) ),
            as.integer( nrow( x ) ),
            as.double( w [ lcdsort ] ),
            options = as.double( opts ),
            minvalue = double( 1 ),
            Jtol = as.double( Jtol ),
            chopts = as.character( chopts ),
            nouter = as.integer( nouter ),
            PACKAGE="LogConcDEAD" )

  y[ lcdsort ] <- out$yvalue
  r <- getinfolcd( x, y, w, chtol, out$minvalue, out$options[9:11])
  return( r )
}


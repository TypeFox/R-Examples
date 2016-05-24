censored <- function(x, 
                     y, 
                     a, 
                     delta, 
                     propen, 
                     phi, 
                     logY = TRUE,
                     intercept = TRUE) {

  #------------------------------------------------------------------#
  # x should be data.matrix internally                               #
  #------------------------------------------------------------------#
  if( is(x,"data.frame") ) x <- data.matrix(x)
  if( is(colnames(x),"NULL") ) {
    colnames(x) <- paste("x",1:ncol(x),sep="")
  }

  #------------------------------------------------------------------#
  # response should be vector internally.                            #
  #------------------------------------------------------------------#
  if( is(y,"data.frame") ) y <- data.matrix(y)
  if( is(y,"matrix") ) y <- y[,1L]

  #------------------------------------------------------------------#
  # treatment should be vector internally.                           #
  #------------------------------------------------------------------#
  if( is(a,"data.frame") ) a <- data.matrix(a)
  if( is(a,"matrix") ) a <- a[,1L]
  #------------------------------------------------------------------#
  # Treatments must be integers.                                     #
  #------------------------------------------------------------------#
  if( is(a, "numeric") ) {
    if( !all.equal(a, round(a,0L)) ) {
      stop("treatment variable must be integer valued.")
    }
    a <- as.integer(round(a,0L))
  } else if( !is(a,"integer") ) {
    stop("treatment variable must be integer valued.")
  }
  txOpts <- unique(a)

  #------------------------------------------------------------------#
  # delta should be vector internally.                               #
  #------------------------------------------------------------------#
  if( is(delta,"data.frame") ) delta <- data.matrix(delta)
  if( is(delta,"matrix") ) delta <- delta[,1L]
  #------------------------------------------------------------------#
  # delta values must be integers                                    #
  #------------------------------------------------------------------#
  if( is(delta, "numeric") ) {
    if( !all.equal(delta, round(delta,0L)) ) {
      stop("delta variable must be integer valued.")
    }
    delta <- as.integer(round(delta,0L))
  } else if( !is(delta,"integer") ) {
    stop("delta variable must be integer valued.")
  }

  n <- nrow(x)

  #------------------------------------------------------------------#
  # propen should be matrix internally.                              #
  #------------------------------------------------------------------#
  if( is(propen, "data.frame") ) propen <- data.matrix(propen)
  if( !is(propen, "matrix") ) {
    propen <- matrix(data = propen, 
                     nrow = n,  
                     ncol = length(propen),  
                     byrow = TRUE)
  }
  if( nrow(propen) == 1L ) {
    propen <- matrix(data = propen, 
                     nrow = n,  
                     ncol = length(propen),  
                     byrow = TRUE)
  }
  #------------------------------------------------------------------#
  # Verify that enough propensity scores are provided.               #
  #------------------------------------------------------------------#
  if( ncol(propen) != length(txOpts) ) {
    if( ncol(propen) != {length(txOpts)-1L} ) {
      stop("Propensity scores must be given for each treatment option.")
    }
    propen <- cbind(1.0-rowSums(propen), propen)
  }
  if( nrow(propen) != n ) {
    stop("Row dimension of propen incorrect.")
  }

  #------------------------------------------------------------------#
  # convert phi to lowercase                                         #
  #------------------------------------------------------------------#
  phi <- tolower(phi)

  #------------------------------------------------------------------#
  # Calculate weights                                                #
  #------------------------------------------------------------------#
  wgt <- weights(y = y, delta = delta)

  #------------------------------------------------------------------#
  # If requested, convert y to log(y)                                #
  #------------------------------------------------------------------#
  if( logY ) y <- log(y)

  #------------------------------------------------------------------#
  # Estimate parameters and optimal treatment                        #
  #------------------------------------------------------------------#
  if( phi == "c" ) {
    case1 <- constantCase(x = x, 
                          y = y,  
                          trt = a,  
                          propen = propen,
                          wgt = wgt,
                          intercept = intercept)
  } else if( phi == "l" ) {
    case1 <- linearCase(x = x, 
                        y = y,  
                        trt = a,  
                        propen = propen,
                        wgt = wgt,
                        intercept = intercept)
  } else {
    stop("phi not recognized. Must be one of {'c','l'}.")
  }

  return(case1)

}

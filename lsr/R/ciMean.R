# file:    ciMean.R 
# author:  Dan Navarro
# contact: daniel.navarro@adelaide.edu.au
# changed: 13 November 2013

# ciMean() computes a confidence interval around the sample mean, under the 
# usual assumption of normality.
ciMean <- function(x, conf = .95, na.rm = FALSE) {
  
  ############  check input  ############  
  
  # check x
  if( missing(x) ) { stop( 'argument "x" is missing, with no default')}
  if( !is( x,"integer") & !is( x,"numeric") & !is(x,"matrix") & !is(x,"data.frame") ) {
    stop( '"x" must be numeric, matrix or data frame' )
  }
  
  # check conf
  if( !is( conf,"numeric" )) { stop( '"conf" must be numeric' )}
  if( length(conf) !=1 ) { stop('"conf" must be a single number') }
  if( conf<0 | conf>1 ) {stop( '"conf" must be between 0 and 1')}
  
  # check na.rm
  if( !is( na.rm, "logical" )) { stop( '"na.rm" must be logical') }
  if( length( na.rm ) !=1 ) { stop('"na.rm" must be of length 1') }
  if( !(na.rm %in% c(TRUE,FALSE)) ) {stop( '"na.rm" must be TRUE or FALSE')} # no NA!
  
  
  ############  function for a single ci  ############ 

  getNames <- function( quantiles ) {
    paste(100*quantiles,'%',sep="")
  }
  
  computeCI <- function( x, conf, na.rm ) {
    
    # remove missing data if requested
    if (na.rm) { x <- x[!is.na(x)] }
    
    # key quantities
    quantiles <- c( (1-conf)/2 , (1+conf)/2 ) # quantiles of t 
    n <- length(x) # sample size
    
    # calculate CI
    if (length(x) < 2 | any(is.na(x))) { 
      CI <- c(NA,NA) # undefined ci
    } else{
      CI <- mean(x) + qt(p = quantiles, df = n-1) * sd(x) / sqrt(n) # normal CI  
      if( sd(x)==0 ) warning( "data have zero variance")
    }
 
    # assign names and return
    names(CI) <- getNames(quantiles)
    return(CI)
    
  }
  
  ############  handle different types of input  ############ 
  
  
  # handle numeric input
  if( is(x,"integer") | is(x,"numeric" ) ) {   
    ci <- computeCI( x, conf, na.rm ) # compute the ci
    out <- matrix( NA, 1, 2) # set up output matrix
    out[1,] <- ci # insert data
    colnames(out) <- names(ci) # column names 
    nn <- match.call()[[2]] # get the x argument from the call
    if( class(nn) == "name" ) { # if there is a variable name...
      rownames(out) <- as.character(nn) # ... add it as a row
    }
    return( out )
  }
  
  # handle matrix input
  if( is(x,"matrix") ) {
    if( mode(x)!="numeric" ) {
      stop("matrix input must be numeric")
    }
    d <- dim(x)
    out <- matrix( NA, nrow=d[2], ncol=2 )
    for( v in 1:d[2] ) {
      ci <- computeCI( x[,v], conf, na.rm )
      out[ v, ] <- ci
    }
    rownames(out) <- colnames(x)
    colnames(out) <- names(ci)
    return( out )
  }
  
  # handle data frame input
  if( is(x,"data.frame" )) {
    d <- dim(x)
    out <- matrix( NA, nrow=d[2], ncol=2 )
    nn <- names(x)
    for( v in 1:d[2] ) {
      if( class(x[[v]] ) %in% c("numeric","integer") ) {
        ci <- computeCI( x[[v]], conf, na.rm )
        out[ v, ] <- ci
      } else {
        nn[v] <- paste0( nn[v],"*" )
      }
    }
    rownames(out) <- nn
    quantiles <- c( (1-conf)/2 , (1+conf)/2 )  
    colnames(out) <- getNames(quantiles)
    return( out )
  }  
 
  # throw error
  stop("hey, how did I get here?")
  
}

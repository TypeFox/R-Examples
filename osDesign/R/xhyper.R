##

	
##
##########################################
## BASED ON A 2-LEVEL EXPOSURE VARIABLE ##
##########################################
##	
dXhyper <- function( theta, data, N1x ){
  value <- .C("dXhyper",
              theta = as.double( theta ),
              Mx = as.integer( c(data$M0,data$M1) ), Ny = as.integer( c(data$N0,data$N1) ),
              N1x = as.integer( N1x ),
              value = as.double( 0 ), PACKAGE="osDesign")
  return( value$value )
}

##	
ddXhyper <- function( theta, data, N1x ){
  value <- .C("ddXhyper",
              theta = as.double( theta ),
              Mx = as.integer( c(data$M0,data$M1) ), Ny = as.integer( c(data$N0,data$N1) ),
              N1x = as.integer( N1x ),
              value = as.double( 0 ), PACKAGE="osDesign")
  return( value$value )
}


##	
dddXhyper <- function( theta, data, N1x ){
  value <- .C("dddXhyper",
              theta = as.double( theta ),
              Mx = as.integer( c(data$M0,data$M1) ), Ny = as.integer( c(data$N0,data$N1) ),
              N1x = as.integer( N1x ),
              value = as.double( 0 ), PACKAGE="osDesign")
  return( value$value )
}


##
modeXhyper <- function( theta, data ){
  value <- .C("modeXhyper",
              theta = as.double( theta ),
              Mx = as.integer( c(data$M0,data$M1) ), Ny = as.integer( c(data$N0,data$N1) ),
              value = as.integer( 0 ), PACKAGE="osDesign")
  return( value$value )
}

##
eXhyper <- function( theta, data ){
  value <- .C("eXhyper",
              theta = as.double( theta ),
              Mx = as.integer( c(data$M0,data$M1) ), Ny = as.integer( c(data$N0,data$N1) ),
              value = as.double( 0 ), PACKAGE="osDesign")
  return( value$value )
}

##
vXhyper <- function( theta, data ){
  value <- .C("vXhyper",
              theta = as.double( theta ),
              Mx = as.integer( c(data$M0,data$M1) ), Ny = as.integer( c(data$N0,data$N1) ),
              value = as.double( 0 ), PACKAGE="osDesign")
  return( value$value )
}

##
rXhyper <- function( theta, data, number = 1 ){
  support <- max(0, data$N1-data$M0):min(data$N1, data$M1)
  pmf <- .C("pmfXhyper",
            theta = as.double( theta ),
            Mx = as.integer( c(data$M0,data$M1) ), Ny = as.integer( c(data$N0,data$N1) ),
            value = as.double( rep(0, length(support)) ), PACKAGE="osDesign")
  cmf <- cumsum( pmf$value )

  value <- rep( 0, number )
  for( i in 1:number ) value[i] <- min( support[cmf > runif(1)] )
 
  return( value )
}

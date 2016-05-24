
# this function is internal, not meant for the user

# R code to process user input, then call
# the corresponding function of the C++ object


.cr <- function(z, CL =0.95, method ="clr", incr =NULL, output ="G") {
  if( is.null(incr) )  incr<- -1  else 
    if( incr <= 0 )  stop("'incr' must be positive")
  method <- toupper(method)
  met <- integer(1)
  if( method=="CLR" )  met <- 1  else  {
    if( method=="AF" )  met <- 2  else
      stop( "'method' must be \"CLR\" or \"AF\"" )
  }
#  if(missing(output) && .Device=="null device")  output <- "T"
  output <- toupper(output)
  if( output=="T" )
    (z$CppObj)$cr3( CL, met, incr )
  else  {
    bounds <- (z$CppObj)$cr4( CL, met, incr, as.integer(FALSE) )
    if( output=="V" )
      return( bounds )
    else  {
      if( output=="G" )  {
        nbd <- nrow(bounds)
        cl <- as( round(100*CL,0), "character" )
        title <- paste( cl, 
          "% conf. region for changepoint by ", method, sep="")
        if( (z$xint && length(z$coef)==4) ||
              (!z$xint && length(z$coef)==3) )  {
#  simple model
          n <- length(z$x1)
          x <- y <- matrix( NA, max(n,nbd), 3 )
          x[1:n,1] <- z$x1
          y[1:n,1] <- z$y
          x[1:nbd,2:3] <- bounds[,1]
          y[1:nbd,2:3] <- bounds[,2:3]
          matplot( x, y,
            type=c('p','l','l'), pch=4, lty='solid', col='black',
            main=title, xlab=z$x1nm, ylab=z$ynm )
        } else {
#  multiple model
          x <- y <- matrix( NA, nbd, 2 )
          x[,1:2] <- bounds[,1]
          y[,1:2] <- bounds[,2:3]
          xnm <- paste( "theta (", z$x1nm, ")", sep="")
          matplot( x, y,
            type=c('l','l'), lty='solid', col='black',
            main=title, xlab=xnm, ylab="alpha" )
        }
      } else
        stop("'output' must be \"G\", \"T\" or \"V\"")
    }
  }
}


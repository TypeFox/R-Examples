########################################
##global functions 
##Date: 2012-07-19
#######################################
#library(GeneSelector)
#library(qvalue)
#library(hash)


##check big matrix
.checkadjmatrix <- function( mat, backingpath = NULL, descriptorfile = NULL ) {

  
  descfile <- NULL
  if( is.big.matrix(mat) ) {
    if( is.null(backingpath) | is.null(descriptorfile) ) {
      stop("Error: backingpath and descriptorfile should be specified for big.matrix")
    }
    descfile <- paste( backingpath, descriptorfile, sep="")
  }
  
  if( is.null(rownames(mat)) ) {
    stop("Error: no rownames for mat")
  }
  
  descfile
}





##get system time for seed and then generate random index
randomSeed <- function() {
  curtime <- format(Sys.time(), "%H:%M:%OS4")
  XXX <- unlist(strsplit(curtime, ":"))
  curtimeidx <- (as.numeric(XXX[1])*3600 + as.numeric(XXX[2])*60 + as.numeric(XXX[3]))*10000
  curtimeidx
}




#Function: get indexes of one vector from another vector
##add 20130326
.getIndex <- function( VecTarget, VecTo ){
  res <- match( VecTo, VecTarget )
  names(res) <- VecTo
  res
}



##coefficient of variance
.cv <- function( vec ) {
  return( sd(vec)/mean(vec) )
}


##expression z-score
.matZScore <- function( mat ) {
  mean <- apply( mat, 1, mean )
  sd <- apply( mat, 1, sd )
  zMat <- sweep( mat, MARGIN=1, mean, '-' )
  zMat <- sweep( zMat, MARGIN=1, sd, '/' )
  zMat
}


##Exponential function
.expFun <- function( x, base = exp(1)) {
  
  if( base == exp(1) ) {
     return( exp(x) )
  }else{
     return( base^x )
  }
}



##plot lines
.plotLines <- function( xylist, xlim= NULL, ylim = NULL, colors, legends, legends.xpos, legends.ypos = NULL, lwd = 2, xlab = "", ylab = "", title = "", ... ) {
  
  call <- match.call()

  if( !is.list(xylist)  ) 
    stop("Error: xylist should be a list. For each component, there are two vectors respectively for x and y.")

  lineNum <- length(xylist)
  if( length(lwd) == 1 ) 
     lwd <- rep( lwd, lineNum )
  if( length(colors) != lineNum ) 
     stop("Error: conflict number between xylist and colors.")
   
  yvalues <- c()
  for( ii in 1:lineNum ) {
     yvalues = c(yvalues, xylist[[ii]]$y )
  }
  if( is.null(ylim) )
    ylim = range(0, range(yvalues) )

  for( ii in 1:lineNum ) {
     curData <- xylist[[ii]]
    
    if( ii == 1 ) {
      plot(curData$x, curData$y, ylim = ylim, lwd = lwd[ii],  xlab=xlab, ylab=ylab, type="l", col= colors[ii], main=title, ...)
    }else {
      lines(curData$x, curData$y, lwd = lwd[ii], col= colors[ii], ...)
    }
  }
  
  legend(x = legends.xpos, y=legends.ypos, legend=legends, lty=rep(1, lineNum), col = colors, lwd = lwd, ... )
}



##
.EDist <- function( v1, v2, sizeNormalized = FALSE ) {

  score <- sqrt( sum( (v1-v2)^2) )
  if( sizeNormalized == TRUE )
    score <- score/sqrt( length(v1) )

  return(score)
}






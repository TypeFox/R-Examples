padMatrix <-
function(X, by.row=TRUE, type="reflect", pad.direction="both", 
                      replaceLinearTrend=FALSE) {
  ##Check X is a matrix##
  if(is.matrix(X)!=TRUE){
    stop("X must be a matrix")
  }
  else if(dim(X)[1]==0 || dim(X)[2]==0){
    stop("X must have entries")
  }
  else if(mode(X)!="numeric"){
    stop("X must be of type numeric")
  }
  else if(anyNA(X)==TRUE){
    stop("X must not have any missing values.  Please impute first.")
  }
  
  ##Check type is reflect, periodic, or mean##
  if( !(type %in% c("reflect","mean","periodic")) ){
    stop("invalid value for type")
  }

  ##Check pad.direction is both, front, or rear##
  if( !(pad.direction %in% c("both","front","rear")) ){
    stop("invalid value for pad.direction") 
  }

  ##Check by.row is logical##
  if (!(by.row %in% c(TRUE, FALSE)) ){
    stop("invalid value for by.row")
  }

  ##Check replaceLinearTrend is logical##
  if( !(replaceLinearTrend %in% c(TRUE, FALSE)) ){
    stop("replaceLinearTrend must be logical.")
  }
  

  ##Transpose X if by.row=FALSE##
  if(by.row==FALSE){
    X <- t(X)
  }

  ##Find the smallest J so that length(X)<=2^J
  J <- 0  #start with length=2^0 and increase in loop below
  while(dim(X)[1]>2^J){
    J <- J+1
  }


  ##Estimate a linear trend for each data series in X
  origRegResid <- matrix( nrow=dim(X)[1], ncol=dim(X)[2])
  origRegParam <- matrix( nrow=2, ncol=dim(X)[2])
  for( j in 1:dim(X)[2] ){
    origReg <- lm(X[ , j]~seq(1, dim(X)[1]) )
    origRegParam[ , j] <- as.vector(origReg$coef)  #estimated intercept, slope consecutively
    origRegResid[ , j] <- as.vector(origReg$residuals)
  }
  
  xPad <- matrix( nrow=2^J, ncol=dim(X)[2])
  ##now, pad the matrix accordingly
  if(dim(X)[1]==2^J){
    xPad <- X
    origSeries <- c(1, dim(X)[1]) 
  }else if(pad.direction=="both"){
    if( ((2^J-dim(X)[1])%%2)==0){
      padLength <- (2^J-dim(X)[1])/2  #how much to pad on each side
      origSeries <- c(padLength+1, padLength+dim(X)[1])
      for(j in 1:dim(X)[2]){
        if(type=="reflect"){
          xPad[ , j] <- c( rev(origRegResid[2:(padLength+1), j]), origRegResid[ , j], 
                           rev(origRegResid[(dim(X)[1]-padLength):(dim(origRegResid)[1]-1), j]) )
        }else if(type=="mean"){
          xPad[ , j] <- c( rep(mean(origRegResid[ , j]), padLength), origRegResid[ , j], 
                           rep(mean(origRegResid[ , j]), padLength) )
        }else if(type=="periodic"){
          xPad[ , j] <- c( origRegResid[(dim(origRegResid)[1]-padLength+1):dim(origRegResid)[1], j], 
                           origRegResid[ , j], origRegResid[1:padLength, j] )
        }
      }
    }else{
      padLength <- floor((2^J-dim(X)[1])/2) #how much to pad at front (back is +1)
      origSeries <- c(padLength+1, padLength+dim(X)[1])
      for(j in 1:dim(X)[2]){
        if(type=="reflect"){
          if(padLength==0){
            xPad[ , j] <- c( origRegResid[ , j], origRegResid[(dim(origRegResid)[1]-1), j] )
          }else{
            xPad[ , j] <- c( rev(origRegResid[2:(padLength+1), j]), origRegResid[ , j], 
                             rev(origRegResid[(dim(origRegResid)[1]-padLength-1):(dim(origRegResid)[1]-1), j]) )
        }
        }else if(type=="mean"){
          xPad[ , j] <- c( rep(mean(origRegResid[ , j]), padLength), origRegResid[ , j], 
                           rep(mean(origRegResid[ , j]), (padLength+1)) )
        }else if(type=="periodic"){
          if(padLength==0){
            xPad[ , j] <- c( origRegResid[ ,j], origRegResid[1, j] )
          }else{
            xPad[ , j] <- c( origRegResid[(dim(origRegResid)[1]-padLength+1):dim(origRegResid)[1], j], 
                             origRegResid[ , j], origRegResid[1:(padLength+1), j])
          }
        }
      }
    }
  }else if(pad.direction=="front"){
    padLength <- (2^J - dim(X)[1])
    origSeries <- c(padLength+1, 2^J)
    for(j in 1:dim(X)[2]){
      if(type=="reflect"){
        xPad[ , j] <- c( rev(origRegResid[(2:(padLength+1)), j]), origRegResid[ , j] )
      }else if(type=="mean"){
        xPad[ , j] <- c( rep(mean(origRegResid[ , j]), padLength), origRegResid[ , j] )
      }else if(type=="periodic"){
        xPad[ , j] <- c( origRegResid[(dim(origRegResid)[1]-padLength+1):dim(origRegResid)[1], j], origRegResid[ , j] )
      }
    }
  }else if(pad.direction=="rear"){
    padLength <- (2^J - dim(X)[1])
    origSeries <- c(1, dim(X)[1])
    for(j in 1:dim(X)[2]){
      if(type=="reflect"){
        xPad[ , j] <- c( origRegResid[ , j], rev(origRegResid[(dim(origRegResid)[1]-padLength):(dim(origRegResid)[1]-1), j]) )
      }else if(type=="mean"){
        xPad[ , j] <- c( origRegResid[ , j], rep(mean(origRegResid[ , j]), padLength) )
      }else if(type=="periodic"){
        xPad[ , j] <- c( origRegResid[ , j], origRegResid[1:padLength, j] )
      }
    }
  }
  
  if(replaceLinearTrend==TRUE){
    xPad <- ( matrix(1, nrow=2^J, ncol=1) %*% origRegParam[1,] + 
              matrix(seq(1, 2^J), nrow=2^J, ncol=1) %*% origRegParam[2, ] + xPad )
  }
  if(by.row==FALSE){
    xPad <- t(xPad)
  }
  return(list(xPad=xPad, origSeriesIndex=origSeries, linearParam=origRegParam, by.row=by.row))
}

padVector <-
function(X, type="reflect", pad.direction="both", replaceLinearTrend=FALSE) {
  ##Check X is a numeric vector##
  if(is.atomic(X)!=TRUE){
    stop("X must be an atomic object")
  }
  else if(length(X)==0){
    stop("X must have length greater than 0")
  }
  else if( !is.null(dim(X)) ){
    stop("X must be a vector")
  }
  else if(mode(X)!="numeric"){
    stop("X must be of type numeric")
  }
  else if(anyNA(X)==TRUE){
    stop("X must not have any missing values.  Please impute first.")
  }
  
  ##Check type is reflect, periodic or mean##
  if( !(type %in% c("reflect","mean","periodic")) ){
    stop("Invalid value for type")
  }

  ##Check pad.direction is both, front, or rear##
  if( !(pad.direction %in% c("both","front","rear")) ){
    stop("Invalid value for pad.direction") 
  }

  ##Check replaceLinearTrend is logical##
  if( !(replaceLinearTrend %in% c(TRUE, FALSE)) ){
    stop("replaceLinearTrend must be logical.")
  }
  
  ##Find the smallest J so that length(x)<=2^J
  J <- 0  #start with length=2^0 and increase in loop below
  while(length(X)>2^J){
    J <- J+1
  }

  ##Estimate the slope and intercept of un-padded data
  origReg <- lm(X~seq(1, length(X)) )
  origRegParam <- as.vector(origReg$coef)  #estimated intercept, slope consecutively
  origRegResid <- as.vector(origReg$residuals) 

  ##now, pad the vector accordingly
  if(length(X)==2^J){
    xPad <- origRegResid
    origSeries <- c(1, length(X))
  }else if(pad.direction=="both"){
    if((2^J-length(X))%%2==0){
      padLength <- (2^J-length(X))/2  #how much to pad on each side
      origSeries <- c(padLength+1, padLength+length(X))
      if(type=="reflect"){
        xPad <- c( rev(origRegResid[2:(padLength+1)]), origRegResid, 
                   rev(origRegResid[(length(origRegResid)-padLength):(length(origRegResid)-1)]) )
      }else if(type=="mean"){
        xPad <- c( rep(mean(origRegResid), padLength), origRegResid, rep(mean(origRegResid), padLength) )
      }else if(type=="periodic"){
        xPad <- c( origRegResid[(length(origRegResid)-padLength+1):length(origRegResid)], 
                   origRegResid, origRegResid[1:padLength])
      }
    }else{
      padLength <- floor((2^J-length(X))/2) #how much to pad at front (back is +1)
      origSeries <- c(padLength+1, padLength+length(X))
      if(type=="reflect"){
        if(padLength==0){
          xPad <- c( origRegResid, origRegResid[(length(origRegResid)-1)] )
        }else{
          xPad <- c( rev(origRegResid[2:(padLength+1)]), origRegResid, 
                     rev(origRegResid[(length(origRegResid)-padLength-1):(length(origRegResid)-1)]) )
        }
      }else if(type=="mean"){
        xPad <- c( rep(mean(origRegResid), padLength), origRegResid, rep(mean(origRegResid), (padLength+1)) )
      }else if(type=="periodic"){
        if(padLength==0){
          xPad <- c( origRegResid, origRegResid[1] )
        }else{
          xPad <- c( origRegResid[(length(origRegResid)-padLength+1):length(origRegResid)], 
                     origRegResid, origRegResid[1:(padLength+1)])
        }
      }
    }
  }else if(pad.direction=="front"){
    padLength <- (2^J - length(X))
    origSeries <- c(padLength+1, 2^J)
    if(type=="reflect"){
      xPad <- c( rev(origRegResid[2:(padLength+1)]), origRegResid )
    }else if(type=="mean"){
      xPad <- c( rep(mean(origRegResid), padLength), origRegResid )
    }else if(type=="periodic"){
      xPad <- c( origRegResid[(length(origRegResid)-padLength+1):length(origRegResid)], origRegResid)
    }
  }else if(pad.direction=="rear"){
    padLength <- (2^J - length(X))
    origSeries <- c(1, length(X))
    if(type=="reflect"){
      xPad <- c( origRegResid, rev(origRegResid[(length(origRegResid)-padLength):(length(origRegResid)-1)]) )
    }else if(type=="mean"){
      xPad <- c( origRegResid, rep(mean(origRegResid), padLength) )
    }else if(type=="periodic"){
      xPad <- c(origRegResid, origRegResid[1:padLength])
    }
  }

  ##Put the estimated trend back in the padded series if specified
  if(replaceLinearTrend==TRUE){
    xPad <- origRegParam[1] + seq(1, 2^J)*origRegParam[2] + xPad
  }

  return(list(xPad=xPad, origSeriesIndex=origSeries, linearParam=origRegParam))
}

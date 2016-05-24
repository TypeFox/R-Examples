smoothTimeSeries <-
function(X, wavFam="DaubLeAsymm", wavFil=8, wavBC="periodic",
                             plotLevels="none", ...){
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
  else if(!(log(length(X), base=2)%%1==0)){
    stop("X must have a length which is a power of 2")
  }
  else if(length(X)<=8){
    stop("X should have a length greater than 8.  This method is not useful otherwise.")
  }

  J <- log(length(X), base=2)   #keep the power of 2 associated with X length 

  ##Check boundary condition##
  if(!(wavBC %in% c("periodic","symmetric")) ){
    stop("wavBC may only take values periodic or symmetric.")
  }

  ##Check wavelet filter and family##
  if(!( (wavFam=="DaubLeAsymm" && wavFil %in% seq(4,10)) || 
        (wavFam=="DaubExPhase" && wavFil %in% seq(1,10)) )){
    stop("wavFam and wavFil combination not allowed in wavethresh.")
  }

  ##Check plot option
  if(!(plotLevels %in% c("none", "all", seq(0, J-1))) ){
    stop("plotLevels is not specified correctly.")
  }
  dots <- list(...)


  
  XWaveDecomp <- wd(X, filter.number=wavFil, family=wavFam, type="wavelet", bc=wavBC)

  ##Start a matrix to keep smoothed series in##
  smoothSeries <- matrix( nrow=2^J, ncol=(J+1) )
  smoothSeries[ , 1] <- X

  ##Names for data frame that we return##
  dataNames <- "origSeries"
  
  holdWave <- XWaveDecomp
  for(j in 1:J){
    coefSeq <- rep(0, 2^(J-j))
    holdWave <- putD(holdWave, level=(J-j), coefSeq)
    smoothSeries[ , (j+1)] <- wr(holdWave)
    dataNames <- c(dataNames, paste("J0plusOne",(J-j), sep='') )
  }

  smoothSeries <- data.frame(smoothSeries)
  colnames(smoothSeries) <- dataNames

  if(plotLevels!="none"){
    if(plotLevels!="all"){
      dots$x <- seq(1, 2^J)
      dots$y <- smoothSeries[,1]
      if(is.null(dots$main)==TRUE){dots$main <- bquote("Observed Data and Smooth " ~ J[0] ~ "+1=" ~ .(plotLevels))}
      if(is.null(dots$xlab)==TRUE){dots$xlab <- ""}
      if(is.null(dots$ylab)==TRUE){dots$ylab <- ""}
      if(is.null(dots$col)==TRUE){dots$col <- "lightgray"}
      if(is.null(dots$type)==TRUE){dots$type <- "l"}
      do.call("plot", dots)
      lines(seq(1, 2^J), smoothSeries[, colnames(smoothSeries)%in%c(paste("J0plusOne",plotLevels,sep=''))],
           lwd=2)
    }else{
      par(mfrow=c(ceiling(sqrt(J)), ceiling(J/ceiling(sqrt(J))) ), mai = c(0.4, 0.4, 0.4, 0.4))
      for(i in 1:J){
        dots$x <- seq(1, 2^J)
        dots$y <- smoothSeries[,1]
        if(is.null(dots$main)==TRUE || dotsMainPop==1){
          dots$main <- bquote("Observed Data and Smooth " ~ J[0] ~ "+1=" ~ .(i-1))
          dotsMainPop <- 1  #to change the plot title each time (dots$main will not be null after first loop)
        }else{ dotsMainPop <- 0}
        if(is.null(dots$xlab)==TRUE){dots$xlab <- ""}
        if(is.null(dots$ylab)==TRUE){dots$ylab <- ""}
        if(is.null(dots$col)==TRUE){dots$col <- "lightgray"}
        if(is.null(dots$type)==TRUE){dots$type <- "l"}
        do.call("plot", dots)
        lines(seq(1, 2^J), smoothSeries[, J+2-i], lwd=2)
      }
    }
  }

  return(smoothSeries)
}

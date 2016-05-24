#' Create functional boxplot using 'bagplot', 'KDE' or 'pointwise' methodology
#'
#' Using an FPCA object create a functional box-plot based on the function scores.
#' The green line corresponds to the functional median, the dark grey area to the area spanned
#' by the curves within the 25th and 75-th percentile and the light gret to the area spanned
#' by the curves within the 2.5th and 97.5-th percentile. 
#'
#' @param fpcaObj An object of class FPCA returned by the function FPCA().
#' @param optns A list of options control parameters specified by \code{list(name=value)}. See `Details'.
#' @param ... Additional arguments for the 'plot' function. 
#' 
#' @details Available control options are 
#' \describe{
#' \item{ifactor}{inflation ifactor for the bag-plot defining the loop of bag-plot or multiplying ifactor the KDE pilot bandwidth matrix. (see ?aplpack::compute.bagplot; ?ks::Hpi respectively; default: 2.58; 2 respectively).}
#' \item{variant}{string defining the method used ('KDE', 'pointwise' or 'bagplot') (default: 'bagplot')}
#' \item{unimodal}{logical specifying if the KDE estimate should be unimodal (default: FALSE, relavant only for variant='KDE')}
#' \item{addIndx}{vector of indeces corresponding to which samples one should overlay (Default: NULL)}
#' \item{k}{integer number of the first k components used for the representation. (default: length(fpcaObj$lambda ))} 
#' }
#'  
#' @examples
#' set.seed(1)
#' n <- 20
#' pts <- seq(0, 1, by=0.05)
#' sampWiener <- Wiener(n, pts)
#' sampWiener <- Sparsify(sampWiener, pts, 10)
#' res <- FPCA(sampWiener$yList, sampWiener$tList, 
#'             list(dataType='Sparse', error=FALSE, kernel='epan', verbose=TRUE))
#' CreateFuncBoxPlot(res, list(addIndx=c(1:3)) )
#' @references
#' \cite{P. J. Rousseeuw, I. Ruts, J. W. Tukey (1999): The bagplot: a bivariate boxplot, The American Statistician, vol. 53, no. 4, 382-387}
#'
#' @export

CreateFuncBoxPlot <- function(fpcaObj, optns = list() , ...){
   
  if(is.null(optns$k)){
    k = length(fpcaObj$lambda)
  } else {
    k = optns$k
  }
  
  if(is.null(optns$addIndx)){
    addIndx = NULL
  } else {
    addIndx = optns$addIndx
  }
  
  if(is.null(optns$ifactor)){
    ifactor = NULL
  } else {
    ifactor = optns$ifactor
  }
  
  if(is.null(optns$outlierList)){
    outlierList = NULL
  } else {
    outlierList = optns$outlierList
  }
  
  if(is.null(optns$unimodal)){
    unimodal = FALSE
  } else {
    unimodal = optns$unimodal
  }
  
  if(is.null(optns$variant)){
    variant = 'bagplot'
  } else {
    variant = optns$variant
  }
  if ( !is.null(unimodal) && !is.logical(unimodal) ){
    stop("The variable 'unimodal' must be logical.")
  }  
  args1 <- list( xlab='s', ylab='')  
  inargs <- list(...)
  args1[names(inargs)] <- inargs
  if( is.na( any(match( variant, c('pointwise', 'bagplot', 'KDE') )) ) ){
    stop("This plotting utility function can only implement a 'KDE', 'bagplot' or 'pointwise' mapping.")
    return(NULL)
  }  
  if ( variant == 'bagplot' && !is.element('aplpack', installed.packages()[,1])){
    warning('Cannot use bagplot because aplpack::compute.bagplot is unavailable; reverting to point-wise');
    variant = 'pointwise'
  }  
  if ( variant == 'KDE' && !is.element('ks', installed.packages()[,1])){
    warning('Cannot use KDE because ks::kde is unavailable; reverting to point-wise');
    variant = 'pointwise'
  }
  
  fittedCurves <- fitted(fpcaObj, k = k, ...)   
  s <- fpcaObj$workGrid
  N <- nrow(fittedCurves)
  
  do.call(plot, c(list(type='n'), list(x=s), list(y=s), 
                  list(ylim=range(fittedCurves)), args1))
  grid()   
  
  if ( length(fpcaObj$lambda) <2) {
    warning('There is a single component used. We will use a standard boxpot on the FPC scores.');
    bgObj = boxplot(plot=FALSE, fpcaObj$xiEst[,1] )
    fittedCurvesFence = fittedCurves[ (fpcaObj$xiEst > bgObj$stats[1]) & (fpcaObj$xiEst < bgObj$stats[5]),];
    fittedCurvesBag = fittedCurves[ (fpcaObj$xiEst > bgObj$stats[2]) & (fpcaObj$xiEst < bgObj$stats[4]),];
    polygon(x=c(s, rev(s)), y = c(apply(rbind(fittedCurvesFence, fittedCurvesBag),2, min), 
                                  rev(apply(rbind( fittedCurvesFence, fittedCurvesBag),2, max))), col= 'lightgrey',border=0)
    polygon(x=c(s, rev(s)), y = c(apply(fittedCurvesBag,2, min), 
                                  rev(apply(fittedCurvesBag,2,max))), col= 'darkgrey',border=1)  
    lines(x=s, y= apply(fittedCurves,2, mean) , col='red')
  }
  
  if ( length(fpcaObj$lambda)> 1) {
    
    if ( !is.null(ifactor) && (1 >= ifactor) ){
      warning("It is nonsensical for an inflation factor to be <= 1. 'ifactor' set to 1.1.")
      ifactor = 1.1;
    } 
    if ( variant == 'bagplot' ){
      
      if ( is.null (ifactor) ){
        ifactor = 2.58
      }       
      
      bgObj = aplpack::compute.bagplot(x= fpcaObj$xiEst[,1], y= fpcaObj$xiEst[,2], approx.limit=3333, factor = ifactor)     
      fittedCurvesFence = fittedCurves[ is.element( rowSums(fpcaObj$xiEst[,1:2]), rowSums(bgObj$pxy.outer) ),]; 
      fittedCurvesBag = fittedCurves[ is.element( rowSums(fpcaObj$xiEst[,1:2]), rowSums(bgObj$pxy.bag) ),];
      
      Y95 = c(apply(rbind( fittedCurvesFence, fittedCurvesBag),2, min), 
              rev(apply(rbind( fittedCurvesFence, fittedCurvesBag),2, max)))
      Y50 = c(apply(fittedCurvesBag,2, min), rev(apply(fittedCurvesBag,2,max))) 
      medianPoint = which.min( apply( fpcaObj$xiEst[,1:2], 1, function(x) sqrt(sum( (x- bgObj$center)^2)) ) )
    } else if (variant== 'pointwise'){ 
      Y95 = c(apply(fittedCurves,2, quantile, 0.025), 
              rev(apply(fittedCurves,2, quantile, 0.975))) 
      Y50 = c(apply(fittedCurves,2, quantile, 0.25), 
              rev(apply(fittedCurves,2, quantile, 0.75)))  
      mediansCurve = apply(fittedCurves,2, quantile, 0.50);
      medianPoint = NULL
    } else if (variant == 'KDE') {
      if (is.null(ifactor)){
        ifactor = 2
      }
      fhat <- ks::kde(x=fpcaObj$xiEst[,1:2], gridsize = c(400,400), compute.cont = TRUE, 
                      H = ks::Hpi( x=fpcaObj$xiEst[,1:2], binned=TRUE,  pilot="dscalar"  ) *  ifactor) 
      zin = fhat$estimate
      
      if( is.null(unimodal) || unimodal ){
        maxIndex = which( zin == max(zin), arr.ind = TRUE)
        zin = monotoniseMatrix( fhat$estimate, maxIndex[1], maxIndex[2])
      }   
      
      qq = quickNNeval(xin = fhat$eval.points[[1]], yin = fhat$eval.points[[2]], zin = zin, 
                       xout = fpcaObj$xiEst[,1], yout = fpcaObj$xiEst[,2] ) 
      curves0to50= which(qq >=  fhat$cont[50])
      curves50to95 = which(qq >  fhat$cont[95] & qq <=  fhat$cont[50])
      
      fittedCurvesBag = fittedCurves[ c(curves0to50  ),]; 
      fittedCurvesFence = fittedCurves[  c( curves50to95,curves0to50), ];
      Y95 = c(  apply(rbind(fittedCurvesFence, fittedCurvesBag),2, min), 
                rev(apply(rbind(fittedCurvesFence, fittedCurvesBag),2, max)))
      Y50 = c(apply(fittedCurvesBag,2, min), rev(apply(fittedCurvesBag,2,max))) 
      medianPoint = which.min( apply( fpcaObj$xiEst[,1:2], 1, function(x)  sqrt(sum( (x- c(0,0))^2)) ) )
      
    } else  {
      stop('Additional variants are not yet implemented')
    }
    polygon(x=c(s, rev(s)), y = Y95, col= 'lightgrey',border=1)
    polygon(x=c(s, rev(s)), y = Y50, col= 'darkgrey', border=1)  
    #n = length(Y95)*0.5 
    #lines(x=s, y= apply( cbind(Y95[1:n],Y50[1:n]),1, min) , col='black', lwd=1)
    #lines(x=s, y= apply( cbind(rev(Y95[n+(1:n)]), rev(Y50[n+(1:n)])),1, max) , col='black', lwd=1)
    
    if(!is.null(medianPoint)){
      mediansCurve = fittedCurves[medianPoint,]
    }
    
    lines(x=s, col='green', lwd=2, y = mediansCurve)
  } 
  
  yList = fpcaObj$inputData$y
  tList = fpcaObj$inputData$t 
  
  #add sample lines
  if (!is.null(addIndx) && !is.null(yList) && !is.null(tList)  ){
    for (i in 1:length(addIndx) ) {
      lines(x = tList[[addIndx[i]]] , y= yList[[addIndx[i]]], lwd = 1.5, type='o', pch=0)
    } 
  }
}

quickNNeval <- function(xin,yin, zin, xout, yout){
  xindeces = sapply( xout, function(myArg) which.min( abs( xin - myArg) ), simplify = TRUE)
  yindeces = sapply( yout, function(myArg) which.min( abs( yin - myArg) ), simplify = TRUE )
  return( zin[ cbind(xindeces,yindeces)] )
}

monotonise <- function(x, maxIndex = NULL){
  xq = x;
  if (is.null(maxIndex)){
    maxIndex = which.max(x);
  }
  
  if( maxIndex != length(x) ){
    for (i in 1:( length(x) - maxIndex)){
      if( xq[ i + maxIndex] > xq[maxIndex + i - 1] ){
        xq[ i + maxIndex] = xq[maxIndex + i - 1]
      }
    }
  }
  if (maxIndex >= 3){
    for (i in 1:(maxIndex - 2 )){
      if( xq[ - 1 - i + maxIndex] > xq[maxIndex - i] ){
        xq[ - 1- i + maxIndex] = xq[maxIndex - i]
      }
    }
  }
  return(xq)
} 

monotoniseMatrix = function(zin, xmaxind, ymaxind){
  if(is.null(xmaxind) && is.null(ymaxind)){
    maxIndx = which( max(zin) == zin, arr.ind = TRUE)
    xmaxind = maxIndx[1]
    ymaxind = maxIndx[2]
  }
  zq = zin;
  for (j in 1:dim(zin)[2]){
    for (i in 1:dim(zin)[1]){
      if (i == 1 || j == 1 || j == dim(zin)[1] || i == dim(zin)[2]){
        sizeOut = max( abs(xmaxind - i) +1, abs(ymaxind - j) +1 )
        xcoord = round(   ( seq(i, xmaxind , length.out = sizeOut) ) )
        ycoord = round(   ( seq(j, ymaxind , length.out = sizeOut) ) ) 
        zq[ cbind(xcoord,ycoord) ] = monotonise( zq[ cbind(xcoord,ycoord) ]) 
      }
    }
  }
  return(zq)
}

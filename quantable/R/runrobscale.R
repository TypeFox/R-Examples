#' running robust scaling of arefw
#' @export
#' @param arefw - data array to scale
#' @param k - windows
#' @param scale - should also scaling be applied
#' @return list with scaled data runmed used to center the data and runmad the running MAD used for scaling
#' @examples
#' res = c(rnorm(1000,0,1),rnorm(2000,4,3))
#' res2 = runrobscale(res)
#' par(mfrow=c(2,1))
#' plot(res,type="p",pch="x",col=1,cex=0.3)
#' lines(res2$runmed,col=3)
#'
#' y = runFun( res2$scaled, k=51, func=mad )
#' #hist(y)
#' y2 = runFun(res2$scaled,k=51,func=median)
#' plot(res2$scaled,pch="*")
#' lines(y2,col=2,lwd=3)
#' lines(y2+y,col=3,lwd=3)
#' lines(y2-y,col=3,lwd=3)
#'

runrobscale = function(arefw,k=101,scale=TRUE){
  medianref = runmed( arefw , k=k ,endrule="constant")
  # adjust intensities
  # center around 0
  scalefactor = medianref
  mediana1w = arefw - scalefactor
  madref <- NULL
  if(scale){
    madref <- runFun(arefw,k=k)
    mediana1w <- mediana1w / madref
  }
  return(list("scaled"  = mediana1w, "runmed" = medianref, "runmad" = madref))
}
#' running total ion count scaling (TIC)
#' @export
#' @param arefw a series to scale
#' @param k - the smoothing window
#' @return list with fields scaled - contains scaled data and mean - averages of window k
#' @examples
#' res = c(rnorm(1000,3,2),rnorm(2000,8,1))
#' res2 = runTICscale(res)
#' plot(res,type="p",pch=".",col=1,cex=0.5)
#' lines(1:length(res),res2$mean,col=3)
#' points(res2$scaled, pch=".",cex=3,col=2)
#'
#' @seealso correctIntRTv2 for context
runTICscale = function(arefw,k=101){
  madref <- runFun(arefw,k=k,mean)
  mediana1w <- arefw / madref
  return(list("scaled"  = mediana1w, "mean" = madref))
}

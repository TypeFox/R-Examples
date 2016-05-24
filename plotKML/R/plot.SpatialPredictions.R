# Purpose        : Standard plots of a geostatistical model;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ;
# Dev Status     : pre-Alpha
# Note           : ;

plot.SpatialPredictions <- function(x, plot.predictions=TRUE, vertical=FALSE, col=SAGA_pal[[1]], ...){
  if(plot.predictions==TRUE){    
    par(mfrow=c(2,2))
  } else { if(vertical==TRUE){
       par(mfrow=c(3,1))
     } else {
       par(mfrow=c(1,3))    
     }
  }
  par(mar=c(4.5,4.5,.8,.8))
  
  if(plot.predictions==TRUE){
    if(any(names(x@predicted) %in% x@variable)){
      plot(raster(x@predicted[x@variable]), col=col)
    } else {
      plot(raster(x@predicted[1]), col=col)
    }
    points(x@observed, pch="+", cex=.8)
  }
  if(sum(!is.na(x@validation$residual))>2){
    plot(y=x@validation$var1.pred, x=x@validation$observed, pch=19, asp=1, col = "red", xlab='observed (CV)', col.main = rgb(0.99,0.99,0.99), ylab='predicted (CV)', ...)
    abline(a=0, b=1, lwd=2)
  } else {
    plot(y=x@observed@data[,paste(x@variable, "modelFit", sep=".")], x=x@observed@data[,paste(x@variable)], pch=19, asp=1, col = "red", xlab='observed', col.main = rgb(0.99,0.99,0.99), ylab='predicted', ...)
    abline(a=0, b=1, lwd=2)
  }
  vv <- variogram(as.formula("observed~1"), x@validation)
  plot(x=vv$dist, y=vv$gamma, pch="+", col = "green", xlab='distance', cex=1.4, ylab='gamma', ylim = c(0, max(vv$gamma)))
  lines(x=c(0, vv$dist), y=rep(var(x@validation$observed), length(vv$dist)+1), lty=2)
  ## Residuals (only if available):
  vgmmodel = x@vgmModel
  class(vgmmodel) <- c("variogramModel", "data.frame")
  vvres <- variogram(as.formula(paste(paste(x@variable, "residual", sep="."), "~ 1")), x@observed[paste(x@variable, "residual", sep=".")]) # model residuals
  if(sum(!is.na(x@validation$residual))>2){
    vve <- variogram(residual ~ 1, x@validation["residual"]) # validation residuals
    plot(x=vvres$dist, y=vvres$gamma, pch="+", col = "green", xlab='distance', cex=1.4, ylab='gamma', ylim = c(0, max(vv$gamma)))
    points(x=vve$dist, y=vve$gamma, pch="+", col = "red", cex=1.4)
    vline <- variogramLine(vgmmodel, maxdist=max(vvres$dist), n=length(vvres$dist))
    lines(x=vline$dist, y=vline$gamma, col = "green", lwd=2)
  } else {
    plot(x=c(0, vv$dist), y=rep(var(x@validation$observed), length(vv$dist)+1), xlab='distance', ylab='gamma', lty=2, type="l", ylim = c(0, max(vv$gamma)))
  }
}

setMethod("plot", signature("SpatialPredictions"), plot.SpatialPredictions)
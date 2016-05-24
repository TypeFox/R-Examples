# 23.12.2009 version 1, mvdl
qqFitPlot <- function(y, L, title=NA, xlab=NA, ylab=NA, fat=FALSE)
{

   if ( L$method=="Method I" ){
     X <- getOutliersII(y, FLim=c(L$Fmin, L$Fmax), distribution=L$distribution)
     L$residuals <- X$residuals
     L$yMin <- X$yMin
     L$yMax <- X$yMax
     }

   yHat <- switch(L$distribution,
      normal = y-L$residuals,
      exponential = y-L$residuals,
      exp(log(y)-L$residuals))

   lgAxis <- switch(L$distribution,
      normal = "",
      exponential = "",
      "xy")

   vLeft <- head(yHat[y==L$yMin],1)
   vRight <- head(yHat[y==L$yMax],1)

   # default titles
   if ( is.na(title) )
    title <- paste("QQ plot with outliers,", L$method, "\n", 
                  L$distribution,"distribution, R2 =",paste(round(L$R2,4)))
   if ( is.na(xlab) )
    xlab <- "Predicted"
   if ( is.na(ylab) )
    ylab <- "Observed"
   
   # fat cex when asked for
   axisCex <- titleCex <- labCex <- lineWidth <- charCex <- 1
   outlierColor <- "red"
   if ( fat ){
      axisCex <- 1.5
      titleCex <- 2
      labCex <- 1.5
      lineWidth <- 2
      charCex <- 1.5
      outlierColor<-"black"
      }

   myTitle <- list(title, cex=titleCex)
   xLab <- list(xlab, cex=labCex)
   yLab <- list(ylab, cex=labCex)
    
   iInlier <- c(L$iLeft,L$iRight)
   iInlier <- iInlier[!is.na(iInlier)]
   if ( length(iInlier) == 0 )
      iInlier <- -seq(1,length(y))

   plot(yHat, y,
      log=lgAxis,
      main=myTitle,
      xlab=xLab,
      ylab=yLab,
      col="white")

   points(yHat[-iInlier], y[-iInlier],
      col="black",
      cex=charCex)

   abline(0,1,lw=lineWidth)
   abline(v=vLeft,lty=2,lw=lineWidth)
   abline(v=vRight,lty=2,lw=lineWidth)
   points(yHat[L$iLeft],y[L$iLeft], 
      pch=8, 
      cex=charCex,
      col=outlierColor)
   points(yHat[L$iRight], y[L$iRight],
      pch=8,
      cex=charCex,
      col=outlierColor)
      if ( L$method=="Method I" ){
         abline(h=L$limit[1],lty=2,lw=lineWidth)
         abline(h=L$limit[2],lty=2,lw=lineWidth)
      }
}


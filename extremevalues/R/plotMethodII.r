# 23.12.2009 version 1, mvdl
plotMethodII <- function(y, L, title=NA, xlab=NA, ylab=NA, fat=FALSE)
{
   yMax <- L$yMax
   yMin <- L$yMin
   Lplus <- L$limit[2]
   Lmin <- L$limit[1]

   # decide logarithmic x-axis
   lgAxis <- "x"  
   if ( L$distribution=="normal" )
     lgAxis <- ""

    

   # default titles
   if ( is.na(title) )
    title <- paste("Residual plot with outliers,",L$method,"\n", 
      L$distribution,"distribution,", "R2 =",paste(round(L$R2,4)))
   if ( is.na(xlab) )
    xlab <- "Observed value"
   if ( is.na(ylab) )
    ylab <- "Residual"


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
    

   xLim <- log10(c(min(y),max(y)))
   yLim <- c(min(L$residuals),max(L$residuals))
   xRange<- diff(xLim)
   yRange<- diff(yLim)
   plot(y, L$residuals,
      cex=charCex,
      cex.axis=axisCex,
      log=lgAxis, 
      usr=c(xLim,yLim),
      main=myTitle,
      xlab=xLab,
      ylab=yLab) 
   abline(h=0,lw=lineWidth)
   if ( !is.na(L$alphaConf[2]) ){
      rect(yMax,Lplus,10^(xLim[2]+0.039*xRange),yLim[2]+0.039*yRange, 
         col="#DDDDDD",lty=0)
      abline(h=Lplus, lty=2, lw=lineWidth)
      points(y[L$iRight], L$residuals[L$iRight], 
         pch=8, 
         col=outlierColor, 
         cex=charCex)
      }
   if ( !is.na(L$alphaConf[1]) ){
      rect(10^(xLim[1] - 0.039*xRange), yLim[1] - 0.039*yRange, 
        yMin, Lmin, col="#DDDDDD", lty=0)
      abline(h=Lmin, lty=2, lw=lineWidth)
      points(y[L$iLeft], L$residuals[L$iLeft], 
         pch=8,
         col= outlierColor,
         cex=charCex)
      }
   abline(v=L$yMax,lty=2, lw=lineWidth)
   abline(v=L$yMin,lty=2, lw=lineWidth)
   iInlier <- c(L$iLeft,L$iRight)
   iInlier <- iInlier[!is.na(iInlier)]
   points(y[-iInlier],L$residuals[-iInlier],cex=charCex)
}


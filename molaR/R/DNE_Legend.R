#' Make legend for DNE3d plot
#'
#' plotting subfunction
#'
#' @param DNELabels values for the labels
#' @param scaled logical indicating whether the values are scaled
#' @param edgeMask logical indicating whether of not edges are being masked
#' and that information to be included in the legend
#' @param outlierMask logical indicating whether outliers are masked
#' @param logColors logical indicating colors are on log scale
#' @param size legend scaling factor
#'
#' @details This is an internal function which builds a better DNE plot legend
#'
#' The legend will reflect the elements used in the plot. This is an internal
#' function. Users will have little need or call to interact with it.
#' DNE_Legend()

DNE_Legend <- function(DNELabels, scaled = F, edgeMask = F, outlierMask = F,
                              logColors = F, size=1){
  par(ann=F, mar=c(0,0,0,0))
  layout(matrix(1:2,ncol=2), widths = c(0.75, 0.25))
  plot(1,1, type='n', axes=F)
  start <- 0
  end <- 0.575
  lineSize <- 2 * size
  textSize <- 1.75 * size
  rectSize <- size
  if(logColors==TRUE){
    ExpNumbers <- exp(seq(start,end,0.001))
    adjust <- (ExpNumbers-min(ExpNumbers))
    adjust <- adjust/max(adjust)
    adjust <- adjust*end
    colorslist <- hsv(adjust)
  }
  if(logColors==FALSE){
    colorslist <- hsv(seq(start, end, 0.001))
  }
  legend_gradient <- as.raster(matrix(colorslist, ncol=1))
  edgeblack <- '#000000'
  outliergray <- '#505050'
  XPos1 <- 1.15                     #X location of center (center of title, left edge legend text)
  XPos2 <- XPos1-(0.55*rectSize)    #X location of left edge rectangle
  XPos3 <- XPos1-(0.05*rectSize)    #X location of right edge rectangle
  XPos4 <- XPos1-(0.45*rectSize)    #X location of inside of left tic marks
  XPos5 <- XPos1-(0.15*rectSize)    #X location of inside of right tic marks
  YPos0 <- 0.5                      #Y location of center
  YPos1 <- YPos0+(0.3*rectSize)     #Y location of top of rectangle / legend text
  YPos2 <- YPos0-(0.25*rectSize)    #Y location of bottom of rectangle / legend text
  YPos3 <- YPos0+(0.375*rectSize)   #Y location of title
  YPos4 <- YPos0-(0.27*rectSize)    #Y location of top of box below rectangle
  YPos5 <- YPos0-(0.32*rectSize)    #Y location of bottom of box below rectangle
  YPos6 <- YPos0-(0.295*rectSize)   #Y location of text for box below rectangle
  YPos7 <- YPos0-(0.34*rectSize)    #Y location of top of second box below rectangle
  YPos8 <- YPos0-(0.39*rectSize)    #Y location of bottom of second box below rectangle
  YPos9 <- YPos0-(0.365*rectSize)   #Y location of text for second box below rectangle
  plot(c(0,2),c(0,1), type = 'n', axes = F,xlab = '', ylab = '')
  PlotLabels <- DNELabels
  PlotLabels <- format(PlotLabels, scientific=FALSE)
  text(x=XPos1, y = seq(YPos2,YPos1,l=10), labels = PlotLabels, adj=c(0, NA), cex=textSize)
  rasterImage(legend_gradient, XPos2, YPos2, XPos3,YPos1)
  rect(XPos2, YPos2, XPos3, YPos1, lwd=lineSize)
  segments(x0=rep(XPos2, 10), y0=seq(YPos2,YPos1,l=10), x1=rep(XPos4, 10), y1=seq(YPos2,YPos1,l=10), lwd=lineSize)
  segments(x0=rep(XPos5, 10), y0=seq(YPos2,YPos1,l=10), x1=rep(XPos3, 10), y1=seq(YPos2,YPos1,l=10), lwd=lineSize)
  if(scaled==FALSE && logColors==FALSE){
    text(x=XPos1, y=YPos3, labels=c("DNE Value\nPer Face"), cex=textSize)
  }
  if(scaled==TRUE && logColors==FALSE){
    text(x=XPos1, y=YPos3, labels=c("Scaled\nDNE Value\nPer Face"), cex=textSize)
  }
  if(scaled==FALSE && logColors==TRUE){
    text(x=XPos1, y=YPos3, labels=c("Log DNE Value\nPer Face"), cex=textSize)
  }
  if(scaled==TRUE && logColors==TRUE){
    text(x=XPos1, y=YPos3, labels=c("Scaled\nLog DNE Value\nPer Face"), cex=textSize)
  }
  if(edgeMask==TRUE && outlierMask==FALSE){
    rect(XPos2, YPos5, XPos3, YPos4, lwd=lineSize, col=edgeblack, border="black")
    text(x=XPos1, y=YPos6, labels="Edges", cex=textSize, adj=c(0,NA))
  }
  if(edgeMask==FALSE && outlierMask==TRUE){
    rect(XPos2, YPos5, XPos3, YPos4, lwd=lineSize, col=outliergray, border = "black")
    text(x=XPos1, y=YPos6, labels="Outliers",cex=textSize, adj=c(0,NA))
  }
  if(edgeMask==TRUE && outlierMask==TRUE){
    rect(XPos2, YPos5, XPos3, YPos4, lwd=lineSize, col=edgeblack, border="black")
    text(x=XPos1, y=YPos6, labels="Edges", cex=textSize, adj=c(0,NA))
    rect(XPos2, YPos8, XPos3, YPos7, lwd=lineSize, col=outliergray, border="black")
    text(x=XPos1, y=YPos9, labels="Outliers", cex=textSize, adj=c(0,NA))  
  }
}
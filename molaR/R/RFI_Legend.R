#' function for building a legend for RFI
#' 
#' crucial plotting subfunction for RFI3d
#' @param surfCol color for the 3D surface defaults to gray
#' @param footCol color for the 2D footprint defualts to red
#' @param size cex style scaling parameter
#' @param opac sets the value for the opacity of the tooth surface when that is
#' engaged
#' RFI_Legend()

RFI_Legend <- function(surfCol = "gray", footCol = "red", size = 1, opac = 1)
{
  par(ann=F, mar=c(0,0,0,0))
  layout(matrix(1:2,ncol=2), widths = c(0.75, 0.25))
  plot(1,1, type='n', axes=F)
  textSizeFactor <- 1.75 * size
  lineSizeFactor <- 2 * size
  XPos1 <- 1.1                 #X location of center
  XPos2 <- XPos1-(0.45*size)   #X location of left edge legend text
  XPos3 <- XPos1-(0.5*size)    #X location of right edge color box
  XPos4 <- XPos1-(0.8*size)    #X location of left edge color box
  YPos0 <- 0.5                 #Y location of center
  YPos1 <- YPos0+(0.048*size)  #Y location of top of top box
  YPos2 <- YPos0+(0.01*size)   #Y location of bottom of top box
  YPos3 <- YPos0-(0.048*size)  #Y location of bottom of bottom box
  YPos4 <- YPos0-(0.01*size)   #Y location of top of bottom box
  YPos5 <- YPos0+(0.029*size)  #Y location of top box text
  YPos6 <- YPos0-(0.029*size)  #Y location of bottom box text
  plot(c(0,2),c(0,1), type = 'n', axes = F, xlab = '', ylab = '')
  SurfaceColor <- col2rgb(surfCol, alpha=T)
  SurfaceColor[4] <- SurfaceColor[4]*opac
  text(x=XPos2, y=YPos5, labels="3D Surface", cex=textSizeFactor, adj=c(0,NA))
  rect(XPos4, YPos1, XPos3, YPos2, lwd=lineSizeFactor, col=rgb(red=SurfaceColor[1],
                                                         green=SurfaceColor[2],
                                                         blue=SurfaceColor[3],
                                                         alpha=SurfaceColor[4], 
                                                         maxColorValue=255),
       border="black")
  text(x=XPos2, y=YPos6, labels="2D Footprint", cex=textSizeFactor, adj=c(0,NA))
  rect(XPos4, YPos3, XPos3, YPos4, lwd=lineSizeFactor, col=footCol, border="black")
}
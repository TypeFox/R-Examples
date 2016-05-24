#' function for building a legend in OPC plots
#' 
#' crucial graphics subfunction
#' @param binColors number sequence for bins and their colors
#' @param binNumber numeric number of different bins
#' @param maskDiscard logical determining whether faces will be blacked out because
#' they are discarded
#' @param textCol color for the text in the circle legend
#' @param size scaling factor for the legend size
#' @param lineCol color for the lines in the legend
#' OPC_Legend()

OPC_Legend <- function(binColors=c(1:8), binNumber=8, maskDiscard=F, textCol="black",
                              size = 1, lineCol="black"){
  textSizeFactor <- 1.75*size
  lineSizeFactor <- 2*size
  par(ann=F, mar=c(0,0,0,0))
  plot(1,1, type='n', axes=F)
  thetas <- seq(pi, 3*pi, length = binNumber+1)
  labelthetas <- thetas+((2*pi)/(2*binNumber))
  XPos1 <- 1.3                  #X location of center of circle: Min = 0.725, Max = 1.3
  XPos2 <- XPos1+(0.115*size)   #X location of +X label
  XPos3 <- XPos1-(0.115*size)   #X location of -X label
  XPos4 <- XPos1-(0.025*size)   #X location of left of Discarded box
  XPos5 <- XPos1-(0.005*size)   #X location of right of Discarded box
  YPos1 <- 1.0125               #Y location of center of circle: Min = 0.725, Max = 1.3
  YPos2 <- YPos1+(0.15*size)    #Y location of title
  YPos3 <- YPos1+(0.115*size)   #Y location of +Y label
  YPos4 <- YPos1-(0.115*size)   #Y location of -Y label
  YPos5 <- YPos1-(0.132*size)   #Y location of top of Discarded box
  YPos6 <- YPos1-(0.152*size)   #Y location of bottom of Discarded box
  YPos7 <- YPos1-(0.142*size)   #Y location of Discarded text
  radius <- 0.1*size            #Radius size of circle
  start <- pi
  end <- (2*pi/binNumber)+pi
  for(i in 1:binNumber){
    Edge <- seq(from=start, to=end, length=200)
    EdgeX <- cos(Edge)*radius+XPos1
    EdgeY <- sin(Edge)*radius+YPos1
    polygon(x=c(XPos1, EdgeX, XPos1), y=c(YPos1, EdgeY, YPos1), col=binColors[i], border=NA)
    start <- start+(2*pi/binNumber)
    end <- end+(2*pi/binNumber)
  }
  CircEdge <- seq(from=0, to=(2*pi), length=2000)
  CircEdgeX <- cos(CircEdge)*radius+XPos1
  CircEdgeY <- sin(CircEdge)*radius+YPos1
  polygon(x=CircEdgeX, y=CircEdgeY, lwd=lineSizeFactor, border=lineCol)
  segments(x0=XPos1, y0=YPos1, x1=((cos(thetas)*radius)+XPos1),
           y1=((sin(thetas)*radius)+YPos1), lwd=lineSizeFactor, col=lineCol)
  text(x=XPos1, y=YPos2, labels=c("Orientation Bins"), cex=textSizeFactor)
  text(x=XPos1, y=YPos3, labels=c("+Y"), cex=0.75*textSizeFactor)
  text(x=XPos1, y=YPos4, labels=c("-Y"), cex=0.75*textSizeFactor)
  text(x=XPos2, y=YPos1, labels=c("+X"), cex=0.75*textSizeFactor)
  text(x=XPos3, y=YPos1, labels=c("-X"), cex=0.75*textSizeFactor)
  text(x=(cos(labelthetas)*(0.75*radius)+XPos1), y=(sin(labelthetas)*(0.75*radius)+YPos1),
       labels=c(1:binNumber), cex=0.75*textSizeFactor, col=textCol)
  if(maskDiscard==TRUE){
    rect(XPos4, YPos5, XPos5, YPos6, lwd=lineSizeFactor, col="black")
    text(x=XPos1, y=YPos7, labels="Discarded", cex=textSizeFactor, adj=c(0,NA))
  }
}
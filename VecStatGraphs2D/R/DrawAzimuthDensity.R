DrawAzimuthDensity <- function (azimuths, Period = 15, SVGf = 0) 

{

  # Plots a circular graphic when the data azimuths densities are represented by color bands  
  # A vector is drawn representing the mode azimuth
  #
  # Args:
  #   azimuths : azimuths (degrees)
  #   Period: number of moving average terms
  #   SVG: 
  #     0: the plot is showed only un the graphic window 
  #     1: the plot is saved as SVG graphic 
  #
  # Returns:
  #   Plot the graphic and/or save the graphic as SVG
  #

 if (SVGf == 1) {
   SVGfiles <- list.files(path=getwd(), pattern=glob2rx("figure*.SVG"), ignore.case=TRUE)
   nfiles = length(SVGfiles)
   if (nfiles == 0) {
     fileSVGname = "figure1.SVG"
   }
   else {
	 newnum = nfiles + 1
     fileSVGname = paste("figure",newnum,".SVG", sep = "")
   }
 svg(fileSVGname, width= 8, height=8)
 }

  plot.new()
  par(fg = "white", mar = c(1,1,1,1))  

  his = Histogram(azimuths, 1)
  den = filter(his[ ,2], rep(1 / Period, Period), sides = 2, circular = TRUE)
  minden = min(den)  
  maxden = max(den)
  intden = maxden / 15
  
  center_x = 0
  center_y = 0
  d1 = 12
  plot(range(d1, -d1), range(d1, -d1), type = "n", axes = FALSE, ann = FALSE)

  # vertical and horizontal lines and labels
  d2 = d1 - 2
 
	text(center_x, center_y + d2, "0"   , adj = c(0, 0), pos = 3, col = "black", cex = 1.25)
	text(center_x - d2, center_y, "270" , adj = c(0, 0), pos = 2, col = "black", cex = 1.25)
	text(center_x, center_y - d2, "180" , adj = c(0, 0), pos = 1, col = "black", cex = 1.25)
	text(center_x + d2, center_y, "90"  , adj = c(0, 0), pos = 4, col = "black", cex = 1.25)
 
	segments(0, center_y - d2 - 0.25 , center_x, center_y + d2 + 0.25, col = "black", lwd = 2)
	segments(0 + d2 + 0.25, center_y, center_x - d2 - 0.25, center_y, col = "black", lwd = 2)

  # diagonal lines and labels

	x1 = cos(ToRadians(45)) * d2 + 0.25
	y1 = sin(ToRadians(45)) * d2 + 0.25
	
	segments(-x1, -y1, x1, y1, col = "black", lwd = 2) 
	segments(x1, -y1, -x1, y1, col = "black", lwd = 2)
  
    text(x1, y1, "45", adj = c(-0.25,-0.25), col = "black", cex = 1.25)
    text(x1, -y1, "135", adj = c(-0.25,1), col = "black", cex = 1.25)
    text(-x1, y1, "315", adj = c(1.25,-0.25), col = "black", cex = 1.25)
    text(-x1, -y1,"225", adj = c(1.25,1.25), col = "black", cex = 1.25)

  # color segments
   
  ramp <- colorRamp(c("white", "yellow", "orange", "red"))
  color = rgb(ramp(seq(0, 1, length = 15)), maxColorValue = 255)	
  colclas = round(den / intden) 
	
  for (i in 0:360) {    
    radian = ToRadians(90 - i)  
    x1 = cos(radian) * 10
    y1 = sin(radian) * 10
	x2 = cos(radian) * 8
    y2 = sin(radian) * 8
	segments(x1, y1, x2, y2, lwd = 6, lend = 1, col = "GhostWhite")
	segments(x1, y1, x2, y2, lwd = 6, lend = 1, col = color[colclas[i]])
  }     
 
  # scale
  xini = 0.5
  y1 = -2
  y2 = -3
  for (i in 0:15) {
	segments(xini + (i / 3), y1, xini + (i / 3), y2, lwd = 10, lend = 1, col = color[i])
  }
  text(xini + 2.75, y1 + 0.5, "Scale", cex = 1, col = "Black")
  text(xini + 1, y2 - 0.5, round(minden, 4), cex = 1, col = "Black")
  text(xini + 5 , y2 - 0.5, round(maxden, 4), cex = 1, col = "Black")

  # circumferences
  DrawCircle(0, 0, 10, border = "black") 
  DrawCircle(0, 0, 8, border = "black") 

  # mode arrow
  for (i in 1:359) {
    if (round(den[i], 6) == round(maxden, 6)) { 
      clmode = i 
    }
  }
 
  radian = ToRadians(90 - clmode)    
  x1 = cos(radian) * 8
  y1 = sin(radian) * 8
  arrows(0, 0, x1, y1, length = 0.15, lwd = 2, code = 2, col = 2)
  text(cos(radian) * 6, sin(radian) * 6, clmode, col = "Red")

  # Label, mode value, and sample size

  n = NumberOfElements(azimuths)
  labp = paste("Density plot aplying moving averages of period", Period)
  text(-11.7, 12, labels = labp, pos = 4, col = "black", cex = 1)
  labp = paste("Mode azimuth (red arrow) =", clmode)
  text(-11.7, 11.2, labels = labp, pos = 4, col = "black", cex = 1)
  labp = paste("Sample size, n =", n)
  text(-11.7, 10.4, labels = labp, pos = 4, col = "black", cex = 1)

  # save SVG
  if (SVGf == 1) {
   dev.off()
   print(paste("Plot has been saved as SVG graphic file '",fileSVGname,"' in",getwd()))
  }
}

DrawDistribution <- function (azimuths, SVGf = 0) 

{

  # Plots a circular graphic when the data azimuths are represented by stacked points  
  # The mean vector and condfidence interval are plotted if the Von Mises parameter < 0.9
  #
  # Args:
  #   azimuths : azimuths (degrees)
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
  max_ = max(his[, 1])
  cbase = max_ %/% 33 + 1       # number of elements for each point in the plot
  d1 = 21
  length_ = 22.5		# fixed circumference radius and plot width and heigth

  center_x = 0
  center_y = 0
  plot(range(length_, -length_), range(length_, -length_), type = "n", axes = FALSE, ann = FALSE)

  DrawCircle(0, 0, d1, border = "black") 
  points(0, 0, pch = 16, col = "black")

  # vertical and horizontal lines and labels

	text(center_x, center_y + d1, "0"   , adj = c(0, 0), pos = 3, col = "black", cex = 1.25)
	text(center_x - d1, center_y, "270" , adj = c(0, 0), pos = 2, col = "black", cex = 1.25)
	text(center_x, center_y - d1, "180" , adj = c(0, 0), pos = 1, col = "black", cex = 1.25)
	text(center_x + d1, center_y, "90"  , adj = c(0, 0), pos = 4, col = "black", cex = 1.25)
 
	segments(center_x, center_y - d1, center_x, center_y + d1, col = "black", lwd = 2)
	segments(center_x + d1, center_y, center_x - d1, center_y, col = "black", lwd = 2)

  # diagonal lines and labels

	x1 = cos(ToRadians(45)) * d1 
	y1 = sin(ToRadians(45)) * d1 
	
	segments(-x1, -y1, x1, y1, col = "black", lwd = 2) 
	segments(x1, -y1, -x1, y1, col = "black", lwd = 2)

    text(x1, y1, "45", adj = c(-0.25,-0.25), col = "black", cex = 1.25)
    text(x1, -y1, "135", adj = c(-0.25,1), col = "black", cex = 1.25)
    text(-x1, y1, "315", adj = c(1.25,-0.25), col = "black", cex = 1.25)
    text(-x1, -y1,"225", adj = c(1.25,1.25), col = "black", cex = 1.25)

  # scale circles and labels
	
  DrawCircle(0, 0, d1/3, border = "black")      
  DrawCircle(0, 0, 2*d1/3, border = "black")   
  lb33 = 33 * (cbase + 1) / 3
  text(center_x + d1/3, center_y, 2 * lb33, adj = c(-0.2, -0.5), col = "black", cex = 1)
  text(center_x + 2*d1/3, center_y, lb33, adj = c(-0.2, -0.5),  col = "black", cex = 1)
  text(center_x, center_y + d1/3, 2 * lb33, adj = c(-0.2, -0.5), col = "black", cex = 1)
  text(center_x, center_y + 2*d1/3, lb33, adj = c(-0.2, -0.5),  col = "black", cex = 1)

  # accumulated points
   
  for (i in 0:359) {
    h = his[i + 1, 1] / cbase  # elements/point as a function of absolute frequency of first classes
    if (h > 0) {
      for (g in 1:h - 1) { 
        radian = ToRadians(90 - i)  
        x = cos(radian) * (d1 - ((d1 * 0.025) * g))
        y = sin(radian) * (d1 - ((d1 * 0.025) * g))
        points(x, y, cex = 0.8, pch = 16, col = "skyblue3")
      }
    }
  }

  azimuth = MeanAzimuth(azimuths)
  radian = ToRadians(90 - azimuth)

  x = cos(radian) * (d1 + (d1 * 0.1))
  y = sin(radian) * (d1 + (d1 * 0.1))

  # number of elements/point and sample size

  n = length(azimuths)
  labp = paste("Each point represents", cbase, "element(s)")
  text(center_x - d1 - 3.5, center_y + d1 + 1.3, labels = labp, pos = 4, col = "black", cex = 1)
  labp = paste("Sample size, n =", n)
  text(center_x - d1 - 3.5, center_y + d1 + 0.1, labels = labp, pos = 4, col = "black", cex = 1)

  # mean azimuth arrow and condidence interval arc
  # the arrow is drawn only for VonMisesParameter < 0.90

  vm = VonMisesParameter(azimuths)

  if (vm >= 0.9)  {
    arrows(center_x, center_y, x, y, length=0.2, code=2, col="red", lty=par("lty"), lwd=2)
    } 
  else   {
    text(center_x - d1 - 3.5, center_y + d1 - 42, "Concentration is low,", pos=4, col="black", cex=1)
    text(center_x - d1 - 3.5, center_y + d1 - 43, "the mean azimuth is not drawn", pos=4, col="black", cex=1)
    }

  # confidence interval

  if (vm >= 0.9)  {
    module = MeanModule(azimuths)
    ci = ConfidenceInterval(n, azimuth, module, vm)

    xmin = cos(ToRadians(90 - ci[1])) * (d1 + (d1 * 0.1))
    xmax = cos(ToRadians(90 - ci[2])) * (d1 + (d1 * 0.1))
    ymin = sin(ToRadians(90 - ci[1])) * (d1 + (d1 * 0.1))
    ymax = sin(ToRadians(90 - ci[2])) * (d1 + (d1 * 0.1))
    DrawArc(center_x, center_y, radius = d1 + (d1 * 0.1), 
        angle1 = ToRadians(90 - ci[1]), angle2 = ToRadians(90 - ci[2]), n = 35, col = "red", lwd = 2)
      
    points(xmin, ymin, col = "red", pch = 19)
    points(xmax, ymax, col = "red", pch = 19)
    }
	
  if (SVGf == 1) {
   dev.off()
   print(paste("Plot has been saved as SVG graphic file '",fileSVGname,"' in",getwd()))
  }
}

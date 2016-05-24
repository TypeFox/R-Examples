DrawHistogram <- function (azimuths, ClassSize = 15, SVGf = 0) 
{

  # Plots a circular histogram of azimuth values
  # The mean vector and confidence interval are plotted if the Von Mises parameter < 0.9
  #
  # Args:
  #   azimuths : a vector with azimuths of the vector data (degrees)
  #   ClassSize : class width ()
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
  par(fg="white", mar=c(0,0,0,0)) 

  his  = Histogram(azimuths, ClassSize)
  max_ = max(his[, 2]) * 105        

  d1 = round(max_)
  d2 = round(d1 - (max_/4))
  d3 = round(d2 - (max_/4))
  d4 = round(d3 - (max_/4))

  length_ = d1 * 1.2
  center_x = 0
  center_y = 0
  plot(range(length_, -length_), range(length_, -length_), type = "n", axes = FALSE, ann = FALSE)

  DrawCircle(0, 0, d1, border = "black")
  DrawCircle(0, 0, d2, border = "black")
  DrawCircle(0, 0, d3, border = "black")
  DrawCircle(0, 0, d4, border = "black")

  # number of elements

  n = length(azimuths)
  labp = paste("Sample size, n =",n)
  text(center_x - (d1 * 1.15), center_y + d1, labels = labp, pos = 4, col = "black", cex = 1)
  
  # diagonal lines and labels

  x1 = cos(ToRadians(45)) * d1 
  y1 = sin(ToRadians(45)) * d1 

  segments(-x1, -y1, x1, y1, col = "black", lwd = 2) 
  segments(x1, -y1, -x1, y1, col = "black", lwd = 2)

  text(x1, y1, "45",  adj = c(-0.25,-0.25), col = "black", cex = 1.25)
  text(x1, -y1, "135", adj = c(-0.25,1), col = "black", cex = 1.25)
  text(-x1, y1, "315", adj = c(1.25,-0.25), col = "black", cex = 1.25)
  text(-x1, -y1, "225", adj = c(1.25,1.25), col = "black", cex = 1.25)

  # vertical and horizontal lines, labels and sectors

  segments(center_x, center_y - d1, center_x, center_y + d1, col="black", lwd = 2)
  segments(center_x + d1, center_y, center_x - d1, center_y, col="black", lwd = 2)

  # primary labels

  text(center_x, center_y + d1, "0", adj = c(0.5,-0.75), cex = 1.25, col = "black")
  text(center_x - d1, center_y, "270", adj = c(1.25, 0.25), cex = 1.25, col = "black")
  text(center_x, center_y - d1, "180", adj = c(0.5,1.5), cex = 1.25, col = "black")
  text(center_x + d1, center_y, "90", adj = c(-0.5,0.25), cex = 1.25, col = "black")

  # scale labels

  text(center_x, center_y + d1, paste(d1,"%"), adj=c(-0.1,1.2), cex = 1, col = "black")
  text(center_x, center_y + d2, paste(d2,"%"), adj=c(-0.1,1.2), cex = 1, col = "black")
  text(center_x, center_y + d3, paste(d3,"%"), adj=c(-0.1,1.2), cex = 1, col = "black")
  text(center_x, center_y + d4, paste(d4,"%"), adj=c(-0.1,1.2), cex = 1, col = "black")

  text(center_x, center_y - d1, paste(d1,"%"), adj=c(1.1,-0.5), cex = 1, col = "black")
  text(center_x, center_y - d2, paste(d2,"%"), adj=c(1.1,-0.5), cex = 1, col = "black")
  text(center_x, center_y - d3, paste(d3,"%"), adj=c(1.1,-0.5), cex = 1, col = "black")
  text(center_x, center_y - d4, paste(d4,"%"), adj=c(1.1,-0.5), cex = 1, col = "black")

  # label d1 is not plotted (overlap)
  text(center_x + d2, center_y, paste(d2,"%"), adj=c(-0.1,-0.3), cex = 1, col = "black")
  text(center_x + d3, center_y, paste(d3,"%"), adj=c(-0.1,-0.3), cex = 1, col = "black")
  text(center_x + d4, center_y, paste(d4,"%"), adj=c(-0.1,-0.3), cex = 1, col = "black")

  text(center_x - d1, center_y, paste(d1,"%"), adj=c(-0.1,-0.3), cex = 1, col = "black")
  text(center_x - d2, center_y, paste(d2,"%"), adj=c(-0.1,-0.3), cex = 1, col = "black")
  text(center_x - d3, center_y, paste(d3,"%"), adj=c(-0.1,-0.3), cex = 1, col = "black")
  text(center_x - d4, center_y, paste(d4,"%"), adj=c(-0.1,-0.3), cex = 1, col = "black")

   # angular sectors   

  angle = 360 / (length(his) / 2)
  for (i in 1:(length(his)/2)) {
    c_ini = (i * angle) - angle
    c_fin = i * angle
    c_len = his[i, 2] * 100
    coor_x = c(0, c_len * sin(ToRadians(c_fin)), c_len * sin(ToRadians(c_ini)), 0)
    coor_y = c(0, c_len * cos(ToRadians(c_fin)), c_len * cos(ToRadians(c_ini)), 0)
    polygon(coor_x, coor_y, border = "Black", col = "SkyBlue3")
  }

  # mean azimuth

  vm = VonMisesParameter(azimuths)
  if (vm >= 0.53)  {	
    azimuth = MeanAzimuth(azimuths)
    radian = ToRadians(90 - azimuth)  
    x = cos(radian) * (d1 * 1.1)
    y = sin(radian) * (d1 * 1.1)
    arrows(center_x, center_y, x, y, length = 0.2, code = 2, col = "red", lty = par("lty"), lwd = 2)
  }
  else {
    text(center_x  - (d1 * 1.15), center_y - d1,"Concentration is low,", 
      pos = 4, col = "black", cex = 1)
    text(center_x -  (d1 * 1.15), center_y - (d1 * 1.06),"the mean azimuth is not drawn", 
      pos = 4, col = "black", cex = 1)
  }

 if (SVGf == 1) {
   dev.off()
   print(paste("Plot has been saved as SVG graphic file '",fileSVGname,"' in",getwd()))
 }

}

DrawModuleAndAzimuthDistribution <- function (data_x, data_y, SVGf = 0) 
{

  # Plots a graphic of all the vectors from a common origin (0,0)  
  #
  # Args:
  #   data_x : vector components (increments) over the X axis
  #   data_y : vector components (increments) over the Y axis
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

  #rectangular_vectors <- array(c(data_x, data_y), dim=c(length(data_x), 2))
  num_data = length(data_x)
  x = data_x
  y = data_y
  module = sqrt(x * x + y * y)
  max_ = max(module) + 1

  scale_factor = 24 / max_

  d1 = round(max_ * scale_factor)
  d2 = round(d1 * 0.75)
  d3 = round(d1 * 0.50)
  d4 = round(d1 * 0.25)

  length_ = d1 + 10
  center_x = 0
  center_y = 0
  plot(range(length_, -length_), range(length_, -length_), type = "n",
    axes = FALSE, ann = FALSE)

  DrawCircle(0, 0, d1, border = "black", lty = 1, lwd = 1)
  DrawCircle(0, 0, d2, border = "black", lty = 1, lwd = 1)
  DrawCircle(0, 0, d3, border = "black", lty = 1, lwd = 1)
  DrawCircle(0, 0, d4, border = "black", lty = 1, lwd = 1)

  # diagonal lines and labels

  x1 = cos(ToRadians(45)) * d1 
  y1 = sin(ToRadians(45)) * d1 

  segments(-x1, -y1, x1, y1, col = "black", lwd = 2) 
  segments(x1, -y1, -x1, y1, col = "black", lwd = 2)

  text(x1, y1, "45",  adj = c(-0.25, -0.25), col = "black", cex = 1.25)
  text(x1, -y1, "135", adj = c(-0.25, 1), col = "black", cex = 1.25)
  text(-x1, y1, "315", adj = c(1.25, -0.25), col = "black", cex = 1.25)
  text(-x1, -y1, "225", adj = c(1.25, 1.25), col = "black", cex = 1.25)

  # vertical and horizontal lines and labels

  segments(center_x, center_y - d1, center_x, center_y + d1, col = "black", lwd = 2)
  segments(center_x + d1, center_y, center_x - d1, center_y, col = "black", lwd = 2)

  # primary labels

  text(center_x, center_y + d1, "0", adj = c(0.5,-1.0), cex = 1.25, col = "black")
  text(center_x - d1, center_y, "270", adj = c(1.25, 0.25), cex = 1.25, col = "black")
  text(center_x, center_y - d1, "180", adj = c(0.5, 2.0), cex = 1.25, col = "black")
  text(center_x + d1, center_y, "90", adj = c(-0.5, 0.25), cex = 1.25, col = "black")

   # scale labels

  r1 = round(d1 / scale_factor)
  r2 = round(d2 / scale_factor)
  r3 = round(d3 / scale_factor)
  r4 = round(d4 / scale_factor)

    text(center_x, center_y + d1, r1, adj = c(-0.1,1.2), cex = 1, col = "black")
    text(center_x, center_y + d2, r2, adj = c(-0.1,1.2), cex = 1, col = "black")
    text(center_x, center_y + d3, r3, adj = c(-0.1,1.2), cex = 1, col = "black")
    text(center_x, center_y + d4, r4, adj = c(-0.1,1.2), cex = 1, col = "black")

    text(center_x, center_y - d1, r1, adj = c(1.1,-0.5), cex = 1, col = "black")
    text(center_x, center_y - d2, r2, adj = c(1.1,-0.5), cex = 1, col = "black")
    text(center_x, center_y - d3, r3, adj = c(1.1,-0.5), cex = 1, col = "black")
    text(center_x, center_y - d4, r4, adj = c(1.1,-0.5), cex = 1, col = "black")

    #text(center_x + d1, center_y, d1)  # do not plot
    text(center_x + d2, center_y, r2, adj = c(-0.1,-0.3), cex = 1, col = "black")
    text(center_x + d3, center_y, r3, adj = c(-0.1,-0.3), cex = 1, col = "black")
    text(center_x + d4, center_y, r4, adj = c(-0.1,-0.3), cex = 1, col = "black")

    text(center_x - d1, center_y, r1, adj = c(-0.1,-0.3), cex = 1, col = "black")
    text(center_x - d2, center_y, r2, adj = c(-0.1,-0.3), cex = 1, col = "black")
    text(center_x - d3, center_y, r3, adj = c(-0.1,-0.3), cex = 1, col = "black")
    text(center_x - d4, center_y, r4, adj = c(-0.1,-0.3), cex = 1, col = "black")

  # vectors and mean vector

  for (i in 1:(num_data))  { 
    arrows(center_x, center_y, data_x[i] * scale_factor, 
       data_y[i] * scale_factor, length = 0.0, code = 2, 
       col = "skyblue3", lwd = 1)

    }

  xmean = sum(data_x) / num_data
  ymean = sum(data_y) / num_data
  arrows(0, 0, xmean * scale_factor, ymean * scale_factor, length = 0.15, code = 0, col = "red", 
       lty = 1, lwd = 2)

  # label and number of elements

  labp = paste("Red line is the mean vector")
  text(center_x - d1 - 3.5, center_y + d1 + 3.3, labels = labp, pos = 4, col = "black", cex = 1)
  labp = paste("Sample size, n =",num_data)
  text(center_x - d1 - 3.5, center_y + d1 + 1.3, labels = labp, pos = 4, col = "black", cex = 1)

  if (SVGf == 1) {
   dev.off()
   print(paste("Plot has been saved as SVG graphic file '",fileSVGname,"' in",getwd()))
  }

}

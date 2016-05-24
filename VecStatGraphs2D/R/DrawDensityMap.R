DrawDensityMap <- function (data_x, data_y, PercentageOutliers = 5, PaintPoint = FALSE, 
    Div = 250, HarmonicMean = FALSE, PaintAxis = FALSE, SVGf = 0) 
{

  # Plots a graphic of density of the vector endnodes placing the initial
  # initial nodes at the origin (0,0)  
  #
  # Args:
  #   data_x : vector components (increments) over the X axis
  #   data_y : vector components (increments) over the Y axis
  #   PercentageOutliers: the defined %  most distant nodes are identified 
  #   PaintPoint: TRUE > the end nodes are plotted as points
  #   Div : defines the graphic resolution of the color map
  #   HarmonicMean : TRUE > uses the harmonic mean to interpolate values
  #                  FALSE > uses the arithmetic mean 
  #   PaintAxs: plots complementary axis and labels
  # SVG: 
  #     0: the plot is showed only in the graphic window 
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

    #require(MASS)
    data_ = array(c(data_x, data_y), dim = c(length(data_x), 2))

    plot.new()
    par(fg="white", mar=c(4,4,4,4), xaxt="s", yaxt="s") 

   # axis adjust in plot area

    minx = min(data_[,1])
    miny = min(data_[,2])
    maxx = max(data_[,1])
    maxy = max(data_[,2])

    xminp = min(0, min(minx, -maxx))
    xmaxp = max(0, max(maxx, -minx))
    yminp = min(0, min(miny, -maxy))
    ymaxp = max(0, max(maxy, -miny))

    axisl = max(xmaxp, ymaxp) * 1.05

    f1 <- kde2d(data_[,1], data_[,2], n=Div, lims=c(-axisl, axisl, -axisl, axisl))

   # end of axis adjust

    ramp <- colorRamp(c("white", "yellow", "orange", "red"))
    color = rgb(ramp(seq(0, 1, length = 13)), maxColorValue = 255)
	
    image(f1, col=color, xlab="", ylab="")

    grid(lty=1,col="gray")
    abline(v=0, lwd=2, col="black")
    abline(h=0, lwd=2, col="black")
    abline(0,1,lwd=1, col="gray")
    abline(0,-1,lwd=1, col="gray")

   # title and sample size

    cpx = ArithmeticMean(data_x)
    cpy = ArithmeticMean(data_y)
    
    pu <- par()$usr 
    n = NumberOfElements(data_x)
  
    if (cpy < 0 | cpx > 0 ) { 
    text(pu[1] - pu[1]/10, pu[4] - pu[4]/8, "Density map", cex=1.2, pos=4, col="black") 
    text(pu[1] - pu[1]/10, pu[4] - pu[4]/9 - pu[4]/10, paste("Sample size, n =",n), pos=4, col="black")
	}
    else {
    text(pu[1] - pu[1]/10, pu[3] + pu[4]/9 + pu[4]/10, "Density map", cex=1.2, pos=4, col="black")
    text(pu[1] - pu[1]/10, pu[3] + pu[4]/8, paste("Sample size, n =",n), pos=4, col="black")  
	}

   # axis

    if (PaintAxis == TRUE) {
        axis(side=1, pos=0, labels=TRUE,las=1, col="black", ticks=FALSE, font=3)
        axis(side=2, pos=0, las=1, col="black", ticks=FALSE, font=3)
    }

    if (PaintPoint == TRUE) {
        if (HarmonicMean == FALSE) {
            Vec1 = VectorsToPolar(data_)
            Vec = Vec1[, 1]
        }
        else {
            Vec = AllHarmonicMean(data_)
        }

        cant = as.integer(length(data_[, 1]) * (PercentageOutliers/100))
        data_ = data_[order(Vec, decreasing = TRUE), ]
        for (i in 1:length(data_[, 1])) {
            if (i > cant) {
                points(data_[i, 1], data_[i, 2], cex = 0.5, col = "SkyBlue3", 
                  pch = 19)
            }
            else {
                points(data_[i, 1], data_[i, 2], cex = 0.5, col = "Red", 
                  pch = 19)
            }
        }
    }
	
  if (SVGf == 1) {
   dev.off()
   print(paste("Plot has been saved as SVG graphic file '",fileSVGname,"' in",getwd()))
  }
}

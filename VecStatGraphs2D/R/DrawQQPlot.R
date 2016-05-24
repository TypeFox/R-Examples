DrawQQPlot <- function (azimuths, SVGf = 0) 
{
  # Plots a Q-Q graphic of azimuths  
  #
  # Args:
  #   azimuths : azimuths (degrees)
  #   SVG: 
  #     0: the plot is showed only un the graphic window 
  #     1: the plot is saved as SVG graphic 
  #   
  # Returns:
  #   Plot the graphic or save the graphic as SVG
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

    percent <- (length(azimuths) * 20) / 100
    sort_vectors = ToRadians(sort(azimuths))
    v = (sort_vectors/(2 * pi))
    x <- array(NA, c(2, length(azimuths)))
    x[2, ] <- v
    x[1, ] <- (seq(1, length(azimuths)) / (length(azimuths) + 1))
    y <- array(NA, c(2, 2 * percent))
    y[2, 1:percent] <- v[1:percent] + 1
    y[2, (percent + 1):(2 * percent)] <- v[(length(azimuths) - 
        (percent - 1)):(length(azimuths))] - 1
    y[1, 1:percent] <- x[1, 1:percent] + 1
    y[1, (percent + 1):(2 * percent)] <- x[1, (length(azimuths) - 
        (percent - 1)):(length(azimuths))] - 1
    z <- array(NA, c(2, 2 * percent + length(azimuths)))
    z[, ] <- c(x[, ], y[, ])
    
    plot.new()
    par(fg="white", mar=c(5,5,5,5)) 

    plot(z[1,], z[2,], type="p", col="skyblue3", cex.lab=1.25, 
	ylab="Data Quantiles", pch=16, xlab="Uniform Theoretical Quantiles")

    box(lty=1, lwd=2, col="black")
    grid(lty=1)
    abline(a=0, b=1, lwd=2, col="black")

  # Title and sample size

	pu <- par()$usr 
	text(pu[1] + 0.1, pu[4] - 0.2, "Q-Q plot of azimuths",cex=1.2, pos=4, col="black") 
	n = NumberOfElements(azimuths)
   	text(pu[1] + 0.1, pu[4] - 0.3, paste("Sample size, n =",n), pos=4, col="black")

		
  if (SVGf == 1) {
   dev.off()
   print(paste("Plot has been saved as SVG graphic file '",fileSVGname,"' in",getwd()))
  }

}

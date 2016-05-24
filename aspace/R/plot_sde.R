"plot_sde" <-
function(plotnew=TRUE, plotSDEaxes=FALSE, plotweightedpts=FALSE, weightedpts.col='black', weightedpts.pch=19, plotpoints=TRUE, points.col='black', 
         points.pch=1, plotcentre=TRUE, centre.col='black', centre.pch=19, titletxt="Title", xaxis="Easting (m)", yaxis="Northing (m)", 
		 sde.col='black', sde.lwd=2, jpeg=FALSE, ...) {
		 
  #=======================================================
  #
  #  TITLE:     STANDARD DEVIATION ELLIPSE (SDE) PLOT FUNCTION
  #  FUNCTION:  plot_sde()
  #  AUTHOR:    RANDY BUI, RON BULIUNG, TARMO K. REMMEL
  #  DATE:      March 28, 2011
  #  CALLS:     jpeg()
  #  NOTES:     THE r.SDE OBJECT IS REQUIRED (GENERATED FROM THE CALC_SDE FUNCTION)
  #             TO PLOT THE SDE. THE PAR(...) OPTION ALLOWS 
  #             FOR ADDITIONAL GRAPHICAL PARAMETERS.
  #
  #=======================================================

	# ALLOW THE USER TO ACCESS AND MODIFY THE LIST OF GRAPHICAL PARAMETERS OF THE CURRENT PLOT.
	par(...)

	  # DEFINE PLOT EXTENTS
	  min.x <- min( (r.SDE$CENTRE.x - max(r.SDE$Sigma.x,r.SDE$Sigma.y)), min(r.SDE$points[,1]) )
	  max.x <- max( (r.SDE$CENTRE.x + max(r.SDE$Sigma.x,r.SDE$Sigma.y)), max(r.SDE$points[,1]) )
	  min.y <- min( (r.SDE$CENTRE.y - max(r.SDE$Sigma.x,r.SDE$Sigma.y)), min(r.SDE$points[,2]) )
	  max.y <- max( (r.SDE$CENTRE.y + max(r.SDE$Sigma.x,r.SDE$Sigma.y)), max(r.SDE$points[,2]) )

	  # JPEG OUTPUT OPTION
	  if(jpeg) {
		jpeg(filename = paste("SDE",r.SDE$id,".jpg", sep=""), width = 600, height = 600, pointsize = 12, quality = 90, bg = "white", res = NA)
	  }
	  
	  # PLOT THE STANDARD DEVIATION ELLIPSE
	  if(plotnew) {
	  plot(1, type="n", asp=1, xlab=xaxis, ylab=yaxis, xlim=c(min.x, max.x), ylim=c(min.y, max.y))
	  }
	  
	  lines(r.SDE$coordsSDE, col=sde.col, lwd=sde.lwd)
	  title(paste(titletxt, sep=""))
	  
	  if(plotSDEaxes) {    
	    # PLOT HALF-AXES (SIGMA-Y) FOR ELLIPSE IN GREEN
	    xprime <- r.SDE$CENTRE.x + ( r.SDE$Sigma.y * cos_d(90-r.SDE$theta) )
	    yprime <- r.SDE$CENTRE.y + ( r.SDE$Sigma.y * sin_d(90-r.SDE$theta) )
	    segments(r.SDE$CENTRE.x, r.SDE$CENTRE.y, xprime, yprime, col=3)
	    
	    # PLOT HALF-AXES (SIGMA-X) FOR ELLIPSE IN BLUE
	    xprime <- r.SDE$CENTRE.x + ( r.SDE$Sigma.x * cos_d(r.SDE$theta) )
	    yprime <- r.SDE$CENTRE.y - ( r.SDE$Sigma.x * sin_d(r.SDE$theta) )
	    segments(r.SDE$CENTRE.x, r.SDE$CENTRE.y, xprime, yprime, col=4)
	  }
	  
	  if(plotweightedpts) {
		# PLOT THE WEIGHTED POINTS	
		points(r.SDE$points, cex=r.SDE$weights, col=weightedpts.col, pch=weightedpts.pch)
	  }	  
	  
	  if(plotpoints) {
	    # PLOT THE POINTS	
        points(r.SDE$points, col=points.col, pch=points.pch)
      }
		  
	  if(plotcentre) {
		# ADD THE CENTRE POINT (NON-WEIGHTED/WEIGHTED/USER-DEFINED)
		points(r.SDE$CENTRE.x, r.SDE$CENTRE.y, col=centre.col, pch=centre.pch)
	  }

	  if(jpeg) {
		# TURN OFF JPEG DEVICE WHEN DONE
		dev.off()
	  }
  }
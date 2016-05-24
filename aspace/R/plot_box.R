"plot_box" <-
function(plotnew=TRUE, plothv=FALSE, plotweightedpts=FALSE, weightedpts.col='black', weightedpts.pch=19, plotpoints=TRUE, 
         points.col='black', points.pch=1, plotcentre=TRUE, centre.col='black', centre.pch=19, titletxt="Title", 
		 xaxis="Easting (m)", yaxis="Northing (m)", box.col='black', box.lwd=2, jpeg=FALSE, ...) {

  #=======================================================
  #
  #  TITLE:     STANDARD DEVIATION BOX PLOT FUNCTION 
  #  FUNCTION:  plot_box()
  #  AUTHOR:    RANDY BUI, RON BULIUNG, TARMO REMMEL
  #  DATE:      March 28, 2011
  #  NOTES:     THE r.BOX OBJECT IS REQUIRED (GENERATED FROM 
  #             CALC_BOX FUNCTION) TO PROPERLY PLOT THE SD BOX. 
  #             THE PAR(...) OPTION ALLOWS FOR ADDITIONAL 
  #             GRAPHICAL PARAMETERS.
  #
  #=======================================================

	# ALLOW THE USER TO ACCESS AND MODIFY THE LIST OF GRAPHICAL PARAMETERS OF THE CURRENT PLOT.
	par(...)
		
	      # DEFINE PLOT EXTENTS
	      min.x <- min( (r.BOX$CENTRE.x - r.BOX$SDD), min(r.BOX$points[,1]) )
	      max.x <- max( (r.BOX$CENTRE.x + r.BOX$SDD), max(r.BOX$points[,1]) )
	      min.y <- min( (r.BOX$CENTRE.y - r.BOX$SDD), min(r.BOX$points[,2]) )
	      max.y <- max( (r.BOX$CENTRE.y + r.BOX$SDD), max(r.BOX$points[,2]) )

	  # JPEG OUTPUT OPTION
	  if(jpeg) {
		jpeg(filename = paste("BOX",r.BOX$id,".jpg", sep=""), width = 600, height = 600, pointsize = 12, quality = 90, bg = "white", res = NA)
	  }
	  
	  # PLOT THE STANDARD DEVIATION BOX
	  if(plotnew) {
	    plot(1, type="n", asp=1, xlab=xaxis, ylab=yaxis, xlim=c(min.x, max.x), ylim=c(min.y, max.y))
		}
	  
	  title(paste(titletxt, sep=""))
		  
      # PLOT STANDARD DEVIATION BOX CENTRED ON (EITHER USER-DEFINED CENTRE, MEAN CENTRE, OR WEIGHTED MEAN CENTRE - DEPENDING ON WHICH IS SELECTED)
      segments(r.BOX$NW.coord[1], r.BOX$NW.coord[2], r.BOX$NE.coord[1], r.BOX$NE.coord[2], col=box.col, lwd=box.lwd)
      segments(r.BOX$SW.coord[1], r.BOX$SW.coord[2], r.BOX$SE.coord[1], r.BOX$SE.coord[2], col=box.col, lwd=box.lwd)
      segments(r.BOX$NE.coord[1], r.BOX$NE.coord[2], r.BOX$SE.coord[1], r.BOX$SE.coord[2], col=box.col, lwd=box.lwd)
      segments(r.BOX$NW.coord[1], r.BOX$NW.coord[2], r.BOX$SW.coord[1], r.BOX$SW.coord[2], col=box.col, lwd=box.lwd)

	  if(plothv) {
		# DRAW HORIZONTAL AND VERTICAL LINES THROUGH CENTRE (NON-WEIGHTED/WEIGHTED/USER-DEFINED)
		abline(h=r.BOX$CENTRE.y, col=1, lty = 2)
		abline(v=r.BOX$CENTRE.x, col=1, lty = 2)
	  }

	  if(plotweightedpts) {
		# PLOT THE WEIGHTED POINTS	
		points(r.BOX$points, cex=r.BOX$weights, col=weightedpts.col, pch=weightedpts.pch)
	  }	  
	  
	  if(plotpoints) {
		# PLOT THE POINTS	
		points(r.BOX$points, col=points.col, pch=points.pch)
	  }
	  
	  if(plotcentre) {
		# ADD THE CENTRE POINT (NON-WEIGHTED/WEIGHTED/USER-DEFINED)
		points(r.BOX$CENTRE.x, r.BOX$CENTRE.y, col=centre.col, pch=centre.pch)
	  }
	  
	  if(jpeg) {
		# TURN OFF JPEG DEVICE WHEN DONE
		dev.off()
	  }	  
  }


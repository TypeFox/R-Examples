"plot_centres" <-
function(plotnew=FALSE, plotSDE=FALSE, xaxis="Easting (m)", yaxis="Northing (m)", robject=NULL, plotweightedpts=FALSE, weightedpts.col='black', weightedpts.pch=19, plotpoints=TRUE, points.col='black', 
         points.pch=1, plotcentre=FALSE, centre.col='black', centre.pch=19, plotcentral=FALSE, central.col='green', central.pch=19, 
		 plotCF2PTS=FALSE, CF2PTS.col='orange', CF2PTS.pch=19, plotmedian=FALSE, median.col='blue', median.pch=17, plotCMD=FALSE, 
		 CMD.col='red', CMD.pch=17, ...) {

  #=======================================================
  #
  #  TITLE:     CENTRES PLOT FUNCTION
  #  FUNCTION:  plot_sdd()
  #  AUTHOR:    RANDY BUI, RON BULIUNG, TARMO K. REMMEL
  #  DATE:      December 08, 2010
  #  NOTES:     THE r.SDD/r.SDE/r.BOX OBJECT IS REQUIRED (GENERATED FROM THE 
  #             CALC_SDD/CALC_SDE/CALC_BOX FUNCTION) TO PLOT. 
  #             THE PAR(...) OPTION ALLOWS FOR ADDITIONAL GRAPHICAL 
  #             PARAMETERS.
  #
  #=======================================================

	# ALLOW THE USER TO ACCESS AND MODIFY THE LIST OF GRAPHICAL PARAMETERS OF THE CURRENT PLOT.
	par(...)
	
		  # PLOT THE STANDARD DEVIATION DISTANCE CIRCLE
		  if(plotnew) {
		  
		  # DEFINE PLOT EXTENTS
	      min.x <- min( (robject$CENTRE.x - robject$SDD), min(robject$points[,1]) )
	      max.x <- max( (robject$CENTRE.x + robject$SDD), max(robject$points[,1]) )
	      min.y <- min( (robject$CENTRE.y - robject$SDD), min(robject$points[,2]) )
	      max.y <- max( (robject$CENTRE.y + robject$SDD), max(robject$points[,2]) )
		  
		  if(plotSDE){
		  min.x <- min( (r.SDE$CENTRE.x - max(r.SDE$Sigma.x,r.SDE$Sigma.y)), min(r.SDE$points[,1]) )
		  max.x <- max( (r.SDE$CENTRE.x + max(r.SDE$Sigma.x,r.SDE$Sigma.y)), max(r.SDE$points[,1]) )
		  min.y <- min( (r.SDE$CENTRE.y - max(r.SDE$Sigma.x,r.SDE$Sigma.y)), min(r.SDE$points[,2]) )
		  max.y <- max( (r.SDE$CENTRE.y + max(r.SDE$Sigma.x,r.SDE$Sigma.y)), max(r.SDE$points[,2]) )
		  }

	      plot(1, type="n", asp=1, xlab=xaxis, ylab=yaxis, xlim=c(min.x, max.x), ylim=c(min.y, max.y))
		  }
		  
	      if(plotweightedpts) {
	        # PLOT THE WEIGHTED POINTS	
	        points(robject$points, cex=robject$weights, col=weightedpts.col, pch=weightedpts.pch)
	      }
		  
	      if(plotpoints) {
	        # PLOT THE POINTS	
	        points(robject$points, col=points.col, pch=points.pch)
	      }
		  
	      if(plotcentre) {
	        # ADD THE CENTRE POINT (NON-WEIGHTED/WEIGHTED/USER-DEFINED)
	        points(robject$CENTRE.x, robject$CENTRE.y, col=centre.col, pch=centre.pch)
	      }
		  
		  if(plotcentral) {
	        # IDENTIFY THE CENTRAL FEATURE
	        points(r.CF$CF.x, r.CF$CF.y, col=central.col, pch=central.pch)
	      }

		  if(plotmedian) {
	        # ADD MEDIAN CENTRE POINT
	        points(r.median$median.x, r.median$median.y, col=median.col, pch=median.pch)
	      }
		  
		  if(plotCMD) {
	        # ADD THE CENTRE OF MINIMUM DISTANCE POINT
	        points(r.CMD$CMD.x, r.CMD$CMD.y, col=CMD.col, pch=CMD.pch)
	      }
		  
		  if(plotCF2PTS) {
	        # IDENTIFY THE CENTRAL FEATURE
	        points(r.CF2PTS$CF2PTS.x, r.CF2PTS$CF2PTS.y, col=CF2PTS.col, pch=CF2PTS.pch)
	      }
	}
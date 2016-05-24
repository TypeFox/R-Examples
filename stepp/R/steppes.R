#################################################################
#
# steppes.R
#
#############################
# stepp estimate/statistics #
#############################
setClass("steppes",
	   representation(subpop   = "stsubpop",	# stepp subpopulation
				model	   = "stmodel",	# stepp model
				effect   = "ANY",		# list of absolute effect est
				result   = "ANY",		# test statistics
				nperm	   = "numeric"	# number of permutations 0-n
				),
	   prototype=c(NULL, NULL, NULL, NULL, NULL, NULL)
	   )


setMethod("estimate",
	    signature="steppes",
	    definition=function(.Object, sp, model){
		.Object@subpop <- sp
		.Object@model  <- model
		.Object@effect <- estimate(model, sp)
		return(.Object)
	    }
)

setMethod("test",
	    signature="steppes",
	    definition=function(.Object, nperm=100, showstatus=TRUE){
		if (is.null(.Object@subpop) | is.null(.Object@model)){
		  print("You have to estimate the effect first before you can test for interaction !")
		}
		else {
		  .Object@nperm  <- nperm
		  .Object@result <- test(.Object@model, nperm, .Object@subpop, .Object@effect, showstatus=TRUE)
		}
		return(.Object)
	    }
)

setMethod("summary",
	    signature="steppes",
	    definition=function(object){
		summary(object@subpop)
	    }
)

setMethod("print",
	    signature="steppes",
	    definition=function(x, estimate=TRUE, cov=TRUE, test=TRUE, ...){
		n0 <- sum(x@model@coltrt==x@model@trts[1])
		n1 <- sum(x@model@coltrt==x@model@trts[2])
		n  <- n0 + n1
		cat("\n")
		write(paste("Sample size in the 1st treatment : ", n0), file="")
		write(paste("Sample size in the 2nd treatment : ", n1), file="")
		write(paste("Total sample size (excluding missing data) : ", n), file="")
		print(x@model, x, estimate, cov, test, ...)
	    }
)

#
# Internal worker routine for plot
.Stepp.plot <- function(x, y, legendy, pline, at, color, ylabel, xlabel, ncex, tlegend,
				nlas, alpha, pointwise, diff, ci, pv, showss, ylimit, dev,
				together, noyscale, rug, ...){

    			skmObs1 <- x@effect$sObs1
    			skmObs2 <- x@effect$sObs2
    			skmSE1  <- x@effect$sSE1
    			skmSE2  <- x@effect$sSE2
    			nsubpop <- x@subpop@nsubpop
    			npatsub <- x@subpop@npatsub
    			medians <- x@subpop@medianz
    			r1      <- x@subpop@win@r1
    			r2      <- x@subpop@win@r2
    			pvalue  <- x@result$pvalue

    			n <- 3
    			if (class(x@model) == "stmodelGLM" ) noyscale <- TRUE 
    			if (together) n <- 1

			if (dev == "") graphics.off()

    			for(i in 1:n) {
			  if (dev == "") dev.new()
			  else
			  if (dev == "postscript") {
	  		    fname = paste("SteppPlot",as.character(i),".ps",sep="")
	  		    postscript(file=fname)
			  } 
			  else

			  if (dev == "eps") {
	  		    fname = paste("SteppPlot",as.character(i),".eps",sep="")
	  		    postscript(file=fname)
			  } 
			  else 
			  if (dev == "pdf") {
	  		    fname = paste("SteppPlot",as.character(i),".pdf",sep="")
	  		    pdf(file=fname)
			  } 
			  else 
			  if (dev == "png") {
	  		    fname = paste("SteppPlot",as.character(i),".png",sep="")
	  		    png(filename=fname)
			  } 
			  else 
			  if (dev == "bmp") {
	  		    fname = paste("SteppPlot",as.character(i),".bmp",sep="")
	  		    bmp(filename=fname)
			  } 
			  else 
			  if (dev == "tiff") {
	  		    fname = paste("SteppPlot",as.character(i),".tif",sep="")
	  		    tiff(filename=fname)
			  } 
			  else 
			  if (dev == "jpeg") {
	  		    fname = paste("SteppPlot",as.character(i),".jpeg",sep="")
	  		    jpeg(filename=fname)
			  } 
    		      }

    		 	devlst <- dev.list()

			#   generate the first stepp plot
			#     STEPP analysis of treatment effect as measured
			#     by KM/HR or cumulative incidence.
			#
			dev.set(devlst[1])

    			skmObs  <- rep(0, nsubpop * 2)
    			xvalues <- rep(0, nsubpop * 2)
    			skmObs  <- c(skmObs1, skmObs2)
    			xvalues <- c(medians, medians)
    			if (!noyscale) skmObs <- skmObs * 100
    			group   <- rep(1:0, each = nsubpop)
    			lbls    <- rep(" ", nsubpop)
    			ssize   <- rep(" ", nsubpop)
    			for (i in 1:nsubpop)
			  ssize[i] <- paste(c("(n=", npatsub[i], ")"), collapse = "")
    			p 	  <- paste(c("Supremum p-value = ", pvalue), collapse = "")
    			par(mfrow = c(1, 1), omi = c(0.4, 0.4, 0.1, 0.1))

    			if (noyscale) {
			  yl      <- c(min(skmObs), max(skmObs))
			  legendy <- max(skmObs) - legendy
    			}
    			else
    			if (length(ylimit) < 2) yl <- c(0,100)
    			else yl   <- ylimit[1:2]

    			if (noyscale)
			  plot(xvalues, skmObs, axes = TRUE, xaxt="n", ylim = yl, ylab = ylabel, xlab = "")
    			else
      		  plot(xvalues, skmObs, axes = FALSE, ylim = yl, ylab = ylabel, xlab = "")

    			points(xvalues[group == 1], skmObs[group == 1], lty = 2, lwd = 2, pch = 19, type = "o", col = color[1], bg = color[1])
    			points(xvalues[group == 0], skmObs[group == 0], lty = 1, lwd = 2, pch = 24, type = "o", col = color[2], bg = color[2])
    			axis(1, at = xvalues, font = 1)

    			if (!noyscale) axis(2, at = c(0, (0:5) * 20), font = 1)
 
    			if (nlas != 3 & nlas != 2) {
        		  if (showss) {
	    		    mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
          		    mtext(xlabel, side = 1, line = 3.5)
	  		  } else 
	    		  mtext(xlabel, side = 1, line = 2)
    		      }
    			if (nlas == 3 | nlas == 2) {
        		  if (showss) {
	    		    mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
          		    mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
	  		  } else
	  	        mtext(xlabel, side = 1, line = 2)
	  		}
    		 
    		    	legend(min(xvalues), legendy, pch = c(19, 24), lty = c(2,1), lwd = 2, col = color, pt.bg = color, legend = tlegend, bty = "n")
    
    		    	if (pv) {
			  if (is.na(at)) mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
			  else mtext(p, side = 1, at = at, line = pline)
			}
    		    	if (rug) rug(xvalues)
   
		    	#
    		    	if(diff){
			  # pointwise is specified, generate two additional plots
			  #
			  if (!together) dev.set(devlst[2])

       		  if(pointwise) zcrit <- qnorm(1-alpha/2)
       		  else zcrit <- qnorm(1-alpha/(2*nsubpop))
       		  skmObs <- rep(0, nsubpop)
       		  xvalues <- rep(0, nsubpop)
       		  skmObs <- skmObs1-skmObs2
       		  se <- sqrt(skmSE1^2+skmSE2^2)
       		  xvalues <- medians
       		  if (!noyscale) skmObs <- skmObs * 100
       		  lbls <- rep(" ", nsubpop)
       		  ssize <- rep(" ", nsubpop)
       		  for (i in 1:nsubpop) ssize[i] <- paste(c("(n=", npatsub[i], ")"), collapse = "")
       		  p <- paste(c("Supremum p-value = ", pvalue), collapse = "")
       		  par(mfrow = c(1, 1), omi = c(0.4, 0.4, 0.1, 0.1))
 
	 		  if (noyscale){
	   		    ext <- (max(skmObs) - min(skmObs))*0.5
	   		    yl <- c(min(skmObs)-ext, max(skmObs)+ext)
	 		  }
	 		  else
       		  if (length(ylimit) < 4) yl <- c(-100,100)
       		  else yl <- ylimit[3:4]

			  if (noyscale)
         		    plot(xvalues, skmObs, axes = TRUE, xaxt="n", ylim = yl, ylab = paste("Difference in",ylabel), xlab = "")
	 		  else 
         		    plot(xvalues, skmObs, axes = FALSE, ylim = yl, ylab = paste("Difference in",ylabel), xlab = "")

       		  points(xvalues, skmObs, lty = 1, lwd = 2, pch = 19, type = "o", col = color[1], bg = color[1])
	 		  if (rug) rug(xvalues)
       		  if(ci){
	    		    if (!noyscale){
            	      cilim <- skmObs - zcrit*se*100
            	      cilim <- pmax(cilim,-100)
	    		    } else cilim <- skmObs - zcrit*se
           
          		  lines(xvalues, cilim, lty = 2, lwd = 2, col = color[1], bg = color[1])
          		  if (!noyscale){
			    cilim <- skmObs + zcrit*se*100
            	    cilim <- pmin(cilim,100)
	    		  } else cilim <- skmObs + zcrit*se

          		  lines(xvalues, cilim, lty = 2, lwd = 2, col = color[1], bg = color[1])
       	      }
       	      lines(c(min(xvalues),max(xvalues)),c(0,0),lty=1)
       	      axis(1, at = xvalues, font = 1)
       	      if (!noyscale) axis(2, at = c(0, (-5:5) * 20), font = 1)
       	      if (nlas != 3 & nlas != 2) {
           	        if (showss) {
	       	    mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
             	    mtext(xlabel, side = 1, line = 3.5)
	     		  } else
	              mtext(xlabel, side = 1, line = 2)
       	      }
       	      if (nlas == 3 | nlas == 2) {
           		  if (showss) {
	       	    mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
             	    mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
	     		  } else
	       	  mtext(xlabel, side = 1, line = 2)
       	      }
       	      if (pv) {
			  if (is.na(at)) mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
			  else mtext(p, side = 1, at = at, line = pline)
			}

			if (!together) dev.set(devlst[3])
	 	      if (class(x@model) == "stmodelKM" | class(x@model) == "stmodelCI" | class(x@model) == "stmodelCOX"){      
    	   		  Rpvalue<- x@result$HRpvalue
    	   		  logR   <- x@effect$logHR
    	   		  logRSE <- x@effect$logHRSE
			  text1  <- "Supremum HR p-value = "
			  text2  <- "Hazard Ratio"
			} else {
    	   		  Rpvalue<- x@result$logRpvalue
    	   		  logR   <- x@effect$logR
    	   		  logRSE <- x@effect$logRSE
			  if (x@model@glm=="gaussian"){
			    text1  <- "Supremum Effect Ratio p-value = "
			    text2  <- paste(ylabel, "Ratio")
			  } else 
			  if (x@model@glm=="binomial"){
			    text1  <- "Supremum Odds Ratio p-value = "
			    text2  <- "Odds Ratio"
			  } else
			  if (x@model@glm=="poisson"){
			    text1  <- "Supremum Risks Ratio p-value = "
			    text2  <- "Risks Ratio"
			  }
			}
         		  p <- paste(c(text1, Rpvalue), collapse = "")
         		  par(mfrow = c(1, 1), omi = c(0.5, 0.5, 0.5, 0.5))
         		  if (length(ylimit) < 6) yl = c(0,3)
	   		  else yl = ylimit[5:6]
         		  plot(xvalues, exp(logR), axes = FALSE, ylim = yl, ylab = text2, xlab = "")
         		  points(xvalues, exp(logR), lty = 1, lwd = 2, pch = 19, type = "o", col = color[1], bg = color[1])
         		  if(ci){
            	    cilim <- logR - zcrit*logRSE
            	    cilim <- exp(cilim)
            	    lines(xvalues, cilim, lty = 2, lwd = 2, col = color[1], bg = color[1])
            	    cilim <- logR + zcrit*logRSE
            	    cilim <- exp(cilim)
            	    lines(xvalues, cilim, lty = 2, lwd = 2, col = color[1], bg = color[1])
         		  }
         		  lines(c(min(xvalues),max(xvalues)),c(1,1),lty=1)
         		  axis(1, at = xvalues, font = 1)
         		  axis(2, at = c(0, (0:15)*.2), font = 1)
         		  if (nlas != 3 & nlas != 2) {
           		    if (showss) {
	       	      mtext(ssize, side = 1, at = xvalues, line = 2, font = 1, cex = ncex, adj = 0.5, las = nlas)
             	      mtext(xlabel, side = 1, line = 3.5)
	     		    } else
		 	    mtext(xlabel, side = 1, line = 2)
         		  }
         		  if (nlas == 3 | nlas == 2) {
           		    if (showss) {
		 	      mtext(ssize, side = 1, at = xvalues, line = 3.5, font = 1, cex = ncex, adj = 0.5, las = nlas)
             	      mtext(xlabel, side = 1, line = 0.8, outer = TRUE)
	     		    } else
		 	    mtext(xlabel, side = 1, line = 2)
         		  }

         		  if (pv) {
			    if (is.na(at)) mtext(p, side = 1, at = (min(xvalues) + 0.2 * (max(xvalues) - min(xvalues))), line = pline)
			    else mtext(p, side = 1, at = at, line = pline)
			  }
	 	      }

    		    for(i in 1:n) {
		      if (dev != "") dev.off(devlst[i])
    		    }
}# end of worker function


#
# S3 method
plot.steppes <- function(x, y, legendy = 30, pline = -2.5, color = c("red", "black"),
			ylabel = "Specify Timepoint & Endpoint",
			xlabel = "Subpopulations by Median Covariate",
			ncex = 0.7, tlegend = c("Specify 1st Treatment", "Specify 2nd Treatment"),
			nlas = 0, alpha = 0.05, pointwise = FALSE, diff = TRUE, ci = TRUE,
			pv = TRUE, showss = TRUE, ylimit = c(0,100,-100,100,0,3), dev = "",
			together = FALSE, noyscale = FALSE, rug = FALSE, at = NA, ...
		 ){
	return(.Stepp.plot(x, y, legendy, pline, at, color, ylabel, xlabel, ncex, tlegend,
				nlas, alpha, pointwise, diff, ci, pv, showss, ylimit, dev,
				together, noyscale, rug, ...)
		)
}
 
#
# S4 method
setMethod("plot",
	    signature="steppes",
	    definition=function(x, y, 
							# graphics parameters - optional
			legendy = 30,		#   the vertical location of the legend according to the units on the y-axis
			pline = -2.5,		#   the vertical location of the p-value, starting at 0 counting outwards
			color = c("red", "black"),
							#   a vector containing the line colors for the 1st and 2nd trt, respectively
			ylabel = "Specify Timepoint & Endpoint",   	# label for the y-axis
			xlabel = "Subpopulations by Median Covariate", 	# label for the x-axis
    			ncex = 0.7,			#   the size of the text for the sample size annotation
			tlegend = c("Specify 1st Treatment", "Specify 2nd Treatment"),
							#   a vector containing the treatment labels, 1st and 2nd trt, respectively
			nlas = 0,			#   the las paramter (0,1,2,3) - the orientation of the sample size annotation
							# plot options - optional
			alpha = 0.05,		#   sig. level
			pointwise = FALSE,	#   pointwise confidence intervals (pointwise=TRUE),
							#   or confidence bands (pointwise=FALSE, default) to be displayed
			diff = TRUE,		#   generate 2 additional plots comparing the diff between measures
			ci = TRUE,			#   display the conf. interval or band
			pv = TRUE,			#   display the supremum pvalue
			showss = TRUE,		#   display sample size on the x-axis
			ylimit = c(0,100,-100,100,0,3),
							#   y limits for the 3 graphs
			dev = "",			#   graphics device for output; default to Null Device
		      together = FALSE,		#   generate the plots together; default to No
			noyscale = FALSE,		#   do not scale y axis to %
			rug = FALSE,		#   put a rug plot for each
			at = NA,			#   centering position for pline
		      ... ){

	return(plot.steppes(x, y, legendy, pline, color, ylabel, xlabel, ncex, tlegend,
				nlas, alpha, pointwise, diff, ci, pv, showss, ylimit, dev,
				together, noyscale, rug, at, ...)
		)
	}
)

# constructor function for steppes
stepp.test <- function(subpop, model, nperm, showstatus=TRUE){
	result <- new("steppes")
	result <- estimate(result, subpop, model)
	result <- test(result, nperm, showstatus=showstatus)
	return(result)

}



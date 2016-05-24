#'Make a PDF report of several YplantQMC objects
#'
#'@description Produce a report containing standard graphs and summaries of YplantQMC
#'objects. See Details for usage.
#'
#'This function produces a number of standard plots and prints of five
#'different YplantQMC objects. The plots cannot at the moment not be
#'customized; please see below for the component functions that generate the
#'plots (these usually have more options for customization). Or, modify the
#'code of \code{ypreport} as you see fit.
#'
#'To create a report, simply use this command: \preformatted{
#'ypreport(plant=myplant, met=asunnyday, phy=eucleaf, hemi=mycanopy,
#'ypsim=eucsim1) } Where myplant, asunnyday, eucleaf, mycanopy and eucsim1 are
#'objects that you have already generated. The function is flexible : you can
#'generate a report on a subset of the objects, for example only on the
#'hemiphoto and the plant: \preformatted{ ypreport(plant=myplant,
#'hemi=spruceforest) }
#'
#'You may want to check out the following functions, which are used in
#'\code{ypreport}: \describe{ \item{list(list("viewplot"))}{Produces a
#'three-panel plot with side and top views of the plant.}
#'\item{list(list("summary.plant3d"))}{Summarizes a plant}
#'\item{list(list("setHemi"))}{Reads a hemiphoto - also describes the plot
#'function.} \item{list(list("fitdistribution"))}{Fits (and plots) a leaf angle
#'distribution.} \item{list(list("setMet"))}{Constructs (and plots) a weather
#'object.} \item{list(list("plot.leaffile"))}{Plots a leaf.}
#'\item{list(list("YplantDay"))}{A daily Yplant simulation (and a standard
#'plot).} }
#'
#'@param plant An object of class 'plant3d', see \code{\link{constructplant}}
#'@param phy An object of class 'ypphy', see \code{\link{setPhy}}
#'@param hemi An object of class 'yphemi', see \code{\link{setHemi}}
#'@param met An object of class 'ypmet', see \code{\link{setMet}}
#'@param ypsim An object of class 'yplantsim', see \code{\link{YplantDay}}
#'@param filename Optional, the name of the output file
#'@return A PDF is generated in the current working directory.
#'@note Note that \code{ypreport} might fail if you have a PDF open with the
#'same name (as you may have generated this report once before on the same
#'day). Make sure to close the PDF before running \code{ypreport}. If a PDF
#'is generated that can't be opened, use this command: \preformatted{ dev.off()
#'} And try again.
#'@author Remko Duursma
#'@seealso
#'\code{\link{plot.plant3d}},\code{\link{viewplot}},\code{\link{YplantDay}}
#'@keywords misc
#'@export
#'@importFrom gplots textplot
#'@importFrom LeafAngle fitdistribution
ypreport <- function(plant=NULL, phy=NULL, met=NULL, hemi=NULL, ypsim=NULL, filename=NA){

	if(is.na(filename)){
		filename <- paste0("YplantQMC_report_",as.Date(Sys.time()),".pdf")
	}

	# open pdf
	pdf(filename, onefile=TRUE)

	emptyplot <- function(){
		par(mar=c(0,0,0,0),xaxs="i",yaxs="i")
		plot(1, type='n', ann=FALSE,axes=FALSE,xlim=c(0,1),ylim=c(0,1))
	}

	# First page.
	emptyplot()
	text(0.5,0.9, expression(bold("YplantQMC Report")), cex=2)
	text(0.1,0.8, paste("Generated on : ",Sys.time()), cex=1.2, font=3, pos=4)
	abline(h=0.75)

	text(0.1,0.6, "Objects included:", font=2, pos=4, cex=1.3)
	# ypos <- c(0.5,0.4,0.3,0.2,0.1)
	ypos <- seq(0.5,0.15, length=5)
	k <- 1
	if(!is.null(plant)){text(0.1,ypos[k],paste("Plant      :", deparse(substitute(plant))), cex=1.2, pos=4);k <- k + 1}
	if(!is.null(phy)){text(0.1,ypos[k],  paste("Physiology :", deparse(substitute(phy))), cex=1.2, pos=4);k <- k + 1}
	if(!is.null(met)){text(0.1,ypos[k],  paste("Weather    :", deparse(substitute(met))), cex=1.2, pos=4);k <- k + 1}
	if(!is.null(hemi)){text(0.1,ypos[k], paste("Hemiphoto  :", deparse(substitute(hemi))), cex=1.2, pos=4);k <- k + 1}
	if(!is.null(ypsim)){text(0.1,ypos[k],paste("Simulation :", deparse(substitute(ypsim))), cex=1.2, pos=4)}

	if(!is.null(plant) && !is.null(ypsim)){
		if(!identical(plant, ypsim$plant))warning("**-Plant object and plant used for simulation are not the same!")
	}
	if(!is.null(plant) && !is.null(hemi)){
		if(!identical(hemi, ypsim$hemi))warning("**-Hemiphoto object and hemiphoto used for simulation are not the same!")
	}
	
	# PLANT
	if(!is.null(plant)){
	
		# Page 1. Print plant object.
		par(mfrow=c(2,2), mar=c(5,5,2,2), cex.axis=0.7)
		textplot(capture.output(print(plant,hint=FALSE)), halign="right", cex=0.9)
		viewplot(plant)
		
		# Page 2. Leaf plot.
		par(mfrow=c(1,1), xaxs="r", yaxs="r", cex.axis=0.9)
		plot(plant$ldata)
		title("Leaf")
		
		# Page 3. Plant summary.
		textplot(capture.output(summary(plant)), halign="left", cex=0.8)
		title("Plant summary", line=1)
	
		# Page 4. Leaf angle distribution.
		ang <- plant$leafdata$ang
		f <- fitdistribution(ang, "twoparbeta")
		par(yaxs="r", xaxs="r")
		plot(f, main="Leaf angle distribution with Beta fit", col="grey")
		
	}
	
	# PHYSIOLOGY
	if(!is.null(phy)){
		# Page 1.
		par(mar=c(5,5,2,2))
		textplot(capture.output(phy), cex=1, halign="center")
		title("Physiology object")
		
		# Page 2.
		# light response curve.
		# parvals <- seq(0, 1800, length=101)
		# 'lightresponse' easy.
		# 'Farquhar' : typical, midday, average met variables (Tair, VPD)?
		# Or one light response curve by timestep? (and label curves with labcurve (package Hmisc)?
		# also put histogram below the x-axis showing PARleaf distribution? (but for what timestep?)
		# timestep with highest available PAR?
	}
	
	# MET
	if(!is.null(met)){
		# Page 1.
		par(mar=c(3,3,2,2))
		textplot(capture.output(met), cex=0.7, halign="center")
		title("Weather (met) object")
		
		# Page 2.
		par(mar=c(3,3,5,2))
		plot(met)
		title("Weather (met) object", line=2)
	}
	
	# HEMI
	if(!is.null(hemi)){
		# Page 1.
		par(mar=c(3,3,2,2))
		if(!is.null(met)){
			plot(hemi, met, warn=FALSE)
		} else {
			plot(hemi)
		}
		title("Hemi photo")
	}
	
	# Yplant Simulation.
	if(!is.null(ypsim)){
		
		# Page 1.
		par(mar=c(3,3,3,2))
		textplot(capture.output(ypsim), cex=0.8, halign="center")
		title("Yplant day simulation")
		
		# Page 2. Sunlit leaf area.
		par(mar=c(5,5,2,2), cex.lab=1.3)
		plot(ypsim, type="LAsunlit", setpar=FALSE)
		
		# Page 3.
		par(xaxs="r", yaxs="r")
		plot(ypsim, openwin=FALSE)
	
	}
	
	dev.off()
	message("Report generation probably successful.")
	message("See file \'",filename,"\' in directory :\n  ",getwd())
	
}


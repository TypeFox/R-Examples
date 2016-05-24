#' Used to plot the RAD, slope and FoG parameter results

#' @description This function creates a pdf figure of plots showing the results of the imageJ analysis for resistance (RAD), tolerance (FoG), and sensitivity (slope).

#' @inheritParams plotRaw
#' @inheritParams twoParamPlot
#' @param slopeMax maximum y axis value for slope (sensitivity) plot

#' @details Basic parameter plotting functions for three parameter plots (RAD, FoG , slope). Input can be the dataframe from either \code{\link{createDataframe}} \code{type="df"} or from \code{\link{aggregateData}} \code{type=="ag"}. The default is to plot tolerance as a barplot and RAD and slope as a dotplot; tolerance can also be plotted as a dotplot with \code{barplot=FALSE} though there is currently not support to plot either RAD or slope as a barplot in this framework. 

#' @return Either a pdf figure figure saved to the 'figures' directory ("projectName_RAD-slope-FoG.pdf" or a figure on screen

#' @seealso \code{\link{oneParamPlot}} for a similar figure with one parameter and \code{\link{twoParamPlot}} for a similar figure with two parameters 

#' @export




threeParamPlot <- function(projectName, type, RAD = "RAD20", FoG = "FoG20", RADmin = 30, slopeMax = 160, tolMax = 100, width = 6, height = 4, xlabels="line", xlabAngle=NA, order=NA, orderFactor = "line", overwrite=TRUE, savePDF= TRUE, popUp = TRUE){
	
	dir.create(paste("figures/", projectName,  sep=""), showWarnings = FALSE)
	t <- file.path("figures", projectName,  paste(projectName, "_RAD-slope-FoG.pdf", sep=""))
	if (!overwrite){
		if (file.exists(t)){
			t <- file.path("figures", projectName, paste(projectName, "_RAD-slope-FoG_2.pdf", sep=""))
			if (file.exists(t)){
				k <- 2
				while(file.exists(t)){
					k <- k+1
					t <- file.path("figures", projectName, paste(projectName, "_RAD-slope-FoG_", k, ".pdf", sep=""))
					}
				}
			}
		}
	if(type == "ag" & !is.na(order[1])){
		data <- eval(parse(text=paste(projectName, ".ag", sep="")))	
		var <- substring(names(data)[length(data)], 1, 2)
		if(order[1]=="factor"){
			ordData<-data[order(data[, orderFactor]),] 
			if(length(xlabels)==1){
		 		xlabels <- as.character(ordData[, xlabels])
			}	
		}
		if(!order[1]=="factor"){
			ordData <-  data[order, ]
			if(length(xlabels)==1){
				 xlabels <- as.character(ordData[, xlabels])
			}
		}
		
	}
	if(is.na(order[1])){
		if(type=="ag"){
			ordData <- eval(parse(text=paste(projectName, ".ag", sep="")))	
			var <- substring(names(ordData)[length(ordData)], 1, 2)	
			if(length(xlabels)==1){
				 xlabels <- as.character(ordData[, xlabels])
			}
		}
		
		if(type=="df"){
			ordData <- eval(parse(text=paste(projectName, ".df", sep="")))
			if(length(xlabels)==1){
				 xlabels <- unique(as.character(ordData[, xlabels]))
			}
		}
	}

	tols <- ordData[, FoG]
	mp <- barplot(t(tols), beside=TRUE, plot=FALSE)	
	if(savePDF){
		 pdf(t, width=width, height=height)
		}	
	par(mfrow=c(3, 1), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
	if(type=="ag"){
		plot(mp[1,], ordData[, RAD], ylim=c(RADmin, 0), yaxt="n", xaxt="n", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
		arrows(mp[1,], ordData[, RAD]-ordData[, paste(var, ".", RAD, sep="")], mp[1,], ordData[, RAD]+ordData[,paste(var, ".", RAD, sep="")], length=0)
	axis(1, at=mp[1,], labels=FALSE)
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, RAD], ylim=c(RADmin, 0), yaxt="n", xaxt="n", pch=19, xlab="", ylab="", col=grey(0.3), cex=1, xlim=c(0.5, length(xlabels)+0.5))
	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
	}
	axis(2, las=2, cex.axis=1)
	mtext("Distance from\n disk (mm)", side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)

	if(type=="ag"){
		plot(mp[1,], ordData[, "slope"], ylim=c(0, slopeMax), yaxt="n", xaxt="n",  pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
		arrows(mp[1,], ordData[, "slope"]-ordData[, paste(var, ".slope", sep="")], mp[1,], ordData[, "slope"]+ordData[,paste(var, ".slope", sep="")], length=0)
	axis(1, at=mp[1,], labels=FALSE)
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, "slope"], ylim=c(0, slopeMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1, xlim=c(0.5, length(xlabels)+0.5))
	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
	}
	axis(2, las=2, cex.axis=1)

	title <- as.list(expression(paste("slope at ", RAD[50], sep="")))
	mtext(do.call(expression, title), side=2, cex=0.8, line = c(3.75,2.5))
	mtext(expression(paste(bold(B), " Sensitivity", sep="")), side=3, adj=0.01)
		
	if(type=="ag"){	
		mp <- barplot(t(tols*100), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
		box()
		 arrows(mp[1,], ordData[,FoG]*100-ordData[,paste(var, ".", FoG, sep="")]*100, mp[1,], ordData[,FoG]*100+ ordData[,paste(var, ".", FoG, sep="")]*100, length=0)
		if(is.na(xlabAngle)) 	axis(1, at=mp[1,], labels=xlabels)
		else{
			axis(1, at=mp[1,], labels=FALSE)
			text(mp[1,],  -10, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, FoG]*100, ylim=c(0, tolMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1, xlim=c(0.5, length(unique(xlabels))+0.5))
		if(is.na(xlabAngle)){
			 axis(1, at=1:length(xlabels), labels=xlabels)
			 print("here")
			 }
		else{
			axis(1, at=1:length(xlabels), labels=FALSE)
			text(1:length(xlabels),  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	axis(2, las=2, at=c(0, 20, 40, 60, 80, 100), cex.axis=1)
	mtext("Growth above\n RAD (%)",  side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(C), " Tolerance", sep="")), side=3, adj=0.01)
	if(savePDF){
		dev.off()
		cat(paste("\tFigure saved: ", t, sep=""))
		if(popUp){
			tt <- paste("open ",t)
			system(tt)
		}
	}

}
#' Used to plot the RAD and FoG results

#' @description This function creates a pdf figure of plots showing the results of the imageJ analysis for resistance (RAD) and tolerance (FoG).

#' @inheritParams plotRaw
#' @param type specify whether the dataset to use is a dataframe with all data ("df") or an aggregated dataframe ("ag")
#' @param RAD specify the RAD (radius) parameter to be plotted ("RAD20", "RAD50" or "RAD80"), default = "RAD20".
#' @param FoG specify the FoG (fraction of growth) parameter to be plotted ("FoG20", "FoG50" or "FoG80"), default = "FoG20".
#' @param RADmin minimum distance from the disk for resistance plot (minimum y axis value), default = 30.
#' @param tolMax maximum y axis value for tolerance plot. Note tolerance is coverted to a percent, default = 100.
#' @param xlabels either a vector containing the desired x-axis labels, or a single value indicating the column name that contains the values to use (likely either the 'line' column or one of the type columns), default = "line".
#' @param xlabAngle indicates whether to print the x axis labels on a angle, if a number is provided this will be the angle used. The defauilt is not to plot on an angle, default = NA.
#' @param order can be either "factor" or "custom". If custom, supply a numberial vector the same length as the dataframe to indicate the desired order. If factor, supply the column name in \code{ordeFactor} to be used to factor. 
#' @param orderFactor if \code{order = "factor"} supply the column name to be used to factor. 
#' @param barplot whether to plot tolerance as a barplot (barplot = TRUE) or dotplot (barplot = FALSE), default = TRUE. Only possible when \code{type = "ag"}

#' @details Basic parameter plotting functions to plot RAD and FoG parameter plots. Input can be the dataframe from either \code{\link{createDataframe}} \code{type="df"} or from \code{\link{aggregateData}} \code{type=="ag"}. The default is to plot RAD as a dotplot and tolerance as a barplot, though tolerance can also be plotted as a dotplot with \code{barplot=FALSE} (currently there is not support to plot RAD as a barplot in this framework). 

#' @return Either a pdf figure figure (projectName_RAD-FoG.pdf) saved to the 'figures' directory or a figure on screen

#' @seealso \code{\link{oneParamPlot}} for a similar figure with one parameter or \code{\link{threeParamPlot}} for a similar figure with three parameters 


#' @export

twoParamPlot <- function(projectName, type, RAD = "RAD20", FoG = "FoG20",  RADmin = 30, tolMax = 100, width = 6, height = 4, xlabels ="line", xlabAngle=NA, order=NA, orderFactor = "line", overwrite=TRUE, savePDF= TRUE, popUp = TRUE, barplot=TRUE){
	if(!(hasArg(type))){
		cont <- readline(paste("Please select whether dataframe is from 'createDataframe' (df) or `aggregateData (ag) ", sep=""))
		type <- cont
	}

	dir.create(paste("figures/", projectName,  sep=""), showWarnings = FALSE)
	t <- file.path("figures", projectName,  paste(projectName, "_RAD-FoG-", type, ".pdf", sep=""))
	if (!overwrite){
		if (file.exists(t)){
			t <- file.path("figures", projectName, paste(projectName, "_RAD-FoG_2-", type, ".pdf", sep=""))
			if (file.exists(t)){
				k <- 2
				while(file.exists(t)){
					k <- k+1
					t <- file.path("figures", projectName, paste(projectName, "_RAD-FoG_", k, "-", type, ".pdf", sep=""))
					}
				}
			}
		}
	
	if(type == "ag"){
		if(!is.na(order[1])){
		data <- eval(parse(text=paste(projectName, ".ag", sep="")))	
		var <- substring(names(data)[length(data)-2], 1, 2)
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
			ordData <- eval(parse(text=paste(projectName, ".ag", sep="")))	
			var <- substring(names(ordData)[length(ordData)-2], 1, 2)	
			if(length(xlabels)==1){
				 xlabels <- as.character(ordData[, xlabels])
			}
		}
	}	
	if(type=="df"){
		if(!is.na(order[1])){
			data <- eval(parse(text=paste(projectName, ".df", sep="")))			
			if(order[1]=="factor"){
				ordData<-data[order(data[, orderFactor]),] 
				if(length(xlabels)==1){
			 		xlabels <- unique(as.character(ordData[, xlabels]))
				}	
			}
			if(!order[1]=="factor"){
				ordData <-  data[order, ]
				if(length(xlabels)==1){
					 xlabels <- unique(as.character(ordData[, xlabels]))
				}
			}
		}
		if(is.na(order[1])){
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

	par(mfrow=c(2, 1), oma=c(4, 4, 1, 1), mar=c(1, 1, 1, 1))
	if(type=="ag"){
		if(barplot == TRUE){		
			plot(mp[1,], ordData[, RAD], ylim=c(RADmin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
			arrows(mp[1,], ordData[, RAD]-ordData[, paste(var, ".", RAD, sep="")], mp[1,], ordData[, RAD]+ordData[,paste(var, ".", RAD, sep="")], length=0)
		axis(1, at=mp[1,], labels=FALSE)
		}
	else{
			plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, RAD], ylim=c(RADmin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), 	xlim=c(0, max(mp)+1), cex=1.4)
			arrows(as.numeric(as.factor(ordData[, orderFactor])), ordData[, RAD]-ordData[, paste(var, ".", RAD, sep="")], as.numeric(as.factor(ordData[, orderFactor])), ordData[, RAD]+ordData[,paste(var, ".", RAD, sep="")], length=0)
		axis(1, at=mp[1,], labels=FALSE)
		}
		
	}
	
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, RAD], ylim=c(RADmin, 0), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1, xlim=c(0.5, length(unique(as.numeric(as.factor(ordData[, orderFactor]))))+0.5))
	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
	}
	
	axis(2, las=2, cex.axis=0.8)
	mtext("Distance\n from disk (mm)", side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(A), " Resistance", sep="")), side=3, adj=0.01)
	
	if(type=="ag"){
		if(barplot == TRUE){	
			mp <- barplot(t(tols*100), ann=FALSE, beside=TRUE, yaxs="i", xaxs="i", ylim=c(0, tolMax), xaxt="n", yaxt="n", xlab="", ylab="", xlim=c(0, max(mp)+1))
			box()
		 	arrows(mp[1,], ordData[,FoG]*100-ordData[, paste(var, ".", FoG, sep="")]*100, mp[1,], ordData[,FoG]*100+ ordData[,paste(var, ".", FoG, sep="")]*100, length=0)
			if(is.na(xlabAngle)) 	axis(1, at=mp[1,], labels=xlabels)
			else{
				axis(1, at=mp[1,], labels=FALSE)
				text(mp[1,],  -10, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
			}
		}
		else{
			plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, FoG], ylim=c(0, tolMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1.4, xlim=c(0.5, length(xlabels)+0.5))
			axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
			arrows(as.numeric(as.factor(ordData[, orderFactor])), ordData[,FoG]*100-ordData[, paste(var, ".", FoG, sep="")]*100, as.numeric(as.factor(ordData[, orderFactor])), ordData[,FoG]*100+ ordData[,paste(var, ".", FoG, sep="")]*100, length=0)
			if(is.na(xlabAngle)) 	axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=xlabels)
			else{
				axis(1, at=as.numeric(as.factor(unique(ordData[, orderFactor]))), labels=FALSE)
				text(as.numeric(as.factor(unique(ordData[, orderFactor]))),  -10, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
			}
		 }
	}
	if(type=="df"){
		plot(as.numeric(as.factor(ordData[, orderFactor])), ordData[, FoG]*100, ylim=c(0, tolMax), yaxt="n", xaxt="n", yaxs="i", xaxs="i", pch=19, xlab="", ylab="", col=grey(0.3), cex=1, xlim=c(0.5, length(unique(as.numeric(as.factor(ordData[, orderFactor]))))+0.5))
		if(is.na(xlabAngle)){
			 axis(1, at=1:length(xlabels), labels=xlabels)
			 }
		else{
			axis(1, at=1:length(xlabels), labels=FALSE)
			text(1:length(xlabels),  -5, xlabels, srt = xlabAngle, xpd=NA, adj=0, cex=0.8)
		}
	}
	axis(2, las=2, at=c(0, 20, 40, 60, 80, 100), cex.axis=0.8)
	mtext("Growth\n above RAD (%)",  side=2, line=2.5, cex=0.8)
	mtext(expression(paste(bold(B), " Tolerance", sep="")), side=3, adj=0.01)
	if(savePDF){
		dev.off()
		cat(paste("\tFigure saved: ", t, sep=""))
		if(popUp){
			tt <- paste("open ",t)
			system(tt)
		}
	}

}


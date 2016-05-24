#Author: Michael Sieger <sieger-michael@web.de>
#Project responsible: Dr. Georg Ohmayer <georg.ohmayer@hswt.de>
#Copyrights: Hochschule Weihenstephan-Triesdorf

#This file is part of the Linkage Maps package.

#The VisuClust package is free software: you can redistribute it and/or modify
#it under the terms of the GNU Lesser General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#The VisuClust package is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU Lesser General Public License for more details.

#You can find a copy of the GNU Lesser General Public License at <http://www.r-project.org/Licenses/> or <http://www.gnu.org/licenses/>.




LinkageMap <- function(xSammon, dist, lineTypes=c("solid","dotted", "dashed"), lineColors=c("red","green","blue"), 
			lineWidths=c(1,1,1), labels = NULL, cluster = NULL, maxValue=0.5, legendDigits = 2, xlab = "", 
			ylab = "", main = "")
	{	 

	# checking input
	xdim = dim(xSammon)
	llen = length(labels)
	clen = length(cluster)
	if(length(xdim) != 2 || xdim[2] != 2)
	{
		print("The dimension of xSammon must be [n, 2]")
	}
	if(length(lineTypes) != length(lineColors) || length(lineTypes) != length(lineWidths))
	{
		print("lineTypes, lineColors, and lineWidths doesnt contain the same number of elements")
	}
	if(llen != 0 && llen != xdim[1])
	{
		print("The dimension of the label field is wrong")
	}
	if(clen != 0 && clen != xdim[1])
	{
		print("The dimension of the cluster field is wrong")
	}
	if(maxValue <= 0 || maxValue > 1)
	{
		print("maxValue is not between 0 and 1")
	}	
	D = as.matrix(dist)

	#defining the variables to prevent them beeing global visible
	maps.s = 1
	maps.nLines = 1
	maps.clusterCols = 1
	maps.circleCols = 1
	maps.n = 1
	maps.nCluster = 1
	maps.xRange = 1
	maps.yRange = 1
	
	#functions
	maps.format <- function(number){
		round(number, legendDigits)
	}

	maps.createLegendString <- function(v1, v2){
		paste(v1, "<= d[i,j] <", v2)
	}

	maps.getSliderValues <- function(){
		sliderValues = rep(NA, maps.nLines)
		for(i in 1:maps.nLines){
			sliderValues[i] = slider(no=i)
		}
		sliderValues
	}

	maps.updateGUI <- function(){
		
		sliderValues = maps.getSliderValues()
		moved = (sliderValues != maps.s)
		for(i in 1:maps.nLines)
		{
			if(moved[i])
			{
				if(i != 1)
				{
					for(j in (i-1):maps.nLines)
					{
						if(slider(no=j) > slider(no=i))
						{
							slider(set.no.value=c(j, slider(no=i)))
						}
					}
				}
				if(i != maps.nLines)
				{
					for(j in (i+1):maps.nLines)
					{
						if(slider(no=j) < slider(no=i))
						{
							slider(set.no.value=c(j, slider(no=i)))
						}
					}
				}

				break;
			}
		}
		maps.s <<- sliderValues
	}

	maps.createCircleColors <- function(colors){
		colLen <- length(colors)
		revCol <- rep("black", colLen)
		#for(i in 1:colLen){
		#	revCol[(i+colLen/2)%%colLen] <- colors[i]
		#}
		revCol
	}

	maps.drawLegend <- function(){
		t <- rep(NA, maps.nLines)
		t[1] = paste("d[i,j] <= ", maps.s[1])
		for(i in 2:maps.nLines){
			t[i] = maps.createLegendString(maps.s[i-1], maps.s[i])
		}
		#legend("topleft",  
		#		t, 
		#		lwd=lineWidths, lty=lineTypes, col=lineColors)
		psize <- par("usr")
		lsize <- legend(0,0,t, horiz=TRUE, plot=FALSE)
		legend(psize[1], psize[4]+lsize$rect$h, t, horiz=TRUE, lwd=lineWidths, lty=lineTypes, col=lineColors)
	}

	maps.draw <- function(){
		maps.updateGUI()
		dev.hold()
		par(xpd=TRUE, ask=FALSE)
		plot(xSammon, pch=21,col=maps.circleCols[cluster], bg = maps.clusterCols[cluster], xlab=xlab, ylab=ylab, main=main)
		for(j in 1:maps.n){		#foreach point
			for(i in 1:j){		#foreach point (only one direction)
			      if(i != j)
			      {
				    for(k in 1:maps.nLines){
					index = maps.nLines-k+1
					if((index == 1 && D[i,j] < maps.s[index]) || (index != 1 && D[i,j] < maps.s[index] && D[i,j] >= maps.s[index-1]))
					{
						segments(xSammon[i, 1], xSammon[i,2], xSammon[j,1], xSammon[j,2], lty=lineTypes[index], col=lineColors[index], lwd=lineWidths[index])
						break;
					}
				    }
			      }
			}
		}
		if(length(labels)!=0)
		{
			relDist = 0.02
			text(xSammon[,1]+maps.xRange*relDist,xSammon[,2]+maps.yRange*relDist,labels)
		}
		maps.drawLegend()
		dev.flush()
	}

	maps.sliderCallback <- function(...){
		maps.draw()
	}

	maps.createSliderNames <- function(){
		res = rep(NA, maps.nLines)
		for(i in 1:maps.nLines){
			res[i] = paste("t[", i, "]")
		}
		res
	}

	maps.init <- function()
	{
		maps.nLines <<- length(lineTypes)
		maps.n <<- length(xSammon[,1])
		if(length(cluster) == 0)
		{
			cluster <<- rep(1,maps.n)
			maps.clusterCols <<- rgb(0,0,0,0)
			maps.nCluster <<- 1
		}
		else
		{
			maps.nCluster <<- length(unique(cluster))
			maps.clusterCols <<- rainbow(maps.nCluster)
		}
		maps.circleCols <<- maps.createCircleColors(maps.clusterCols)
		min <- min(D[D !=0])*0.9	#below smallest nonzero distance
		max <- max(D)*maxValue
		maps.s <<- rep(min, maps.nLines)
		maps.xRange <<- max(xSammon[,1])-min(xSammon[,1])
		maps.yRange <<- max(xSammon[,2])-min(xSammon[,2])

		slider(maps.sliderCallback, 
				sl.names=maps.createSliderNames(),
				sl.mins=rep(min,maps.nLines*2),
				sl.maxs=rep(max,maps.nLines),
				sl.deltas=rep(0.01,maps.nLines), 
				sl.defaults=rep(0,maps.nLines),
				prompt=T,	
				title="Thresholds")
	}

	#start it
	maps.init()
}





# 		PCA tool: provides additional functionality to R's 'base' package 'svd' 
#                 method, part of the 'MetStaT' package  
#		Copyright (C) 2012 Tim Dorscheidt
#		
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.
#		
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#		
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# 		Email: g.zwanenburg@uva.nl ('MetStaT' contact person) or tdorscheidt@gmail.com
###############################################################################

# use R's build in 'svd' method from the 'base' package to calculate the svd, and also add scores and the percentage of variance explained per component
PCA.Calculate <- function(data) {

	# performs a standard value decomposition (pca) on the data
	svd.result <- svd(data)
	svd.result$var.explained <- svd.result$d^2
	svd.result$var.explained <- svd.result$var.explained/(sum(svd.result$var.explained)) 
	svd.result$t <- svd.result$u %*% diag(svd.result$d)
	svd.result$u <- NULL
	svd.result
}

# plot scores of PCA, wehereby points can be distinguished by different colors and/or symbol-types
PCA.PlotScores <- function(pr.object, pcs = c(1,2), labels = "none", custom.labels = NULL, dot.class.vector = NULL, col.class.vector = NULL) {
	prev.mfrow <- par("mfrow") # save current plot-window dimensions to restore this back to default after we're done
	prev.xpd <- par("xpd")# save current plot-window clipping parameter to restore this back to default after we're done
	prev.oma <- par("oma")# save current plot-window margin parameters to restore this back to default after we're done
	if (is.character(pcs)) eval(parse(text=paste(sep="","pcs <- c(",pcs,")")))
	list.of.pc.tuples <- MetStaT.GetPcTuples(length(pcs)) # get tuples of PC-indices that are plotted versus one another
	par(mfrow=c(ceiling(length(list.of.pc.tuples)/ceiling(length(list.of.pc.tuples)^0.5)),ceiling(length(list.of.pc.tuples)^0.5))) # determine dimensions of plot-window
	par(xpd=NA) # diasable clipping
	par(oma=c(0,0,0,4)) # increase right margin for legend plotting
	no.points <- length(pr.object$d)
	tuple.index <- 0
	first.warning <- TRUE
	for (tuple in list.of.pc.tuples) {
		tuple.index <- tuple.index + 1
		pc1 <- pcs[tuple[1]] # 1st element of tuple contains index of PC1 within PCs
		pc2 <- pcs[tuple[2]] # 2nd element of tuple contains index of PC2 within PCs
		main.title <- paste("Scoreplot PC",pc1," vs PC",pc2,sep="") # title contains names of PCs
		plot.type = 'p' # set default plot-type, which can be overwritten by the following switch for label-types
		labels.to.plot <- NULL
		switch(pmatch(labels, c("dots","none","numerical","custom")), # depending on value of 'labels', set plot label variables
				{plot.type = 'p'},{plot.type = 'p'}, # no labels, just dots
				{plot.type = 'n'; labels.to.plot<-1:no.points}, # numerical labels
				{ # custom labels
					plot.type = 'n'; 
					if (!is.null(custom.labels)) labels.to.plot<-unlist(strsplit(as.character(custom.labels),split=","));
					if (length(labels.to.plot)==0){
						if ( first.warning ) warning("No custom labels defined. Using numerical labels.")
						labels.to.plot<-1:no.points; first.warning <- FALSE
					}
					if (no.points>length(labels.to.plot)){
						if ( first.warning ) warning("No. points to plot exceeds no. custom labels defined. Re-using from start.");  first.warning <- FALSE
					}
					
					if (length(labels.to.plot>no.points)){labels.to.plot <- labels.to.plot[1:no.points]} 
				} 
		)
		label.offset <- FALSE; # by default, labels have no offset
		dot.types.to.plot <- 1 # by default, all dots are of type 1
		if (!is.null(dot.class.vector)) {
			plot.type = 'p'; label.offset <- TRUE; # in case dot-plotting was overwritten with custom labels, re-activate plotting of dots and draw labels with offset
			if (length(dot.class.vector)==1) dot.class.vector <- unlist(strsplit(dot.class.vector,split=","))
			dot.types.to.plot <- MetStaT.ConvertToNumericClasses(dot.class.vector,new.classes=1:25) # the number of classes plottable with dot-types cannot exceed 25; will loop otherwise
			dot.types.to.plot <<- dot.types.to.plot
		}
		col.types.to.plot <- 1 # by default, everything's plotted in black
		if (!is.null(col.class.vector)) {
			if (length(col.class.vector)==1) col.class.vector <- unlist(strsplit(col.class.vector,split=","))
			col.types.to.plot <- MetStaT.ConvertToNumericClasses(col.class.vector,new.classes=1:8) # the number of colors plottable with dot-types cannot exceed 8; will loop otherwise
			col.types.to.plot <<- col.types.to.plot
		}
		plot(pr.object$t[,pc1],pr.object$t[,pc2], # prepare plot window
				main=main.title,
				xlab=paste("PC",pc1," (",formatC(pr.object$var.explained[pc1] * 100,digits=2,format="f"),"%)",sep=""),
				ylab=paste("PC",pc2," (",formatC(pr.object$var.explained[pc2] * 100,digits=2,format="f"),"%)",sep=""),
				type = plot.type,
				pch = dot.types.to.plot, 
				col=col.types.to.plot);
		if (!is.null(labels.to.plot)) { # if labels were requested, plot these
			text(pr.object$t[,pc1],pr.object$t[,pc2],
					labels=labels.to.plot,
					cex=0.8,
					col=col.types.to.plot,
					if (label.offset) {offset=-0.3})
		}
		if (tuple.index==ceiling(length(list.of.pc.tuples)^0.5)) { # when plotting the first plot that's completely on the right
			legend.names <- c(); legendColors <- c(); legendPch <- c(); legendLty <- c(); legendLwd <- c()
			if (!is.null(dot.class.vector)) { # if dot-types are classifiers, then plot dot type legend
				noOfLegendItems <- length(unique(dot.class.vector))
				legend.names <- unique(dot.class.vector)
				legendColors <- rep(1,noOfLegendItems)
				legendPch <- rep(1:25,10)[1:noOfLegendItems]
				legendLty <- rep(NA,noOfLegendItems)
				legendLwd <- rep(1,noOfLegendItems)
			}
			if (!is.null(col.class.vector)) { # if color-types are classifiers, then plot color type legend
				noOfLegendItems <- length(unique(col.class.vector))
				legend.names <- c(legend.names,unique(col.class.vector))
				legendColors <- c(legendColors,rep(1:8,10)[1:noOfLegendItems])
				legendPch <- c(legendPch,rep(NA,noOfLegendItems))
				legendLty <- c(legendLty,rep(1,noOfLegendItems))
				legendLwd <- c(legendLwd,rep(4,noOfLegendItems))
			}
			if (!is.null(legend.names)) {
				legend(x = max(pr.object$t[,pc1])+0.05*(max(pr.object$t[,pc1])-min(pr.object$t[,pc1])), y=max(pr.object$t[,pc2]), title = "Legend:",
						legend=legend.names,
						cex=0.8, 
						col=legendColors, 
						pch=legendPch,
						lty=legendLty,
						lwd=legendLwd);
			}
		}
	}
	par(mfrow=prev.mfrow) # restore plot window dimensions back to previous
	par(xpd=prev.xpd) # restore clipping parameter back to previous
	par(oma=prev.oma) # restore margin parameters back to previous
}

# plot loadings of PCA
PCA.PlotLoadings <- function(pr.object, pcs = c(1,2)) {
	if (is.character(pcs)) eval(parse(text=paste(sep="","pcs <- c(",pcs,")")))
	list.of.pc.tuples <- MetStaT.GetPcTuples(length(pcs)) # get tuples of PC-indices that are plotted versus one another
	tuple.index <- 0
	prev.mfrow <- par("mfrow") # save current plot-window dimensions to restore this back to default after we're done
	prev.xpd <- par("xpd")# save current plot-window clipping parameter to restore this back to default after we're done
	prev.oma <- par("oma")# save current plot-figure margin parameters to restore this back to default after we're done
	prev.mar <- par("mar")# save current plot-area margin parameters to restore this back to default after we're done
	par(oma=c(0,0,0,0)) # reduce figure margins to maximize drawing area 
	par(mfrow=c(2,1)) # we'll be plotting two seperate loading-plots above one another
	par(xpd=NA) # diasable clipping
	for (tuple in list.of.pc.tuples) {
		tuple.index <- tuple.index + 1
		pc1 <- pcs[tuple[1]] # 1st element of tuple contains index of PC1 within PCs
		pc2 <- pcs[tuple[2]] # 2nd element of tuple contains index of PC2 within PCs
		
		# top plot
		par(mar=c(0,4.1,4.1,2.1)) # x-axis will only contain ticks, remove lower margins to the point where plots are touching
		main.title <- paste("Loadings PC",pc1," and PC",pc2,sep="") # title contains names of PCs
		ylab.title <- paste("Distance to PC",pc1,sep="") # y axis label
		color.to.Use <- pc1
		plot(1:dim(pr.object$v)[1],pr.object$v[,pc1],type="h",col=color.to.Use,
				main=main.title,
				axes=FALSE, # manually make axes
				xlab="",
				ylab=ylab.title,
				ylim= c(
						min(c(pr.object$v[,pc1],pr.object$v[,pc2])),
						max(c(pr.object$v[,pc1],pr.object$v[,pc2])) 
				) # limit y-axis to min and max of either pc
		)
		axis(1, labels=FALSE)
		axis(2, labels=TRUE)
		
		# bottom plot
		par(mar=c(5.1,4.1,0,2.1)) # remove upper margins to the point where plots are touching
		ylab.title <- paste("Distance to PC",pc2,sep="") # y axis label
		color.to.Use <- pc2
		plot(1:dim(pr.object$v)[1],pr.object$v[,pc2],type="h",col=color.to.Use,
				main=NULL,
				axes=FALSE, # manually make axes
				xlab="Variables",
				ylab=ylab.title,
				ylim= c(
						min(c(pr.object$v[,pc1],pr.object$v[,pc2])),
						max(c(pr.object$v[,pc1],pr.object$v[,pc2])) 
				) # limit y-axis to min and max of either pc
		)
		axis(1, labels=TRUE)
		axis(2, labels=TRUE)
	}
	par(mfrow=prev.mfrow) # restore plot window dimensions back to previous
	par(xpd=prev.xpd) # restore clipping parameter back to previous
	par(oma=prev.oma) # restore figure margin parameters back to previous
	par(mar=prev.mar) # restore plot margin parameters back to previous
}
separationplot <-
function(pred, actual, type="line", line=T, lwd1=0.5, lwd2=0.5, heading="", xlab="", shuffle=T, width=9, height=1.2, col0="#FEF0D9", col1="#E34A33", flag=NULL, flagcol=1, file=NULL, newplot=T, locate=NULL, rectborder=NA, show.expected=F, zerosfirst=T, BW=F){

	# do some error-testing first:
	if (is.vector(pred)==F) stop("The pred argument needs to be a vector")
	if (is.vector(actual)==F) stop("The actual argument needs to be a vector")
	if (length(pred)!=length(actual)) stop("The pred and actual vectors are of different lengths.")
	if (any(is.na(pred))) stop("Missing values in the pred vector.")
	if (any(is.na(actual))) stop("Missing values in the actual vector.")
	
	resultsmatrix<-data.frame(pred, actual, flags=0)
	rows<-nrow(resultsmatrix)
	
	if (!is.null(flag)) resultsmatrix$flags[flag]<-1
	
	if (shuffle==T){
	set.seed(1)
	resultsmatrix<-resultsmatrix[sample(1:rows,rows),]
		} # close if shuffle==T condition
	
	# sort the results matrix
	resultsmatrix<-resultsmatrix[order(resultsmatrix$pred),]
	resultsmatrix<-cbind(resultsmatrix, position=1:rows)
	
	# warn user if there are relatively few discrete values of pred:
	if (shuffle==F & length(unique(pred))<0.9*length(pred)) warning("Your pred vector contains one or more long runs of identical values.  The separation plot needs to be interpreted with care.  Alternatively, try setting shuffle=T.")
	
	if (type=="bands"){width<-6; height<-2} # note that this overrides the user's options for now.
	
	# implement the black & white color scheme:
	if (BW==T){col0="#F0F0F0"; col1="#636363"}
	
	# open the plot space only if newplot==T:
	if (newplot==T){
		if (is.null(file)) dev.new(width=width, height=height)
		if (!is.null(file)) pdf(file=file, width=width, height=height)
		par(mgp=c(3,0,0), lend=2, mar=c(2,2,2,2))	
		} # close newplot==T condition
	
	# set up the plot space:
	if (type!="bands"){	 			
		plot(1:nrow(resultsmatrix),1:nrow(resultsmatrix), xlim=c(0.5, nrow(resultsmatrix)+0.5), ylim=c(-0.1,1),  type="n", bty="n", yaxt="n", xaxt="n", xlab=xlab, ylab="")
		title(main=heading)
		} # close type!="bands"
		
	resultsmatrix$color<-NA
	resultsmatrix$color[resultsmatrix$actual==1]<-col1
	resultsmatrix$color[resultsmatrix$actual==0]<-col0

	events<-resultsmatrix[resultsmatrix$actual==1,]
	nonevents<-resultsmatrix[resultsmatrix$actual==0,]
	
	# add the line segments as lines with the 1s plotted on top of the 0s:
	if (type=="line" & zerosfirst==T){
		if (nrow(nonevents)>0) segments(x0=nonevents$position, x1=nonevents$position, y0=0, y1=1, col=col0, lwd=lwd1)
		if (nrow(events)>0) segments(x0=events$position, x1=events$position, y0=0, y1=1, col=col1, lwd=lwd1)
		# add flags:
		if (!is.null(flag)) segments(x0=resultsmatrix$position[resultsmatrix$flags==1], x1=resultsmatrix$position[resultsmatrix$flags==1], y0=0, y1=1, col=flagcol, lwd=lwd1)
		}

	# add the line segments as lines with the 0s plotted on top of the 1s:
	if (type=="line" & zerosfirst==F){
		if (nrow(events)>0) segments(x0=events$position, x1=events$position, y0=0, y1=1, col=col1, lwd=lwd1)
		if (nrow(nonevents)>0) segments(x0=nonevents$position, x1=nonevents$position, y0=0, y1=1, col=col0, lwd=lwd1)

		# add flags:
		if (!is.null(flag)) segments(x0=resultsmatrix$position[resultsmatrix$flags==1], x1=resultsmatrix$position[resultsmatrix$flags==1], y0=0, y1=1, col=flagcol, lwd=lwd1)
		}
	
	# add the line segments as rectangles:
	if (type=="rect") {
		
		rect(xleft=resultsmatrix$position-0.5, ybottom=0, xright=resultsmatrix$position+0.5, ytop=1, col=resultsmatrix$color, border=rectborder)
		# add flags
		if (!is.null(flag)) rect(xleft=resultsmatrix$position[resultsmatrix$flags==1]-0.5, xright=resultsmatrix$position[resultsmatrix$flags==1]+0.5, ybottom=0, ytop=1, col=flagcol,  border=rectborder)
	
		} # close type=="rect" condition

	
	# Calculate the expected number of events:
				
		# sum under the phat curve:
		expectedevents<-round(sum(resultsmatrix$pred))
				
		# add a marker on the separation plot showing the expected number of events:		
		if (show.expected) points(nrow(resultsmatrix)-expectedevents+0.5, -0.1, pch=24, bg=1, cex=0.7)
		
		# calculate the new cutpoint:
		newcutpoint<-sort(resultsmatrix$pred, decreasing=T)[expectedevents]
		
		# calculate the proportion correctly predicted:
		tp<-length(resultsmatrix$actual[resultsmatrix$pred>=newcutpoint & resultsmatrix$actual==1])
		fp<-length(resultsmatrix$actual[resultsmatrix$pred>=newcutpoint & resultsmatrix$actual==0])
		tn<-length(resultsmatrix$actual[resultsmatrix$pred<newcutpoint & resultsmatrix$actual==0])
		fn<-length(resultsmatrix$actual[resultsmatrix$pred<newcutpoint & resultsmatrix$actual==1])
		
		pcp<-(tp+tn)/length(resultsmatrix$actual)
	
	
	# Try a different type of plot when n is very large:
		
	if (type=="localaverage"){
		
		if (nrow(resultsmatrix)>5000) cat("\nCalculating the moving averages.  This may take a few moments due to the large number of observations.\n")
		
		windowsize<-round(nrow(resultsmatrix)*0.01)
		resultsmatrix$localaverage<-NA
		for (i in 1:nrow(resultsmatrix)){
			lower<-max(c(1, i-windowsize))
			upper<-min(c(nrow(resultsmatrix), i+windowsize))
			resultsmatrix$localaverage[i]<-mean(resultsmatrix$actual[lower:upper])
			
			} # close i loop
		
		lines(1:rows, resultsmatrix$localaverage, lwd=lwd1, col=col1)
		
		} # close local average condition


	# add a line showing the predicted probabilities
	if (line==T & type!="bands")
	lines(1:rows, resultsmatrix$pred, lwd=lwd2)

	
	# bands plot:
	if (type=="bands"){
		
		breaks<-seq(0,0.9,0.1)
		cols<-RColorBrewer::brewer.pal(9,"Reds")
		a<-colorRampPalette(cols)
		cols<-a(10)
		
		phat.events<-events[,1]
		phat.nonevents<-nonevents[,1]	
		
		par(mgp=c(3,0,0), lend=2, mar=c(2.5,2,2.5,2))
		
		layout(matrix(c(1,2,1,2,1,2,1,2,3,3), nrow=2, ncol=5))

		plot(1:length(phat.events),1:length(phat.events), xlim=c(0.5, length(phat.events)+0.5), ylim=c(0,1),  type="n", bty="n", yaxt="n", xaxt="n", xlab="", ylab="")
		
		title(main=paste("y=1 (n=", length(phat.events), ")", sep=""))
		segments(x0=1:length(phat.events), x1=1:length(phat.events), y0=0, y1=1, col=cols[findInterval(phat.events,breaks)])
		
		plot(1:length(phat.nonevents),1:length(phat.nonevents), xlim=c(0.5, length(phat.nonevents)+0.5), ylim=c(0,1),  type="n", bty="n", yaxt="n", xaxt="n", xlab="", ylab="")
		title(main=paste("y=0 (n=", length(phat.nonevents), ")", sep=""))
		segments(x0=1:length(phat.nonevents), x1=1:length(phat.nonevents), y0=0, y1=1, col=cols[findInterval(phat.nonevents,breaks)])	
		
		# add the legend:
		plot.new()
		par(mar=c(1,1,1,1))
		legend("center", legend=c("over 0.9", "0.8 - 0.9", "0.7 - 0.8", "0.6 - 0.7", "0.5 - 0.6", "0.4 - 0.5", "0.3 - 0.4", "0.2 - 0.3", "0.1 - 0.2", "under 0.1"), fill=rev(cols), title="Probabilities:", cex=1.25, bty="n")      
		
		} # close type==bands condition
	
	# close pdf device:
	if (!is.null(file)) dev.off()
	
	# locate the points:
	if (!is.null(locate)) {
		a<-locator(n=locate)
		resultsmatrix[round(a$x),]
		} # close locate condition
	
	# return the resultsmatrix
	invisible(resultsmatrix)
	} # close separationplot function


plot.importance <-
function(x, type=c('dotchart', 'bar'), horiz=TRUE, ...){
	if(class(x)!="importance")
		stop("Wrong class")

	labels <- names(x)
	tmpx <- as.numeric(x); names(tmpx) <- labels

	if(type == 'dotchart')
		dotchart(tmpx, ...)
	
	if(type == 'bar'){

		if(horiz){
			barplot(tmpx, horiz=TRUE, las=1, ...)
		}else{
			mp <- barplot(tmpx, axes=FALSE, axisnames=FALSE, ...)
			text(mp, par("usr")[3], labels=labels, srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=.9)
			axis(2)
		}
	}
}

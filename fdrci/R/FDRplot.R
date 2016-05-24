FDRplot <-
function(plotdat,lowerbound,upperbound,mn,lpos = "bottomleft",outfile=FALSE){
	# plotdat is a table that is returned from fdrTbl()
	# pname contains a string that is the name of the permutation p-value column
	# lowerbound and upperbound are -log10(p-value) bounds for the x-axis of the plot
	# mn = main title

	##### myfdrplot.pdf
	if(outfile != FALSE) pdf(outfile)
	if(sum(!is.na(plotdat[,"ul"]))) yul = max(plotdat[,"ul"],na.rm=TRUE)
	if(abs(yul) > 1) yul = 1
	tpos = yul * 3 / 4
	plot(plotdat$threshold,plotdat$fdr,xlim=c(lowerbound,upperbound),ylim=c(0,yul),
		col='blue',type='l',lty=1,main=mn,
		xlab=expression(paste('p-value cutoff (-',log[10],'(p-value))')),
		ylab=expression(hat(FDR)))
	points(plotdat$threshold,plotdat$ul,col='red',type='l',lty=2)
	points(plotdat$threshold,plotdat$ll,col='red',type='l',lty=2)
	legend(lpos,legend=c(expression(hat(FDR)),"Upper Limit","Lower Limit"),
		col=c("blue","red","red"),cex=.8,lty=c(1,2,2))
	#text(3.5,.37,"FDR = .11, CI = (0.07,0.17)",cex=1.1,col='brown',font=2)
	text(plotdat$threshold,tpos,plotdat$S,cex=.8,col='green',srt=90)
	if(outfile != FALSE) dev.off()
	return("")
} # End plot function


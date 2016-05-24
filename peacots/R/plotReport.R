plotReport <-
function(name="", times=NULL, signal=NULL, report=NULL, plotFile=NULL, dataFile=NULL, sep=" "){
	Nplots = sum(c(!is.null(times) && !is.null(signal), !is.null(report) && !report$error));
	if(Nplots==0) return;
	
	old.par <- par(mfrow=c(1, Nplots))	
	
	if(!is.null(times) && !is.null(signal)){
		#plot time series
		plot(ts(times), ts(signal), xy.label=FALSE, type="l", ylab="signal", xlab="time", main=paste0("Time series", if(name=="") "" else paste0(" (",name,")")), cex=0.8, cex.main=0.9)
	}
	
	if(!is.null(report) && !report$error){
		#plot periodogram
		plot(ts(report$frequencies), ts(report$periodogram), xy.label=FALSE, type="l", ylab="power", xlab="frequency", main=sprintf("Periodogram OUSS analysis%s\n(peak freq=%.3g, P=%.2g, Plocal=%.2g)", (if(name=="") "" else paste0("\n",name)), report$frequencies[report$peakMode],report$P,report$Plocal), col="black", cex=0.8, cex.main=0.9);
	
		#plot fitted OUSS periodogram
		lines(report$frequencies[report$minFitMode:length(report$frequencies)], report$fittedPS[report$minFitMode:length(report$fittedPS)], col="red");
	
		#plot legend
		legend((0.6*report$frequencies[1]+0.4*tail(report$frequencies,1)), (0.85*max(report$periodogram)), c("periodogram", "fitted OUSS"), lty=c(1,1), col=c("black", "red"), bty="n", cex=0.8)

	}
	
	par(old.par)
	
	if(!is.null(plotFile)){
		#save plot as PDF
		dir.create(dirname(plotFile), showWarnings=FALSE);
		dev.copy2pdf(file=plotFile);
	}
	if(!is.null(dataFile)){
		# write data to file
		dir.create(dirname(dataFile), showWarnings=FALSE);
		cat(sprintf("# OUSS periodogram analysis of times series %s\n", name), file=dataFile, append=FALSE);
		if(!is.null(times) && !is.null(signal)){ 
			cat(sprintf("# times%ssignal\n",sep), file=dataFile, append=FALSE);
			write.table(data.frame(times, signal), , file=dataFile, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=sep);
		}
		if(!is.null(report) && !report$error){
			cat(sprintf("\n\n# periodogram\n# frequency%spower\n",sep), file=dataFile, append=TRUE);
			write.table(data.frame(report$frequencies, report$periodogram), , file=dataFile, append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=sep);
		}
	}
}

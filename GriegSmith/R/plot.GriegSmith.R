plot.GriegSmith <-
function(x,main, ...){

	if(missing(main)) {main <- "Orientation Averaged Grieg-Smith MSr"}

	plot(x[,3],ylim=c(max(min(x[,4])-.1,0),max(max(x[,5])+.1,max(x[,3])+.1)),xaxt='n',main=main,type="b",xlab="Block size",ylab="MSr")
	lines(x[,4],lty=2)
	lines(x[,5],lty=2)
	axis(1,at=1:length(x[,3]),labels=rownames(x))



	
}


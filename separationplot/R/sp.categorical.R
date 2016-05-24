sp.categorical <-
function(pred, actual, file=NULL, cex=1.5, ...){
	
	# first determine the number of categories:
	ncat<-ncol(pred)
	
	# open the plot space
	width<-9; height<-ncat
	if (is.null(file)) dev.new(width=width, height=height)
	if (!is.null(file)) pdf(file=file, width=width, height=height)
	par(mgp=c(0,0,0), lend=2, mar=c(0,0,0,0))
	layout(mat=matrix(data=1:(ncat*2), nrow=ncat, ncol=2, byrow=T), widths=c(7,1))
	
	
	# now loop through the categories:
	for (i in 1:ncat){
		pred.current<-pred[,i] # choose the ith column
		actual.current<-as.numeric(actual==colnames(pred)[i]) # find the matches in the actual vector for category i
		separationplot(pred=pred.current, actual=actual.current, newplot=F, ...)
		
		# add the category name in the next plot space:
		plot(1:100,1:100, type="n", bty="n", xaxt="n", yaxt="n")
		text(0, 50, colnames(pred)[i], adj=0, cex=cex)

		}# close i loop
	
	# close pdf device:
	if (!is.null(file)) dev.off()
	
	}# close sp.categorical function


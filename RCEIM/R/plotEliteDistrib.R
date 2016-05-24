plotEliteDistrib <- function(elite) {
# 
# Plot distribution of parameters from elite members at the solution 
# 
	nParams <- dim(elite)[2]-1
	nLines <- floor(sqrt(nParams))
	nCols <- ceiling(nParams/nLines)
	nSize <- 6
	fraction <- nCols/nLines
	dev.new(width=nSize,height=nSize/fraction)
	split.screen(c(nLines, nCols))
	for(i in 1:nParams){
		screen(i)
		plot(density(elite[,i]), xlab=paste("Parameter",i), ylab="Density", main="")

		abline(v=elite[1,i], col="red")
		abline(v=mean(elite[,i]), col="blue")
		posF <- 2
		posA <- 4
		if(elite[i,1] > mean(elite[,i])) {
			posF <- 4
			posA <- 2
		}
		text(x=elite[1,i], y=0, col="red", labels="Fittest", pos=posF)
		text(x=mean(elite[,i]), y=0, col="blue", label="Average", pos=posA)
	}	
	close.screen(all.screens=TRUE)
}
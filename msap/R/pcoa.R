


pcoa <- function(DM, groups, inds, name,surname){
	ntt <- length(levels(groups))
	#PCoA
  pcol <- cmdscale(DM, eig=T)

	var <- pcol$eig / sum(pcol$eig) *100
	var
	var1 <- round(var[1], digits=1)
	var2 <- round(var[2], digits=1)

	pcol$points<-as.data.frame(pcol$points)
	spcoo<-split(pcol$points, groups)
  maxX <- max(pcol$points$V1)
	minX <- min(pcol$points$V1)
	maxY <- max(pcol$points$V2)
	minY <- min(pcol$points$V2)
	fullname <- paste(name,surname, sep='-')
	cat("- Saving PCoA figure for ", surname, " in file: ", paste(fullname,"png",sep='.'),".....")
	png(filename=paste(fullname,"png",sep='.'))
	par(bty = 'n')
	plot(0,0, main=paste(name,surname, sep=": "), type = "n",xlab=paste("C1 (",var1,"%)"),ylab=paste("C2 (",var2,"%)"), xlim=c(minX-10^floor(log10(abs(minX))),maxX+10^floor(log10(abs(maxX)))), ylim=c(minY-10^floor(log10(abs(minY))),maxY+10^floor(log10(abs(maxY)))), frame=TRUE, cex=1.5)
	bgcolors<-rainbow(ntt+1)[-2] #No yellow?
	symbs <- c(21,22,23,24,25) #What if we have > 7 groups???
	for(i in 1:ntt){
		points(spcoo[[i]], pch=21, col="black", bg=bgcolors[i])
	}	
	s.class(pcol$points, groups, cpoint=0, col=bgcolors, add.plot=TRUE)
	dev.off()
  cat("Ok!\n")
	
	prefix <- paste(name,surname,sep='-')
	write.csv(data.frame(Group=groups,ind=inds,pcol$points), file=paste(prefix,"PCoA.coor.csv", sep='-'))
	write.csv(data.frame(pcol$eig), file=paste(prefix,"PCoA.eige.csv", sep='-'))

}

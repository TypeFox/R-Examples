###################
#   plot.recons   #
###################

plotRECON <- function(phy, likelihoods, piecolors=NULL, cex=0.5, pie.cex=0.25, file=NULL, height=11, width=8.5, show.tip.label=TRUE, title=NULL, ...){
#plotRECON <- function(phy, likelihoods, piecolors=NULL, cex=0.5, file=NULL, height=11, width=8.5, show.tip.label=TRUE, title=NULL, ...){
	if(is.null(piecolors)){
		piecolors=c("white","black","red","yellow","forestgreen","blue","coral","aquamarine","darkorchid","gold","grey","yellow","#3288BD","#E31A1C")
	}
	if(!is.null(file)){
		pdf(file, height=height, width=width,useDingbats=FALSE)
	}
	plot(phy, cex=cex, show.tip.label=show.tip.label, ...)

	if(!is.null(title)){
		title(main=title)
	}
#	nodelabels(pie=likelihoods,piecol=piecolors, cex=.25)
	nodelabels(pie=likelihoods,piecol=piecolors, cex=pie.cex)
	states <- colnames(likelihoods)
	legend(x="topleft", states, cex=0.8, pt.bg=piecolors,col="black",pch=21);

	if(!is.null(file)){
		dev.off()
	}
}

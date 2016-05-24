

### demotour ###


#' @export
#' @name demotour
#' @aliases LSD.demotour
#' @title LSD teaser
#' @description A compilation of selected plot examples.
#' @author Bjoern Schwalb
#' @seealso \code{\link{heatscatter}}, \code{\link{clusterplot}}, \code{\link{disco}}, \code{\link{colorpalette}}, \code{\link{daltonize}}
#' @examples demotour()
#' @keywords heatscatter, clusterplot, disco, colorpalette, daltonize


demotour = function()
{
	# LSD.show #
	
	LSD.show = function(it){
		emptyplot()
		text(-0.7,0.15,it,pos=4,cex=pmax(1,6-log(nchar(it))))
		points(-0.7,0,pch=3,cex=8,col="red")
	}
	
	# heatscatter #
	
	par(mfrow=c(1,1))
	LSD.show("heatscatter")
	devAskNewPage(ask = TRUE)
	points = 10^4
	x = c(rnorm(points/2),rnorm(points/2)+4)
	y = x + rnorm(points,sd=0.8)
	x = sign(x)*abs(x)^1.3
	par(mfrow=c(2,2))
	heatscatter(x,y)
	heatscatter(x,y,colpal="bl2gr2rd",main="bl2gr2rd")
	heatscatter(x,y,main="greyscales with add.contour=TRUE",add.contour=TRUE,color.contour="red",greyscale=TRUE)
	heatscatter(x,y,colpal="spectral",main="spectral with add.contour=TRUE",add.contour=TRUE)
	devAskNewPage(ask = TRUE)
	
	# clusterplot #
	
	par(mfrow=c(1,1))
	LSD.show("clusterplot")
	devAskNewPage(ask = TRUE)
	samples = 150
	probes = 75
	at = 1:probes
	clus = matrix(rnorm(probes*samples,sd=1),ncol=probes)
	clus = rbind(t(t(clus)+sin(1:probes/10))+1:nrow(clus)/samples,t(t(clus)+sin(pi/2+1:probes/10))+1:nrow(clus)/samples)
	labs = paste("cluster",kmeans(clus,4)$cluster)
	par(mfrow=c(2,2))
	clusterplot(clus,label=labs,main="Data",colpal=c("standardheat","crazyblue","crazyred","standardtopo"))
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	clusterplot(clus,label=paste("cluster",kmeans(clus,2)$cluster),separate=FALSE,main="Alpha overlay",fromto=c(0.3,0.7),colpal=c("greens","purples"),outer.col="none",ylim=c(-1,2),alpha=50,quartiles.col = c("transparent","black","transparent"))
	devAskNewPage(ask = TRUE)
	
	# daltonize #
	
	par(mfrow=c(1,1))
	LSD.show("daltonize")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	daltonize("rdylgn",cvd = "d")
	devAskNewPage(ask = TRUE)
	
	
	### disco ###
	
	nrcol = 11
	par(mfrow=c(1,1))
	LSD.show("disco")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	ownpals = c("colorblind","standard","crazyred","crazygreen","crazyblue","girly","jamaica","mountain","heat")
	pals = ownpals
	npal = length(ownpals)
	plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
	for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
	for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
	text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
	mtext(paste("palettes from the LSD package"),3,2,cex=1.25)
	mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	brewerpals = c("spectral","rdylgn","rdylbu","rdgy","rdbu","puor","prgn","piyg","brbg")
	pals = brewerpals
	npal = length(brewerpals)
	plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
	for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
	for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
	text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
	mtext(paste("palettes from the RColorBrewer package"),3,2,cex=1.25)
	mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	brewerpals = c("ylorrd","ylorbr","ylgnbu","ylgn","reds","rdpu","purples","purd","pubugn","pubu","orrd","oranges","greys","greens","gnbu","bupu","bugn","blues")
	pals = brewerpals
	npal = length(brewerpals)
	plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
	for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
	for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
	text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
	mtext(paste("palettes from the RColorBrewer package"),3,2,cex=1.25)
	mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	rampspals = c("bl2gr","bl2gr2rd","bl2rd","bl2yl","cy2yl","gr2rd","ma2gr","matlablike","matlablike2")
	pals = rampspals
	npal = length(rampspals)
	plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
	for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
	for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
	text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
	mtext(paste("palettes from the colorRamps package"),3,2,cex=1.25)
	mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	grpals = c("standardterrain","standardtopo","standardheat","standardrainbow","standardcm")
	pals = grpals
	npal = length(grpals)
	plot(1,1,xlim=c(0,nrcol),ylim=c(0,npal),type="n",axes=FALSE,bty="n",xlab="",ylab="",main="")
	for (i in 1:npal){rect(xleft = 0:(nrcol - 1),ybottom = i - 1,xright = 1:nrcol,ytop = i - 0.2,col = colorpalette(pals[i],nrcol),border = NA)}
	for (i in 1:npal){rect(xleft = 0,ybottom = i - 1,xright = nrcol,ytop = i - 0.2,col = "transparent",border = "darkgrey")}
	text(rep(-0.1,npal),(1:npal)-0.6,labels = pals,xpd = TRUE,adj = 1,cex=0.8)
	mtext(paste("palettes from the grDevices package"),3,2,cex=1.25)
	mtext(paste("( character strings can be passed to LSD functions as 'colpal' )"),3,0,col="darkgrey")
	
	# heatboxplot #
	
	par(mfrow=c(1,1))
	LSD.show("heatboxplot")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	f = c(rnorm(200),rnorm(200)+4)
	h = rf(500,15,15)*10
	g = rnorm(300)+1
	heatboxplot(list(f=f,g=g),colpals=c("rdpu","greens"),labels=c("bimodal","unimodal"))
	devAskNewPage(ask = TRUE)
	
	
	# intersphere #
	
	par(mfrow=c(1,1))
	LSD.show("intersphere")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	data = list("A" = sample(1:200,100),"B" = sample(1:200,150),"C" = sample(1:200,50),"D" = sample(1:200,75))
	intersphere(data,colors = c("orange","skyblue","green","purple"),expand.circles = 0.5,expand.lims = 0.5)
	devAskNewPage(ask = TRUE)
		
	# align #
	
	homer = list()
	homer[[1]] = c(0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
	homer[[2]] = c(0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
	homer[[3]] = c(0,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0)
	homer[[4]] = c(1,0,1,1,1,2,2,2,2,2,2,1,1,0,0,0,0,0,0,0,0,0)
	homer[[5]] = c(1,0,1,2,2,2,2,2,2,2,2,2,2,1,0,0,0,0,0,0,0,0)
	homer[[6]] = c(0,1,2,1,2,2,2,2,2,2,2,2,2,2,1,0,0,0,0,0,0,0)
	homer[[7]] = c(0,1,2,2,2,2,2,2,2,2,2,2,2,2,1,0,0,0,0,0,0,0)
	homer[[8]] = c(1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,0,0,0,0,0,0)
	homer[[9]] = c(1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,0,0,0,0)
	homer[[10]] = c(1,2,2,2,2,2,2,2,2,2,1,1,1,1,1,0,0,0,1,0,0,0)
	homer[[11]] = c(1,2,2,2,2,2,2,2,2,1,0,0,0,0,1,0,0,0,0,1,0,0)
	homer[[12]] = c(1,2,2,2,2,2,2,2,1,0,0,0,0,0,0,1,0,1,0,1,0,0)
	homer[[13]] = c(1,2,2,2,2,2,2,2,1,0,0,0,0,0,0,1,0,0,0,1,0,0)
	homer[[14]] = c(0,1,2,2,1,2,2,2,1,0,0,0,1,0,0,1,1,1,1,0,0,0)
	homer[[15]] = c(0,1,2,2,1,1,2,2,2,1,0,0,0,0,1,2,2,2,2,1,0,0)
	homer[[16]] = c(0,1,1,1,2,2,1,2,2,2,1,1,1,1,2,2,2,2,2,1,0,0)
	homer[[17]] = c(0,1,0,1,1,2,2,2,2,2,2,2,2,2,2,1,1,1,1,0,0,0)
	homer[[18]] = c(0,0,0,1,1,1,2,2,2,2,2,2,1,1,1,3,3,3,3,1,0,0)
	homer[[19]] = c(0,0,0,1,2,2,2,2,2,2,2,1,3,3,3,3,3,3,3,3,1,0)
	homer[[20]] = c(0,0,0,1,2,2,1,2,2,2,1,3,3,3,3,3,3,3,3,3,1,0)
	homer[[21]] = c(0,0,0,0,1,1,2,2,2,1,3,3,3,3,3,3,3,3,3,3,3,1)
	homer[[22]] = c(0,0,0,0,0,1,2,2,2,1,3,1,3,3,3,3,3,3,3,3,3,1)
	homer[[23]] = c(0,0,0,0,0,1,2,2,2,1,3,1,1,1,1,1,1,1,1,1,1,0)
	homer[[24]] = c(0,0,0,0,0,1,2,2,2,1,3,3,3,3,3,3,3,1,0,0,0,0)
	homer[[25]] = c(0,0,0,0,0,1,2,2,2,2,1,3,3,3,3,3,1,0,0,0,0,0)
	homer[[26]] = c(0,0,0,0,0,1,2,2,2,2,1,3,3,3,3,3,1,0,0,0,0,0)
	homer[[27]] = c(0,0,0,0,1,1,2,2,2,2,2,1,3,3,3,1,0,0,0,0,0,0)
	homer[[28]] = c(0,0,0,1,0,0,1,1,2,2,2,2,1,1,1,1,1,0,0,0,0,0)
	homer[[29]] = c(0,0,0,1,0,0,0,0,1,1,2,2,2,2,2,1,0,1,0,0,0,0)
	homer[[30]] = c(0,0,1,1,0,0,0,0,0,0,1,1,2,2,2,1,0,0,1,1,0,0)
	homer[[31]] = c(0,1,0,0,1,0,0,0,0,0,0,1,1,2,2,1,0,0,0,1,0,0)
	
	par(mfrow=c(1,1))
	LSD.show("align")
	devAskNewPage(ask = TRUE)
	par(mfrow=c(1,1))
	align(homer,colpal = c("white","black","yellow","wheat3"),main = "D'OH!",asp = 1,axes = FALSE)
	devAskNewPage(ask = TRUE)
	
	# and many more ... #

	par(mfrow=c(1,1))
	LSD.show("and many more ...")
}


### aliases ###


LSD.demotour = demotour




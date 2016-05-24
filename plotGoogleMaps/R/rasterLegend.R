rasterLegend <-
function (attribute,colPalette,legendName="Legend",bgc='#B0C4DEFF') {
png(filename=paste(legendName,'.png',sep=""), width=250, height=nlevels(factor(attribute))*25,units = "px", bg="white")
par(mai=c(0,0,0,0),bg=bgc)
plot(0,xlab="",ylab="",type="n", axes=F,xlim=c(0,3),ylim=c(0,nlevels(factor(attribute))))
	niv  <- levels(factor(attribute))
	  if(is.null(colPalette)){
	  cols<-rainbow(length(niv))
    }else{cols <-colPalette}
	k=1
	for(i in nlevels(factor(attribute)):1) {
	polygon(c(1  ,1   ,0  ,0, 1),
		        c(i-1  ,i    , i , i-1 , i-1),
			col=cols[k], border=NA)
	 text(2,i-.5,niv[k],cex=1)
    k=k+1}
     graph1 <- dev.cur()
     dev.off(graph1)
    }

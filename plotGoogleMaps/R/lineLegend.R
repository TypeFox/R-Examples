lineLegend <-
function (attribute,colPalette,linew=1,legendName="Legend",bgc='white') {
                       
		niv  <- attribute 
if(length(niv)>30){h=720}else{h=360}
png(filename =paste(legendName,'.png',sep=""),  height=h, width=144,units = "px", bg="white")
par(mai=c(0,0,0,0),bg=bgc)
plot(0,xlab="",ylab="",type="n", axes=F,xlim=c(0,3),ylim=c(0,length(niv)))

  if (length(niv)!=length(linew)) {
                           linew=rep(linew[1],length(niv)) }
cols <-colPalette

k=1
	for(i in 1: nlevels(factor(attribute))) {
	segments(0,  i-.5, 1.7,i-.5,
			col=cols[i], lwd=linew[i])
	 text(2,i-.5,niv[i],cex=1)
    k=k+1}
    graph1 <- dev.cur()
     dev.off(graph1)
    }

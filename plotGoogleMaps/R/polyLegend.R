polyLegend <-
function (attribute,colPalette=NULL,legendName="Legend",bgc='#B0C4DEFF',strokeColor="#FFAA00") {

if(strokeColor==""){strokeColor='black'}


                              
niv  <- attribute
		
png(filename =paste(legendName,'.png',sep=""), width=250, height=length(niv)*25,units = "px", bg="white")
par(mai=c(0,0,0,0),bg=bgc)
plot(0,xlab="",ylab="",type="n", axes=F,xlim=c(0,3),ylim=c(0,length(niv)))
		
cols <-colPalette

	k=1
	for(i in length(niv):1) {
	polygon(c(1  ,1   ,0  ,0, 1),
		        c(i-1  ,i    , i , i-1 , i-1),
			col=cols[k],
      border=strokeColor)
	 text(2,i-.5,niv[k],cex=1)
    k=k+1}
     graph1 <- dev.cur()
     dev.off(graph1)
    }

mapCircles <-
function (shpFile,shpId,df,dfId,var,
         fixedNorm=FALSE, shareOfCircles=0.02,radiusMax=0.5,valueMax=max(df[,var],na.rm=TRUE),
         lgdRnd=0,posLeg="bottomleft",
         circleCol="#FD8D3C",baseCol="#FFEDA0",
         title=var,legend=var,author="author",sources="sources",
         scalebar=FALSE,scalebarSize,scalebarText,
         northArrow=FALSE,northArrowSize,
         width=NULL,height=NULL,txtCex=NULL){

# required libraries
library(RColorBrewer)
library(maptools) 
  
# import of the shapefile
fdc<-readShapeSpatial(shpFile)

# data frame with x y coordinates of the centroids of the polygons
pt<-cbind(fdc@data[,shpId],as.data.frame(coordinates(fdc)))

# renaming of the columns' names
colnames(pt)<-c("Code","x","y")

# joint between the shapefile data and the data frame
pt<-merge(pt,df, by.x="Code",by.y=dfId, all.x=TRUE)

# suppression of zero values
pt<-pt[pt[,var]>0,]

# frames management 
opar<-par(mar=c(0,0,0,0))

# frames auto-size
frameSizes <- frameAutoSize(width=width,height=height,fdc=fdc)
width<-frameSizes[1]
height<-frameSizes[2]

# text auto-size
txtCex <- txtAutoSize(txtCex=txtCex,height=height)

# size of the frames
layout(matrix(c(1,2,3), 3, 1, byrow = TRUE),
       widths=c(rep(lcm(width),3)), 
       heights=c(lcm(height/10),lcm(17*height/20),lcm(height/20)))

# title display
plot.new()
par(usr=c(0,1,0,1))
text(0.5,0.5,labels=title,cex=txtCex,col="black")
segments(x0=c(0,0,1),y0=c(0.5,1,1),x1=c(0,1,1),y1=c(1,1,0.5))

# base map display
plot(fdc, border="Grey", col=baseCol)

# Size of the circles
pt<-circlesSize(fdc=fdc,fixedNorm=fixedNorm,pt=pt,var=var,
                shareOfCircles=shareOfCircles,
                radiusMax=radiusMax,valueMax=valueMax)

# reorder of the circles display (the largest first)
pt<-pt[order(pt$varSize,decreasing=TRUE),]

# cicles display
symbols(pt[,c("x","y")],circles=pt$varSize,add=TRUE,bg=circleCol,inches=FALSE,lwd=0.5)

# legend display
lgdDisplayCircles(posLegCircles=posLeg, pt=pt,varSize="varSize",var=var,
                  txtCexLeg=txtCex*0.5, legendCircles = legend , 
                  circleCol=circleCol,lgdRnd = lgdRnd)

# scalebar display
scalebarDisplay(scalebar=scalebar,
                scalebarSize=scalebarSize,scalebarText=scalebarText,
                txtCexScalebar=txtCex*0.5)

# north arrow display 
northArrowDisplay(northArrow=northArrow,northArrowSize=northArrowSize)

# author an sources display
plot.new()
par(usr=c(0,1,0,1))
text(0.05,0.7,labels=author,cex=txtCex*0.4,adj=0,col="black")
text(0.05,0.3,labels=sources,cex=txtCex*0.4,adj=0,col="black")
segments(x0=c(1,1,0),y0=c(1,0,0),x1=c(1,0,0),y1=c(0,0,1))

# re-init to the originals graphicals parameters
par(opar)
}

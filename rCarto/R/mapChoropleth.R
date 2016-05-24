mapChoropleth <-
function (shpFile,shpId,df,dfId,var,
          nclass=6,style="quantile",fixBrks=FALSE,listBrks=NULL,diverg=FALSE,divergBrk=0,
          lgdRnd=2,posLeg="bottomleft",
          palCol="Greens",palColPos="Reds",palColNeg="Blues",NACol="grey",
          title=var,legend=var,author="author",sources="sources",
          scalebar=FALSE,scalebarSize,scalebarText,
          northArrow=FALSE,northArrowSize,
					width=NULL,height=NULL,txtCex=NULL){

# required libraries
library(RColorBrewer)
library(maptools)
library(classInt)

# import of the shapefile
fdc<-readShapeSpatial(shpFile)

# joint between the shapefile data and the data frame
fdc@data <- data.frame(fdc@data, df[match(fdc@data[,shpId], df[,dfId]),])

# discretization of the variable
discretParam<-discretVar(fixBrks=fixBrks,listBrks=listBrks,pt=fdc@data,var=var,nclass=nclass,
                         style=style,palCol=palCol,diverg=diverg,divergBrk=divergBrk,palColPos=palColPos,
                         palColNeg=palColNeg,NACol=NACol,lgdRnd=lgdRnd)

fdc@data<-discretParam[[1]]
colours<-discretParam[[2]]
pdd<-discretParam[[3]]
lblLeg<-discretParam[[4]]

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

# map display
plot(fdc, col=fdc$col)

# Distr legend display
lgdDisplayDistr (posLeg=posLeg, 
                 lblLeg=lblLeg, colours=colours,  
                 na.leg=pdd , txtCexLeg=txtCex*0.5, 
                 legend=legend,NACol=NACol)

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

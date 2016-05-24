### R code from vignette source 'Vquantile.Rnw'

###################################################
### code chunk number 1: ExSetup set options
###################################################
options(prompt = "R> ")
options(width = 75)
options(continue=" ")
pdf("Vplots.pdf")


###################################################
### code chunk number 2: ExSetup set width
###################################################
options(width=60)


###################################################
### code chunk number 3: ExSetup load package
###################################################
library("ModelMap")


###################################################
### code chunk number 4: ExElevReadGDAL
###################################################

library(raster)
elevfn <- paste(getwd(),"/VModelMapData_dem_ELEVM_250.img",sep="")
mapgrid <- raster(elevfn)


###################################################
### code chunk number 5: ExElev
###################################################

opar <- par(mar=c(4,4,3,6),xpd=NA,mgp=c(3, 2, .3))

col.ramp<-terrain.colors(101)

zlim <- c(1500,maxValue(mapgrid))
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-col.ramp[trunc((legend.label/max(legend.label))*100)+1]
legend.label<-paste(legend.label,"m",sep="")

legend.label<-paste((7:3)*500,"m")
legend.colors<-col.ramp[c(100,75,50,25,1)]

image( mapgrid, 
       col = col.ramp,
       xlab="", ylab="", 
       zlim=zlim,
       asp=1, bty="n", main="")

legend( x=xmax(mapgrid),y=ymax(mapgrid),
        legend=legend.label,
        fill=legend.colors,
        bty="n",
        cex=1.2)
mtext("Elevation of Study Region",side=3,line=1,cex=1.5)
par(opar)


###################################################
### code chunk number 6: ExElevFig
###################################################

opar <- par(mar=c(4,4,3,6),xpd=NA,mgp=c(3, 2, .3))

col.ramp<-terrain.colors(101)

zlim <- c(1500,maxValue(mapgrid))
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-col.ramp[trunc((legend.label/max(legend.label))*100)+1]
legend.label<-paste(legend.label,"m",sep="")

legend.label<-paste((7:3)*500,"m")
legend.colors<-col.ramp[c(100,75,50,25,1)]

image( mapgrid, 
       col = col.ramp,
       xlab="", ylab="", 
       zlim=zlim,
       asp=1, bty="n", main="")

legend( x=xmax(mapgrid),y=ymax(mapgrid),
        legend=legend.label,
        fill=legend.colors,
        bty="n",
        cex=1.2)
mtext("Elevation of Study Region",side=3,line=1,cex=1.5)
par(opar)


###################################################
### code chunk number 7: ExSetup Define training and test files
###################################################
qdatafn <- "VModelMapData.csv"
qdata.trainfn <- "VModelMapData_TRAIN.csv"
qdata.testfn  <- "VModelMapData_TEST.csv"


###################################################
### code chunk number 8: ExSetup define folder
###################################################
folder <- getwd()


###################################################
### code chunk number 9: ExSetup split training and test
###################################################
get.test(       proportion.test=0.2,
                qdatafn=qdatafn,
                seed=42,
                folder=folder,
                qdata.trainfn=qdata.trainfn,
                qdata.testfn=qdata.testfn)


###################################################
### code chunk number 10: ExSetup Define predictors
###################################################
predList <- c( "ELEV250",
               "NLCD01_250",
               "EVI2005097",
               "NDV2005097",
               "NIR2005097",
               "RED2005097")
predFactor <- c("NLCD01_250")


###################################################
### code chunk number 11: Ex1 Define Identifier
###################################################
unique.rowname <- "ID"


###################################################
### code chunk number 12: ExSetup update raster LUT
###################################################
rastLUTfn     <- "VModelMapData_LUT.csv"
rastLUTfn     <- read.table( rastLUTfn,
                             header=FALSE,
                             sep=",",
                             stringsAsFactors=FALSE)
rastLUTfn[,1] <- paste(folder,rastLUTfn[,1],sep="/")


###################################################
### code chunk number 13: Ex1 Define model filenames
###################################################
MODELfn.pinyon    <- "VQuantile_QRF_Pinyon"
MODELfn.sage    <- "VQuantile_QRF_Sage"


###################################################
### code chunk number 14: Ex1 Define response
###################################################
response.name.pinyon <- "PINYON"
response.name.sage <- "SAGE"
response.type   <- "continuous"


###################################################
### code chunk number 15: Ex1 Create Model
###################################################
QRF.pinyon <- model.build( model.type="QRF",
                               qdata.trainfn=qdata.trainfn,
                               folder=folder,
                               unique.rowname=unique.rowname,          
                               MODELfn=MODELfn.pinyon,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.pinyon,
                               response.type=response.type,
                               importance=TRUE, 
                               quantiles=c(0.1,0.2,0.5,0.8,0.9))
           
QRF.sage <- model.build( model.type="QRF",
                               qdata.trainfn=qdata.trainfn,
                               folder=folder,
                               unique.rowname=unique.rowname,          
                               MODELfn=MODELfn.sage,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.sage,
                               response.type=response.type,
                               importance=TRUE, 
                               quantiles=c(0.1,0.2,0.5,0.8,0.9))


###################################################
### code chunk number 16: Ex1 Model Diagnostics
###################################################
QRF.pinyon.pred <- model.diagnostics( model.obj=QRF.pinyon,
                                      qdata.testfn=qdata.testfn,
                                      folder=folder,           
                                      MODELfn=MODELfn.pinyon,
                                      unique.rowname=unique.rowname,
                                      quantiles=c(0.1,0.5,0.9),
                             # Model Validation Arguments
                                      prediction.type="TEST",
                                      device.type=c("pdf","png"),
                                      cex=1.2)
           
QRF.sage.pred <- model.diagnostics( model.obj=QRF.sage,
                                      qdata.testfn=qdata.testfn,
                                      folder=folder,           
                                      MODELfn=MODELfn.sage,
                                      unique.rowname=unique.rowname,
                                      quantiles=c(0.1,0.5,0.9),
                              # Model Validation Arguments
                                      prediction.type="TEST",
                                      device.type=c("pdf","png"),
                                      cex=1.2)


###################################################
### code chunk number 17: Ex1CompImpRF_QRF
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

Imp.pinyon<-model.importance.plot( model.obj.1=QRF.pinyon$RF,
                       model.obj.2=QRF.pinyon$QRF,
                       model.name.1="RFmean",
                       model.name.2="QRF median",
                       quantile.2=0.5,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Percent Cover",
                       device.type="none",
                       cex=0.9)
                      
Imp.sage<-model.importance.plot( model.obj.1=QRF.sage$RF,
                       model.obj.2=QRF.sage$QRF,
                       model.name.1="RFmean",
                       model.name.2="QRF median",
                       quantile.2=0.5,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Percent Cover",
                       device.type="none",
                       cex=0.9)
                       
mtext("Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 18: Ex1CompImpFigRF_QRF
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

Imp.pinyon<-model.importance.plot( model.obj.1=QRF.pinyon$RF,
                       model.obj.2=QRF.pinyon$QRF,
                       model.name.1="RFmean",
                       model.name.2="QRF median",
                       quantile.2=0.5,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Percent Cover",
                       device.type="none",
                       cex=0.9)
                      
Imp.sage<-model.importance.plot( model.obj.1=QRF.sage$RF,
                       model.obj.2=QRF.sage$QRF,
                       model.name.1="RFmean",
                       model.name.2="QRF median",
                       quantile.2=0.5,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Percent Cover",
                       device.type="none",
                       cex=0.9)
                       
mtext("Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 19: Ex1CompImpLow_Up
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

Imp.pinyon<-model.importance.plot( model.obj.1=QRF.pinyon$QRF,
                       model.obj.2=QRF.pinyon$QRF,
                       model.name.1="Lower Quantile",
                       model.name.2="Upper Quantile",
                       quantile.1=0.1,
                       quantile.2=0.9,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Percent Cover",
                       device.type="none",
                       cex=0.9)
                      
Imp.sage<-model.importance.plot( model.obj.1=QRF.sage$QRF,
                       model.obj.2=QRF.sage$QRF,
                       model.name.1="Lower Quantile",
                       model.name.2="Upper Quantile",
                       quantile.1=0.2,
                       quantile.2=0.8,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Percent Cover",
                       device.type="none",
                       cex=0.9)
                       
mtext("Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 20: Ex1CompImpFigLow_Up
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

Imp.pinyon<-model.importance.plot( model.obj.1=QRF.pinyon$QRF,
                       model.obj.2=QRF.pinyon$QRF,
                       model.name.1="Lower Quantile",
                       model.name.2="Upper Quantile",
                       quantile.1=0.1,
                       quantile.2=0.9,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Percent Cover",
                       device.type="none",
                       cex=0.9)
                      
Imp.sage<-model.importance.plot( model.obj.1=QRF.sage$QRF,
                       model.obj.2=QRF.sage$QRF,
                       model.name.1="Lower Quantile",
                       model.name.2="Upper Quantile",
                       quantile.1=0.2,
                       quantile.2=0.8,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Percent Cover",
                       device.type="none",
                       cex=0.9)
                       
mtext("Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 21: Ex1InteractionPlots
###################################################

pred.means <- c( ELEV250    = 2300,
                 NLCD01_250 = 42,
                 EVI2005097 = 1800,
                 NDV2005097 = 3100,
                 NIR2005097 = 2000,
                 RED2005097 = 1000)


opar <- par(mfrow=c(3,1),mar=c(2,3,0,2),oma=c(0,0,3,0))
for(quantile in c(0.1,0.5,0.9)){
     model.interaction.plot( QRF.pinyon$QRF,
                             x="ELEV250",
                             y="EVI2005097",
                             main="",
                             plot.type="persp",
                             device.type=c("none"),
                             MODELfn=MODELfn.pinyon,
                             folder=paste(folder,"/interaction",sep=""),
                             quantile=quantile,
                             pred.means=pred.means,
                             zlim=c(0,60))
     mtext(paste( quantile*100,"% Quantile",sep=""),side=2,line=1,font=2,cex=1)
}

mtext("Pinyon Cover by Quantile",side=3,line=1,cex=1.8,outer=TRUE)
par(opar)



###################################################
### code chunk number 22: Ex1InteractionPlots
###################################################

pred.means <- c( ELEV250    = 2300,
                 NLCD01_250 = 42,
                 EVI2005097 = 1800,
                 NDV2005097 = 3100,
                 NIR2005097 = 2000,
                 RED2005097 = 1000)


opar <- par(mfrow=c(3,1),mar=c(2,3,0,2),oma=c(0,0,3,0))
for(quantile in c(0.1,0.5,0.9)){
     model.interaction.plot( QRF.pinyon$QRF,
                             x="ELEV250",
                             y="EVI2005097",
                             main="",
                             plot.type="persp",
                             device.type=c("none"),
                             MODELfn=MODELfn.pinyon,
                             folder=paste(folder,"/interaction",sep=""),
                             quantile=quantile,
                             pred.means=pred.means,
                             zlim=c(0,60))
     mtext(paste( quantile*100,"% Quantile",sep=""),side=2,line=1,font=2,cex=1)
}

mtext("Pinyon Cover by Quantile",side=3,line=1,cex=1.8,outer=TRUE)
par(opar)



###################################################
### code chunk number 23: Ex1 update raster LUT
###################################################
rastLUTfn     <- "VModelMapData_LUT.csv"
rastLUTfn     <- read.table(  rastLUTfn,
                              header=FALSE,
                              sep=",",
                              stringsAsFactors=FALSE)
rastLUTfn[,1] <- paste(folder,rastLUTfn[,1],sep="/")


###################################################
### code chunk number 24: Ex1 Produce Maps
###################################################

model.mapmake(  model.obj=QRF.pinyon,
                folder=folder,           
                MODELfn=MODELfn.pinyon,
                rastLUTfn=rastLUTfn,
                na.action="na.omit",
             # Mapping arguments
                map.sd=TRUE,
                quantiles=c(0.1,0.5,0.9))
            

model.mapmake(  model.obj=QRF.sage,
                folder=folder,           
                MODELfn=MODELfn.sage,
                rastLUTfn=rastLUTfn,
                na.action="na.omit",
             # Mapping arguments
                map.sd=TRUE,
                quantiles=c(0.1,0.5,0.9))


###################################################
### code chunk number 25: Ex1 define color sequence
###################################################
l <- seq(100,0,length.out=101)
c <- seq(0,100,length.out=101)
col.ramp <- hcl(h = 120, c = c, l = l)


###################################################
### code chunk number 26: Ex1RFMap
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.pinyon <- raster(paste(MODELfn.pinyon,"_map_RF.img",sep=""))
mapgrid.sage <- raster(paste(MODELfn.sage,"_map_RF.img",sep=""))

zlim <- c(0,60)
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-col.ramp[trunc((legend.label/max(legend.label))*100)+1]
legend.label<-paste(legend.label,"%",sep="")

image( mapgrid.pinyon,
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage, 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("RF - Mean Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 27: Ex1RFMapFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.pinyon <- raster(paste(MODELfn.pinyon,"_map_RF.img",sep=""))
mapgrid.sage <- raster(paste(MODELfn.sage,"_map_RF.img",sep=""))

zlim <- c(0,60)
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-col.ramp[trunc((legend.label/max(legend.label))*100)+1]
legend.label<-paste(legend.label,"%",sep="")

image( mapgrid.pinyon,
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage, 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("RF - Mean Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 28: Ex1QRFMap
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.pinyon <- brick(paste(MODELfn.pinyon,"_map_QRF.img",sep=""))
mapgrid.sage <- brick(paste(MODELfn.sage,"_map_QRF.img",sep=""))

image( mapgrid.pinyon[[2]],
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage[[2]], 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("QRF - Median (50% quantile) Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 29: Ex1QRFMapFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.pinyon <- brick(paste(MODELfn.pinyon,"_map_QRF.img",sep=""))
mapgrid.sage <- brick(paste(MODELfn.sage,"_map_QRF.img",sep=""))

image( mapgrid.pinyon[[2]],
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage[[2]], 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("QRF - Median (50% quantile) Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 30: Ex1QRFMapLow
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)


image( mapgrid.pinyon[[1]],
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage[[1]], 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("QRF - Lower bound (10% quantile)",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 31: Ex1QRFMapLowFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)


image( mapgrid.pinyon[[1]],
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage[[1]], 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("QRF - Lower bound (10% quantile)",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 32: Ex1QRFMapUp
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

image( mapgrid.pinyon[[3]],
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage[[3]], 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("QRF model - Upper - 90% quantile",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 33: Ex1QRFMapUpFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

image( mapgrid.pinyon[[3]],
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage[[3]], 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("QRF model - Upper - 90% quantile",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 34: Ex1 Define model filenames
###################################################
MODELfn.pinyon    <- "VQuantile_CF_Pinyon"
MODELfn.sage    <- "VQuantile_CF_Sage"


###################################################
### code chunk number 35: Ex1 Define response
###################################################
response.name.pinyon <- "PINYON"
response.name.sage <- "SAGE"
response.type   <- "continuous"


###################################################
### code chunk number 36: Ex1 Create Model
###################################################

  qdata<-read.csv(qdata.trainfn)
	IS.NUM<-sapply(qdata,is.numeric)
	qdata[,IS.NUM]<-sapply(qdata[,IS.NUM],as.numeric)

CF.pinyon  <- model.build( model.type="CF",
                               qdata.trainfn=qdata,
                               folder=folder,
                               unique.rowname=unique.rowname,          
                               MODELfn=MODELfn.pinyon,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.pinyon,
                               response.type=response.type)
           
CF.sage  <- model.build( model.type="CF",
                               qdata.trainfn=qdata,
                               folder=folder,
                               unique.rowname=unique.rowname,          
                               MODELfn=MODELfn.sage,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.sage,
                               response.type=response.type)


###################################################
### code chunk number 37: Ex2CompImpCF_RF
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

Imp.pinyon<-model.importance.plot( model.obj.1=CF.pinyon,
                       model.obj.2=QRF.pinyon$RF,
                       model.name.1="conditional (CF)",
                       model.name.2="unconditional (RF)",
                       cf.conditional.1=FALSE,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Percent Cover",
                       device.type="none",
                       cex=0.9)
                      
Imp.sage<-model.importance.plot( model.obj.1=CF.sage,
                       model.obj.2=QRF.sage$RF,
                       model.name.1="conditional (CF)",
                       model.name.2="unconditional (RF)",
                       cf.conditional.1=FALSE,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Percent Cover",
                       device.type="none",
                       cex=0.9)
                       
mtext("CF versus RF Variable Importance",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 38: Ex2CompImpFigCF_RF
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

Imp.pinyon<-model.importance.plot( model.obj.1=CF.pinyon,
                       model.obj.2=QRF.pinyon$RF,
                       model.name.1="conditional (CF)",
                       model.name.2="unconditional (RF)",
                       cf.conditional.1=FALSE,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Percent Cover",
                       device.type="none",
                       cex=0.9)
                      
Imp.sage<-model.importance.plot( model.obj.1=CF.sage,
                       model.obj.2=QRF.sage$RF,
                       model.name.1="conditional (CF)",
                       model.name.2="unconditional (RF)",
                       cf.conditional.1=FALSE,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Percent Cover",
                       device.type="none",
                       cex=0.9)
                       
mtext("CF versus RF Variable Importance",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 39: Ex1 Produce Maps
###################################################

model.mapmake(  model.obj=CF.pinyon,
                folder=folder,           
                MODELfn=MODELfn.pinyon,
                rastLUTfn=rastLUTfn,
                na.action="na.omit"
             )
            

model.mapmake(  model.obj=CF.sage,
                folder=folder,           
                MODELfn=MODELfn.sage,
                rastLUTfn=rastLUTfn,
                na.action="na.omit"
             )


###################################################
### code chunk number 40: Ex2CFMap
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.pinyon <- raster(paste(MODELfn.pinyon,"_map.img",sep=""))
mapgrid.sage <- raster(paste(MODELfn.sage,"_map.img",sep=""))

zlim <- c(0,60)

image( mapgrid.pinyon,
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage, 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("CF - Mean Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 41: Ex2CFMapFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.pinyon <- raster(paste(MODELfn.pinyon,"_map.img",sep=""))
mapgrid.sage <- raster(paste(MODELfn.sage,"_map.img",sep=""))

zlim <- c(0,60)

image( mapgrid.pinyon,
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.pinyon,side=3,line=1,cex=1.2)
image( mapgrid.sage, 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.sage,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.sage),y=ymax(mapgrid.sage),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("CF - Mean Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 42: Remove Rplots pdf
###################################################
dev.off()

file.remove("Vplots.pdf")



### R code from vignette source 'VModelMap.Rnw'

###################################################
### code chunk number 1: Ex1 set options
###################################################
options(prompt = "R> ")
options(width = 75)
options(continue=" ")
pdf("Vplots.pdf")


###################################################
### code chunk number 2: Ex1 set width
###################################################
options(width=60)


###################################################
### code chunk number 3: Ex1 load package
###################################################
library("ModelMap")


###################################################
### code chunk number 4: Ex1 Define model type
###################################################
model.type <- "RF"


###################################################
### code chunk number 5: Ex1 Define training and test files
###################################################
qdatafn <- "VModelMapData.csv"
qdata.trainfn <- "VModelMapData_TRAIN.csv"
qdata.testfn  <- "VModelMapData_TEST.csv"


###################################################
### code chunk number 6: Ex1 define folder
###################################################
folder <- getwd()


###################################################
### code chunk number 7: Ex1 split training and test
###################################################
get.test(       proportion.test=0.2,
                qdatafn=qdatafn,
                seed=42,
                folder=folder,
                qdata.trainfn=qdata.trainfn,
                qdata.testfn=qdata.testfn)


###################################################
### code chunk number 8: Ex1 Define model filenames
###################################################
MODELfn.a    <- "VModelMapEx1a"
MODELfn.b    <- "VModelMapEx1b"


###################################################
### code chunk number 9: Ex1 Define predictors
###################################################
predList <- c( "ELEV250",
                "EVI2005097",
                "NDV2005097",
                "NIR2005097",
                "RED2005097")
predFactor <- FALSE


###################################################
### code chunk number 10: Ex1 Define response
###################################################
response.name.a <- "PINYON"
response.name.b <- "SAGE"
response.type   <- "continuous"


###################################################
### code chunk number 11: Ex1 set seed
###################################################
seed.a <- 38
seed.b <- 39


###################################################
### code chunk number 12: Ex1 Define Identifier
###################################################
unique.rowname <- "ID"


###################################################
### code chunk number 13: Ex1 Create Model
###################################################
model.obj.ex1a <- model.build( model.type=model.type,
                               qdata.trainfn=qdata.trainfn,
                               folder=folder,
                               unique.rowname=unique.rowname,          
                               MODELfn=MODELfn.a,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.a,
                               response.type=response.type,
                               seed=seed.a)
           
model.obj.ex1b <- model.build( model.type=model.type,
                               qdata.trainfn=qdata.trainfn,
                               folder=folder,
                               unique.rowname=unique.rowname,          
                               MODELfn=MODELfn.b,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.b,
                               response.type=response.type,
                               seed=seed.b)


###################################################
### code chunk number 14: Ex1 Model Diagnostics
###################################################
model.pred.ex1a <- model.diagnostics( model.obj=model.obj.ex1a,
                                      qdata.testfn=qdata.testfn,
                                      folder=folder,           
                                      MODELfn=MODELfn.a,
                                      unique.rowname=unique.rowname,
                             # Model Validation Arguments
                                      prediction.type="TEST",
                                      device.type=c("pdf"),
                                      cex=1.2)
           
model.pred.ex1b <- model.diagnostics( model.obj=model.obj.ex1b,
                                      qdata.testfn=qdata.testfn,
                                      folder=folder,           
                                      MODELfn=MODELfn.b,
                                      unique.rowname=unique.rowname,
                              # Model Validation Arguments
                                      prediction.type="TEST",
                                      device.type=c("pdf"),
                                      cex=1.2)


###################################################
### code chunk number 15: Ex1 Compare Importance Plots
###################################################
model.importance.plot( model.obj.1=model.obj.ex1a,
                       model.obj.2=model.obj.ex1b,
                       model.name.1="Pinyon",
                       model.name.2="Sage",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Variable Importance",
                       device.type="pdf",
                       PLOTfn="VModelMapEx1CompareImportance",
                       folder=folder)
												


###################################################
### code chunk number 16: Ex1CompImpType
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex1a,
                       model.obj.2=model.obj.ex1a,
                       model.name.1="",
                       model.name.2="",
                       imp.type.1=1,
                       imp.type.2=2,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon",
                       device.type="none",
                       cex=0.9) 
                      
model.importance.plot( model.obj.1=model.obj.ex1b,
                       model.obj.2=model.obj.ex1b,
                       model.name.1="",
                       model.name.2="",
                       imp.type.1=1,
                       imp.type.2=2,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage",
                       device.type="none",
                       cex=0.9)
                       
mtext("Comparison of Importance Types",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 17: Ex1CompImpTypeFig
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex1a,
                       model.obj.2=model.obj.ex1a,
                       model.name.1="",
                       model.name.2="",
                       imp.type.1=1,
                       imp.type.2=2,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon",
                       device.type="none",
                       cex=0.9) 
                      
model.importance.plot( model.obj.1=model.obj.ex1b,
                       model.obj.2=model.obj.ex1b,
                       model.name.1="",
                       model.name.2="",
                       imp.type.1=1,
                       imp.type.2=2,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage",
                       device.type="none",
                       cex=0.9)
                       
mtext("Comparison of Importance Types",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 18: Ex1 Interaction Plots
###################################################
model.interaction.plot( model.obj.ex1a,
                        x="NIR2005097",
                        y="RED2005097",
                        main=response.name.a,
                        plot.type="image",
                        device.type="pdf",
                        MODELfn=MODELfn.a,
                        folder=folder)
												
model.interaction.plot( model.obj.ex1b,
                        x="NIR2005097",
                        y="RED2005097",
                        main=response.name.b,
                        plot.type="image",
                        device.type="pdf",
                        MODELfn=MODELfn.b,
                        folder=folder)
														
model.interaction.plot( model.obj.ex1a,
                        x=1,
                        y=3,
                        main=response.name.a,
                        plot.type="image",
                        device.type="pdf",
                        MODELfn=MODELfn.a,
                        folder=folder)
												
model.interaction.plot( model.obj.ex1b,
                        x=1,
                        y=3,
                        main=response.name.b,
                        plot.type="image",
                        device.type="pdf",
                        MODELfn=MODELfn.b,
                        folder=folder)


###################################################
### code chunk number 19: ExElevReadGDAL
###################################################

elevfn <- paste(folder,"/VModelMapData_dem_ELEVM_250.img",sep="")
mapgrid <- raster(elevfn)


###################################################
### code chunk number 20: ExElev
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
### code chunk number 21: ExElevFig
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
### code chunk number 22: Ex1 update raster LUT
###################################################
rastLUTfn     <- "VModelMapData_LUT.csv"
rastLUTfn     <- read.table(  rastLUTfn,
                              header=FALSE,
                              sep=",",
                              stringsAsFactors=FALSE)
rastLUTfn[,1] <- paste(folder,rastLUTfn[,1],sep="/")


###################################################
### code chunk number 23: Ex1 Define numrows
###################################################



###################################################
### code chunk number 24: Ex1 Produce Maps
###################################################

model.mapmake(  model.obj=model.obj.ex1a,
                folder=folder,           
                MODELfn=MODELfn.a,
                rastLUTfn=rastLUTfn,
                na.action="na.omit",
             # Mapping arguments
                map.sd=TRUE)
            

model.mapmake(  model.obj=model.obj.ex1b,
                folder=folder,           
                MODELfn=MODELfn.b,
                rastLUTfn=rastLUTfn,
                na.action="na.omit",
             # Mapping arguments
                map.sd=TRUE)


###################################################
### code chunk number 25: Ex1 define color sequence
###################################################
l <- seq(100,0,length.out=101)
c <- seq(0,100,length.out=101)
col.ramp <- hcl(h = 120, c = c, l = l)


###################################################
### code chunk number 26: Ex1Map
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <- raster(paste(MODELfn.a,"_map.img",sep=""))
mapgrid.b <- raster(paste(MODELfn.b,"_map.img",sep=""))

zlim <- c(0,max(maxValue(mapgrid.a),maxValue(mapgrid.b)))
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-col.ramp[trunc((legend.label/max(legend.label))*100)+1]
legend.label<-paste(legend.label,"%",sep="")

image( mapgrid.a,
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.a,side=3,line=1,cex=1.2)
image( mapgrid.b, 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 27: Ex1MapFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <- raster(paste(MODELfn.a,"_map.img",sep=""))
mapgrid.b <- raster(paste(MODELfn.b,"_map.img",sep=""))

zlim <- c(0,max(maxValue(mapgrid.a),maxValue(mapgrid.b)))
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-col.ramp[trunc((legend.label/max(legend.label))*100)+1]
legend.label<-paste(legend.label,"%",sep="")

image( mapgrid.a,
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.a,side=3,line=1,cex=1.2)
image( mapgrid.b, 
       col=col.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 28: Ex1 define stdev color sequence
###################################################
stdev.ramp   <- hcl(h = 15, c = c, l = l)


###################################################
### code chunk number 29: Ex1StdevMap
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <- raster(paste(MODELfn.a,"_map_stdev.img",sep=""))
mapgrid.b <- raster(paste(MODELfn.b,"_map_stdev.img",sep=""))

zlim <- c(0,max(maxValue(mapgrid.a),maxValue(mapgrid.b)))
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-stdev.ramp[trunc((legend.label/max(legend.label))*100)+1]
legend.label<-paste(legend.label,"%",sep="")

image( mapgrid.a, 
       col=stdev.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.a,side=3,line=1,cex=1.2)
image( mapgrid.b, 
       col=stdev.ramp,xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("Standard Deviation of Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 30: Ex1StdevMapFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <- raster(paste(MODELfn.a,"_map_stdev.img",sep=""))
mapgrid.b <- raster(paste(MODELfn.b,"_map_stdev.img",sep=""))

zlim <- c(0,max(maxValue(mapgrid.a),maxValue(mapgrid.b)))
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-stdev.ramp[trunc((legend.label/max(legend.label))*100)+1]
legend.label<-paste(legend.label,"%",sep="")

image( mapgrid.a, 
       col=stdev.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.a,side=3,line=1,cex=1.2)
image( mapgrid.b, 
       col=stdev.ramp,xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("Standard Deviation of Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 31: Ex1 define coefv color sequence
###################################################
coefv.ramp <- hcl(h = 70, c = c, l = l)


###################################################
### code chunk number 32: Ex1CoefvMap
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <-raster(paste(MODELfn.a,"_map_coefv.img",sep=""),as.image=TRUE)
mapgrid.b <- raster(paste(MODELfn.b,"_map_coefv.img",sep=""),as.image=TRUE)

zlim <- c(0,max(maxValue(mapgrid.a),maxValue(mapgrid.b)))
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-coefv.ramp[trunc((legend.label/max(legend.label))*100)+1]

image( mapgrid.a, 
       col=coefv.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.a,side=3,line=1,cex=1.2)
image( mapgrid.b, 
       col=coefv.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("Coefficient of Variation of Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 33: Ex1CoefvMapFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <-raster(paste(MODELfn.a,"_map_coefv.img",sep=""),as.image=TRUE)
mapgrid.b <- raster(paste(MODELfn.b,"_map_coefv.img",sep=""),as.image=TRUE)

zlim <- c(0,max(maxValue(mapgrid.a),maxValue(mapgrid.b)))
legend.label<-rev(pretty(zlim,n=5))
legend.colors<-coefv.ramp[trunc((legend.label/max(legend.label))*100)+1]

image( mapgrid.a, 
       col=coefv.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.a,side=3,line=1,cex=1.2)
image( mapgrid.b, 
       col=coefv.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=zlim,
       asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
	    	legend=legend.label,
		    fill=legend.colors,
		    bty="n",
		    cex=1.2)
mtext("Coefficient of Variation of Percent Cover",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 34: Ex2 Define model type
###################################################
model.type <- "RF"


###################################################
### code chunk number 35: Ex2 Define training and test files
###################################################
qdatafn <- "VModelMapData.csv"


###################################################
### code chunk number 36: Ex2 define folder
###################################################
folder <- getwd()


###################################################
### code chunk number 37: Ex2 Define model filename
###################################################
MODELfn.a    <- "VModelMapEx2a"
MODELfn.b    <- "VModelMapEx2b"


###################################################
### code chunk number 38: Ex2 Define predictors
###################################################
predList <- c( "ELEV250",
               "NLCD01_250",
               "EVI2005097",
               "NDV2005097",
               "NIR2005097",
               "RED2005097")
predFactor <- c("NLCD01_250")


###################################################
### code chunk number 39: Ex2 Define response
###################################################
response.name.a <- "PINYON"
response.name.b <- "SAGE"
response.type <- "binary"


###################################################
### code chunk number 40: Ex2 set seed
###################################################
seed.a <- 40
seed.b <- 41


###################################################
### code chunk number 41: Ex2 Define Identifier
###################################################
unique.rowname <- "ID"


###################################################
### code chunk number 42: Ex2 update raster LUT
###################################################
rastLUTfn     <- "VModelMapData_LUT.csv"
rastLUTfn     <- read.table( rastLUTfn,
                             header=FALSE,
                             sep=",",
                             stringsAsFactors=FALSE)
rastLUTfn[,1] <- paste(folder,rastLUTfn[,1],sep="/")


###################################################
### code chunk number 43: Ex2 Create Model
###################################################
model.obj.ex2a <- model.build( model.type=model.type,
                               qdata.trainfn=qdatafn,
                               folder=folder,
                               unique.rowname=unique.rowname,               
                               MODELfn=MODELfn.a,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.a,
                               response.type=response.type,
                               seed=seed.a)
           
model.obj.ex2b <- model.build( model.type=model.type,
                               qdata.trainfn=qdatafn,
                               folder=folder,
                               unique.rowname=unique.rowname,               
                               MODELfn=MODELfn.b,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.b,
                               response.type=response.type,
                               seed=seed.b)


###################################################
### code chunk number 44: Ex2 Model Diagnostics
###################################################
model.pred.ex2a <- model.diagnostics( model.obj=model.obj.ex2a,
                                      qdata.trainfn=qdatafn,
                                      folder=folder,           
                                      MODELfn=MODELfn.a,
                                      unique.rowname=unique.rowname,
                             # Model Validation Arguments
                                      prediction.type="OOB",
                                      device.type=c("jpeg","pdf","postscript"),
                                      cex=1.2)
           
model.pred.ex2b <- model.diagnostics( model.obj=model.obj.ex2b,
                                      qdata.trainfn=qdatafn,
                                      folder=folder,           
                                      MODELfn=MODELfn.b,
                                      unique.rowname=unique.rowname,
                              # Model Validation Arguments
                                      prediction.type="OOB",
                                      device.type=c("jpeg","pdf","postscript"),
                                      cex=1.2)


###################################################
### code chunk number 45: Ex2 Optimal Threshold Table
###################################################
opt.thresh.a <- read.table( paste(MODELfn.a,"_pred_optthresholds.csv",sep=""),
                            header=TRUE,
                            sep=",",
                            stringsAsFactors=FALSE)
opt.thresh.a[,-1]<-signif(opt.thresh.a[,-1],2)

opt.thresh.b <- read.table( paste(MODELfn.b,"_pred_optthresholds.csv",sep=""),
                            header=TRUE,
                            sep=",",
                            stringsAsFactors=FALSE)
opt.thresh.b[,-1]<-signif(opt.thresh.b[,-1],2)


pred.prev.a <- read.table( paste(MODELfn.a,"_pred_prevalence.csv",sep=""),
                           header=TRUE,
                           sep=",",
                           stringsAsFactors=FALSE)
pred.prev.a[,-1]<-signif(pred.prev.a[,-1],2)

pred.prev.b <- read.table( paste(MODELfn.b,"_pred_prevalence.csv",sep=""),
                           header=TRUE,
                           sep=",",
                           stringsAsFactors=FALSE)
pred.prev.b[,-1]<-signif(pred.prev.b[,-1],2)


###################################################
### code chunk number 46: Ex2a Optimal Threshold Table Show
###################################################
opt.thresh.a


###################################################
### code chunk number 47: Ex2b Optimal Threshold Table Show
###################################################
opt.thresh.b


###################################################
### code chunk number 48: Ex2a Prevalence Table Show
###################################################
pred.prev.a


###################################################
### code chunk number 49: Ex2b Prevalence Table Show
###################################################
pred.prev.b


###################################################
### code chunk number 50: Ex2 Compare Importance Plots
###################################################
model.importance.plot( model.obj.1=model.obj.ex2a,
                       model.obj.2=model.obj.ex2b,
                       model.name.1="Pinyon",
                       model.name.2="Sage",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Variable Importance",
                       device.type="pdf",
                       PLOTfn="VModelMapEx2CompareImportance",
                       folder=folder)
												


###################################################
### code chunk number 51: Ex2CompImpPA
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex2a,
                       model.obj.2=model.obj.ex2a,
                       model.name.1="Absence",
                       model.name.2="Presence",
                       class.1="0",
                       class.2="1",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Variable Importance",
                       device.type="none",
                       cex=0.9)
	

model.importance.plot( model.obj.1=model.obj.ex2b,
                       model.obj.2=model.obj.ex2b,
                       model.name.1="Absence",
                       model.name.2="Presence",
                       class.1="0",
                       class.2="1",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Variable Importance",
                       device.type="none",
                       cex=0.9)	                       
mtext("Presence-Absence Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 52: Ex2CompImpPAFig
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex2a,
                       model.obj.2=model.obj.ex2a,
                       model.name.1="Absence",
                       model.name.2="Presence",
                       class.1="0",
                       class.2="1",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Variable Importance",
                       device.type="none",
                       cex=0.9)
	

model.importance.plot( model.obj.1=model.obj.ex2b,
                       model.obj.2=model.obj.ex2b,
                       model.name.1="Absence",
                       model.name.2="Presence",
                       class.1="0",
                       class.2="1",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Variable Importance",
                       device.type="none",
                       cex=0.9)	                       
mtext("Presence-Absence Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 53: Ex2 Interaction Plots
###################################################
model.interaction.plot( model.obj.ex2a,
                        x="ELEV250",
                        y="NLCD01_250",
                        main=response.name.a,
                        plot.type="image",
                        device.type="pdf",
                        MODELfn=MODELfn.a,
                        folder=folder)


model.interaction.plot( model.obj.ex2b,
                        x="ELEV250",
                        y="NLCD01_250",
                        main=response.name.b,
                        plot.type="image",
                        device.type="pdf",
                        MODELfn=MODELfn.b,
                        folder=folder)	
																									
model.interaction.plot( model.obj.ex2a,
                        x="ELEV250",
                        y="NLCD01_250",
                        main=response.name.a,
                        plot.type="persp",
                        device.type="pdf",
                        MODELfn=MODELfn.a,
                        folder=folder,
                        theta=300,
                        phi=55)


model.interaction.plot( model.obj.ex2b,
                        x="ELEV250",
                        y="NLCD01_250",
                        main=response.name.b,
                        plot.type="persp",
                        device.type="pdf",
                        MODELfn=MODELfn.b,
                        folder=folder,
                        theta=300,
                        phi=55)


###################################################
### code chunk number 54: Ex2 Produce Maps
###################################################
model.mapmake( model.obj=model.obj.ex2a,
               folder=folder,           
               MODELfn=MODELfn.a,
               rastLUTfn=rastLUTfn,
               na.action="na.omit")
            

model.mapmake( model.obj=model.obj.ex2b,
               folder=folder,           
               MODELfn=MODELfn.b,
               rastLUTfn=rastLUTfn,
               na.action="na.omit")



###################################################
### code chunk number 55: Ex2 define color sequence
###################################################
h=c( seq(10,30,length.out=10),
      seq(31,40,length.out=10),
      seq(41,90,length.out=60),
      seq(91,100,length.out=10),
      seq(101,110,length.out=10))
l =c( seq(25,40,length.out=10),
      seq(40,90,length.out=35),
      seq(90,90,length.out=10),
      seq(90,40,length.out=35),
      seq(40,10,length.out=10))
probpres.ramp <- hcl(h = h, c = 80, l = l)


###################################################
### code chunk number 56: Ex2ProbabilitySurface
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <- raster(paste(MODELfn.a,"_map.img",sep=""))
mapgrid.b <- raster(paste(MODELfn.b,"_map.img",sep=""))

legend.subset<-c(100,80,60,40,20,1)
legend.colors<-probpres.ramp[legend.subset]
legend.label<-c("100%"," 80%"," 60%"," 40%"," 20%","  0%")

image( mapgrid.a, 
       col=probpres.ramp,
       xlab="",ylab="",yaxt="n",main="",zlim=c(0,1),
       asp=1,bty="n",xaxt="n")
mtext(response.name.a,side=3,line=1,cex=1.2)
image( mapgrid.b, 
       col=probpres.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=c(0,1),
       asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
    legend=legend.label,
    fill=legend.colors,
    bty="n",
    cex=1.2)
				    
mtext("Probability of Presence",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 57: Ex2ProbabilitySurfaceFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <- raster(paste(MODELfn.a,"_map.img",sep=""))
mapgrid.b <- raster(paste(MODELfn.b,"_map.img",sep=""))

legend.subset<-c(100,80,60,40,20,1)
legend.colors<-probpres.ramp[legend.subset]
legend.label<-c("100%"," 80%"," 60%"," 40%"," 20%","  0%")

image( mapgrid.a, 
       col=probpres.ramp,
       xlab="",ylab="",yaxt="n",main="",zlim=c(0,1),
       asp=1,bty="n",xaxt="n")
mtext(response.name.a,side=3,line=1,cex=1.2)
image( mapgrid.b, 
       col=probpres.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=c(0,1),
       asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
    legend=legend.label,
    fill=legend.colors,
    bty="n",
    cex=1.2)
				    
mtext("Probability of Presence",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 58: Ex2aPresenceMaps
###################################################
opar <- par(mfrow=c(2,2),mar=c(2.5,3,4,1),oma=c(0,0,4,6),xpd=NA)
mapgrid <- raster(paste(MODELfn.a,"_map.img",sep=""))
criteria <- c("Default","MaxKappa","ReqSens","ReqSpec")
criteria.labels<-c("Default","MaxKappa","ReqSens = 0.9","ReqSpec = 0.9")
for(i in 1:4){
      thresh <- opt.thresh.a$threshold[opt.thresh.a$opt.methods==criteria[i]]
      presencegrid <- mapgrid
      v <- getValues(presencegrid)
	v <- ifelse(v > thresh,1,0) 
      presencegrid <- setValues(presencegrid, v)

	image( presencegrid,
	       col=c("white","forestgreen"),
	       zlim=c(0,1),
	       asp=1,
	       bty="n",
	       xaxt="n", yaxt="n",
             main="",xlab="",ylab="")
	if(i==2){
		legend( x=xmax(mapgrid),y=ymax(mapgrid),
		        legend=c("Present","Absent"),
		        fill=c("forestgreen","white"),
		        bty="n",
		        cex=1.2)}
	mtext(criteria.labels[i],side=3,line=2,cex=1.2)
	mtext(paste("threshold =",thresh),side=3,line=.5,cex=1)
}
mtext(MODELfn.a,side=3,line=0,cex=1.2,outer=TRUE)
mtext(response.name.a,side=3,line=2,cex=1.5,outer=TRUE)
par(opar)


###################################################
### code chunk number 59: Ex2aPresenceMapsFig
###################################################
opar <- par(mfrow=c(2,2),mar=c(2.5,3,4,1),oma=c(0,0,4,6),xpd=NA)
mapgrid <- raster(paste(MODELfn.a,"_map.img",sep=""))
criteria <- c("Default","MaxKappa","ReqSens","ReqSpec")
criteria.labels<-c("Default","MaxKappa","ReqSens = 0.9","ReqSpec = 0.9")
for(i in 1:4){
      thresh <- opt.thresh.a$threshold[opt.thresh.a$opt.methods==criteria[i]]
      presencegrid <- mapgrid
      v <- getValues(presencegrid)
	v <- ifelse(v > thresh,1,0) 
      presencegrid <- setValues(presencegrid, v)

	image( presencegrid,
	       col=c("white","forestgreen"),
	       zlim=c(0,1),
	       asp=1,
	       bty="n",
	       xaxt="n", yaxt="n",
             main="",xlab="",ylab="")
	if(i==2){
		legend( x=xmax(mapgrid),y=ymax(mapgrid),
		        legend=c("Present","Absent"),
		        fill=c("forestgreen","white"),
		        bty="n",
		        cex=1.2)}
	mtext(criteria.labels[i],side=3,line=2,cex=1.2)
	mtext(paste("threshold =",thresh),side=3,line=.5,cex=1)
}
mtext(MODELfn.a,side=3,line=0,cex=1.2,outer=TRUE)
mtext(response.name.a,side=3,line=2,cex=1.5,outer=TRUE)
par(opar)


###################################################
### code chunk number 60: Ex2bPresenceMaps
###################################################
opar <- par(mfrow=c(2,2),mar=c(2.5,3,4,1),oma=c(0,0,4,6),xpd=NA)
mapgrid <- raster(paste(MODELfn.b,"_map.img",sep=""))
criteria <- c("Default","MaxKappa","ReqSens","ReqSpec")
criteria.labels<-c("Default","MaxKappa","ReqSens = 0.9","ReqSpec = 0.9")
for(i in 1:4){
      thresh <- opt.thresh.b$threshold[opt.thresh.b$opt.methods==criteria[i]]
      presencegrid <- mapgrid
      v <- getValues(presencegrid)
	v <- ifelse(v > thresh,1,0) 
      presencegrid <- setValues(presencegrid, v)

	image( presencegrid,
	       col=c("white","forestgreen"),
	       xlab="",ylab="",xaxt="n", yaxt="n",
	       zlim=c(0,1),
	       asp=1,bty="n",main="")
	if(i==2){
		legend( x=xmax(mapgrid),y=ymax(mapgrid),
		        legend=c("Present","Absent"),
		        fill=c("forestgreen","white"),
		        bty="n",
		        cex=1.2)}
	mtext(criteria.labels[i],side=3,line=2,cex=1.2)
	mtext(paste("threshold =",thresh),side=3,line=.5,cex=1)
}
mtext(MODELfn.b,side=3,line=0,cex=1.2,outer=TRUE)
mtext(response.name.b,side=3,line=2,cex=1.5,outer=TRUE)
par(opar)


###################################################
### code chunk number 61: Ex2bPresenceMapsFig
###################################################
opar <- par(mfrow=c(2,2),mar=c(2.5,3,4,1),oma=c(0,0,4,6),xpd=NA)
mapgrid <- raster(paste(MODELfn.b,"_map.img",sep=""))
criteria <- c("Default","MaxKappa","ReqSens","ReqSpec")
criteria.labels<-c("Default","MaxKappa","ReqSens = 0.9","ReqSpec = 0.9")
for(i in 1:4){
      thresh <- opt.thresh.b$threshold[opt.thresh.b$opt.methods==criteria[i]]
      presencegrid <- mapgrid
      v <- getValues(presencegrid)
	v <- ifelse(v > thresh,1,0) 
      presencegrid <- setValues(presencegrid, v)

	image( presencegrid,
	       col=c("white","forestgreen"),
	       xlab="",ylab="",xaxt="n", yaxt="n",
	       zlim=c(0,1),
	       asp=1,bty="n",main="")
	if(i==2){
		legend( x=xmax(mapgrid),y=ymax(mapgrid),
		        legend=c("Present","Absent"),
		        fill=c("forestgreen","white"),
		        bty="n",
		        cex=1.2)}
	mtext(criteria.labels[i],side=3,line=2,cex=1.2)
	mtext(paste("threshold =",thresh),side=3,line=.5,cex=1)
}
mtext(MODELfn.b,side=3,line=0,cex=1.2,outer=TRUE)
mtext(response.name.b,side=3,line=2,cex=1.5,outer=TRUE)
par(opar)


###################################################
### code chunk number 62: Ex3 Define model type
###################################################
model.type <- "RF"


###################################################
### code chunk number 63: Ex3 Define training and test files
###################################################
qdatafn <- "VModelMapData.csv"


###################################################
### code chunk number 64: Ex3 define folder
###################################################
folder <- getwd()


###################################################
### code chunk number 65: Ex3 Define model filename
###################################################
MODELfn    <- "VModelMapEx3"


###################################################
### code chunk number 66: Ex3 Define predictors
###################################################
predList <- c( "ELEV250",
               "NLCD01_250",
               "EVI2005097",
               "NDV2005097",
               "NIR2005097",
               "RED2005097")
predFactor <- c("NLCD01_250")


###################################################
### code chunk number 67: Ex3 Define response
###################################################
response.name <- "VEGCAT"
response.type <- "categorical"


###################################################
### code chunk number 68: Ex3 set seed
###################################################
seed <- 44


###################################################
### code chunk number 69: Ex3 Define Identifier
###################################################
unique.rowname <- "ID"


###################################################
### code chunk number 70: Ex3 update raster LUT
###################################################
rastLUTfn     <- "VModelMapData_LUT.csv"
rastLUTfn     <- read.table( rastLUTfn,
                             header=FALSE,
                             sep=",",
                             stringsAsFactors=FALSE)
rastLUTfn[,1] <- paste(folder,rastLUTfn[,1],sep="/")


###################################################
### code chunk number 71: Ex3 Create Model
###################################################
model.obj.ex3 <- model.build( model.type=model.type,
                               qdata.trainfn=qdatafn,
                               folder=folder,
                               unique.rowname=unique.rowname,               
                               MODELfn=MODELfn,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name,
                               response.type=response.type,
                               seed=seed)
           


###################################################
### code chunk number 72: Ex3 Model Diagnostics
###################################################
model.pred.ex3 <- model.diagnostics( model.obj=model.obj.ex3,
                                      qdata.trainfn=qdatafn,
                                      folder=folder,           
                                      MODELfn=MODELfn,
                                      unique.rowname=unique.rowname,
                             # Model Validation Arguments
                                      prediction.type="OOB",
                                      device.type="pdf",
                                      cex=1.2)
           


###################################################
### code chunk number 73: Ex3 Confusion Matrix Show
###################################################
CMX.CSV <- read.table(   paste( MODELfn,"_pred_cmx.csv",sep=""),
                                header=FALSE,
                                sep=",",
                                stringsAsFactors=FALSE)

CMX.CSV


###################################################
### code chunk number 74: Ex3 Observed and Predicted Show
###################################################
PRED <-read.table(   paste( MODELfn,"_pred.csv",sep=""),
                            header=TRUE,
                            sep=",",
                            stringsAsFactors=TRUE)
head(PRED)


###################################################
### code chunk number 75: Ex3 Prevalence Table Show
###################################################
#
#these lines are needed for numeric categories, redundant for character categories
#
PRED$pred<-as.factor(PRED$pred)
PRED$obs<-as.factor(PRED$obs)
#
#adjust levels so all values are included in both observed and predicted
#
LEVELS<-unique(c(levels(PRED$pred),levels(PRED$obs)))
PRED$pred<-factor(PRED$pred,levels=LEVELS)
PRED$obs<- factor(PRED$obs, levels=LEVELS)
#
#calculate confusion matrix
#
CMX<-table( predicted=PRED$pred, observed= PRED$obs)
CMX


###################################################
### code chunk number 76: Ex 3 Marginals Show
###################################################
CMX.diag  <- diag(CMX)

CMX.OMISSION  <- 1-(CMX.diag/apply(CMX,2,sum))
CMX.COMISSION <- 1-(CMX.diag/apply(CMX,1,sum))

CMX.OMISSION
CMX.COMISSION


###################################################
### code chunk number 77: Ex3 PCC Show
###################################################
CMX.PCC <- sum(CMX.diag)/sum(CMX)
CMX.PCC


###################################################
### code chunk number 78: Ex3 Kappa Show
###################################################
CMX.KAPPA <- PresenceAbsence::Kappa(CMX)
CMX.KAPPA


###################################################
### code chunk number 79: Ex3 MAUC Show
###################################################
VOTE <- HandTill2001::multcap(  response = PRED$obs,
                	predicted= as.matrix(PRED[,-c(1,2,3)]) )
MAUC  <- HandTill2001::auc(VOTE)
MAUC


###################################################
### code chunk number 80: Ex3CompImp
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex2a,
                       model.obj.2=model.obj.ex3,
                       model.name.1="Pinyon",
                       model.name.2="VEGCAT",
                       type.label=FALSE,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Presence vs VEGCAT",
                       device.type="none",
                       cex=0.9)
                      
model.importance.plot( model.obj.1=model.obj.ex2b,
                       model.obj.2=model.obj.ex3,
                       model.name.1="Sage",
                       model.name.2="VEGCAT",
                       type.label=FALSE,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Presence vs VEGCAT",
                       device.type="none",
                       cex=0.9)
                       
mtext("Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 81: Ex3CompImpFig
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex2a,
                       model.obj.2=model.obj.ex3,
                       model.name.1="Pinyon",
                       model.name.2="VEGCAT",
                       type.label=FALSE,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Presence vs VEGCAT",
                       device.type="none",
                       cex=0.9)
                      
model.importance.plot( model.obj.1=model.obj.ex2b,
                       model.obj.2=model.obj.ex3,
                       model.name.1="Sage",
                       model.name.2="VEGCAT",
                       type.label=FALSE,
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Presence vs VEGCAT",
                       device.type="none",
                       cex=0.9)
                       
mtext("Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 82: Ex3CompCatImp
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex3,
                       model.obj.2=model.obj.ex3,
                       model.name.1="SHRUB",
                       model.name.2="TREE",
                       class.1="SHRUB",
                       class.2="TREE",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="VEGCAT - SHRUB vs. TREE",
                       device.type="none",
                       cex=0.9) 
                      
model.importance.plot( model.obj.1=model.obj.ex3,
                       model.obj.2=model.obj.ex3,
                       model.name.1="OTHERVEG",
                       model.name.2="NONVEG",
                       class.1="OTHERVEG",
                       class.2="NONVEG",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="VEGCAT - OTHERVEG vs. NONVEG",
                       device.type="none",
                       cex=0.9)
                       
mtext("Category Specific Variable Importance",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 83: Ex3CompCatImpFig
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex3,
                       model.obj.2=model.obj.ex3,
                       model.name.1="SHRUB",
                       model.name.2="TREE",
                       class.1="SHRUB",
                       class.2="TREE",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="VEGCAT - SHRUB vs. TREE",
                       device.type="none",
                       cex=0.9) 
                      
model.importance.plot( model.obj.1=model.obj.ex3,
                       model.obj.2=model.obj.ex3,
                       model.name.1="OTHERVEG",
                       model.name.2="NONVEG",
                       class.1="OTHERVEG",
                       class.2="NONVEG",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="VEGCAT - OTHERVEG vs. NONVEG",
                       device.type="none",
                       cex=0.9)
                       
mtext("Category Specific Variable Importance",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 84: Ex3 Interaction Plots
###################################################
model.interaction.plot( model.obj.ex3,
                        x="ELEV250",
                        y="NLCD01_250",
                        main=response.name,
                        plot.type="image",
                        device.type="pdf",
                        MODELfn=MODELfn,
                        folder=folder,
                        response.category="SHRUB")

	
																									
model.interaction.plot( model.obj.ex3,
                        x="ELEV250",
                        y="NLCD01_250",
                        main=response.name,
                        plot.type="image",
                        device.type="pdf",
                        MODELfn=MODELfn,
                        folder=folder,
                        response.category="NONVEG")



###################################################
### code chunk number 85: Ex3 Produce Maps
###################################################

model.mapmake( model.obj=model.obj.ex3,
               folder=folder,           
               MODELfn=MODELfn,
               rastLUTfn=rastLUTfn,
               na.action="na.omit")            


###################################################
### code chunk number 86: Ex3 Code key Show
###################################################

MAP.CODES<-read.table( paste(MODELfn,"_map_key.csv",sep=""),
                       header=TRUE,
                       sep=",",
                       stringsAsFactors=FALSE)

MAP.CODES           


###################################################
### code chunk number 87: Ex3 define color sequence
###################################################

MAP.CODES$colors<-c("bisque3","springgreen2","paleturquoise1","green4")
MAP.CODES


###################################################
### code chunk number 88: Ex3 Import data
###################################################

mapgrid <- raster(paste(MODELfn,"_map.img",sep=""))

integergrid <- mapgrid
v <- getValues(mapgrid)
v <- MAP.CODES$row[match(v,MAP.CODES$integercode)]
integergrid <- setValues(integergrid, v)


###################################################
### code chunk number 89: Ex3ProductionMaps
###################################################


opar <- par(mfrow=c(1,1),mar=c(3,3,2,1),oma=c(0,0,3,8),xpd=NA)

image(	integergrid, 
				col = MAP.CODES$colors,
				xlab="",ylab="",xaxt="n",yaxt="n",
				zlim=c(1,nrow(MAP.CODES)),
				main="",asp=1,bty="n")

mtext(response.name,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid),y=ymax(mapgrid),
        legend=MAP.CODES$category,
        fill=MAP.CODES$colors,
        bty="n",
        cex=1.2)
				    
par(opar)


###################################################
### code chunk number 90: Ex3ProductionMapsFig
###################################################


opar <- par(mfrow=c(1,1),mar=c(3,3,2,1),oma=c(0,0,3,8),xpd=NA)

image(	integergrid, 
				col = MAP.CODES$colors,
				xlab="",ylab="",xaxt="n",yaxt="n",
				zlim=c(1,nrow(MAP.CODES)),
				main="",asp=1,bty="n")

mtext(response.name,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid),y=ymax(mapgrid),
        legend=MAP.CODES$category,
        fill=MAP.CODES$colors,
        bty="n",
        cex=1.2)
				    
par(opar)


###################################################
### code chunk number 91: Ex4 Define model type
###################################################
model.type <- "SGB"


###################################################
### code chunk number 92: Ex4 Define training and test files
###################################################
qdatafn <- "VModelMapData.csv"


###################################################
### code chunk number 93: Ex4 define folder
###################################################
folder <- getwd()


###################################################
### code chunk number 94: Ex4 Define model filename
###################################################
MODELfn.a    <- "VModelMapEx4a"
MODELfn.b    <- "VModelMapEx4b"


###################################################
### code chunk number 95: Ex4 Define predictors
###################################################
predList <- c( "ELEV250",
               "NLCD01_250",
               "EVI2005097",
               "NDV2005097",
               "NIR2005097",
               "RED2005097")
predFactor <- c("NLCD01_250")


###################################################
### code chunk number 96: Ex4 Define response
###################################################
response.name.a <- "PINYON"
response.name.b <- "SAGE"
response.type <- "binary"


###################################################
### code chunk number 97: Ex4 set seed
###################################################
seed.a <- 42
seed.b <- 43


###################################################
### code chunk number 98: Ex4 Define Identifier
###################################################
unique.rowname <- "ID"


###################################################
### code chunk number 99: Ex4 update raster LUT
###################################################
rastLUTfn   <- "VModelMapData_LUT.csv"
rastLUTfn     <- read.table( rastLUTfn,
                             header=FALSE,
                             sep=",",
                             stringsAsFactors=FALSE)
rastLUTfn[,1] <- paste(folder,rastLUTfn[,1],sep="/")


###################################################
### code chunk number 100: Ex4 Create Model
###################################################
model.obj.ex4a <- model.build( model.type=model.type,
                               qdata.trainfn=qdatafn,
                               folder=folder,
                               unique.rowname=unique.rowname,               
                               MODELfn=MODELfn.a,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.a,
                               response.type=response.type,
                               seed=seed.a)
           
model.obj.ex4b <- model.build( model.type=model.type,
                               qdata.trainfn=qdatafn,
                               folder=folder,
                               unique.rowname=unique.rowname,               
                               MODELfn=MODELfn.b,
                               predList=predList,
                               predFactor=predFactor,
                               response.name=response.name.b,
                               response.type=response.type,
                               seed=seed.b)


###################################################
### code chunk number 101: Ex4 Model Diagnostics
###################################################
model.pred.ex4a <- model.diagnostics( model.obj=model.obj.ex4a,
                                      qdata.trainfn=qdatafn,
                                      folder=folder,           
                                      MODELfn=MODELfn.a,
                                      unique.rowname=unique.rowname,
                                      seed=44,
                             # Model Validation Arguments
                                      prediction.type="CV",
                                      device.type=c("jpeg","pdf","postscript"),
                                      cex=1.2,
                                      na.action = "na.roughfix")
           
model.pred.ex4b <- model.diagnostics( model.obj=model.obj.ex4b,
                                      qdata.trainfn=qdatafn,
                                      folder=folder,           
                                      MODELfn=MODELfn.b,
                                      unique.rowname=unique.rowname,
                                      seed=45,
                              # Model Validation Arguments
                                      prediction.type="CV",
                                      device.type=c("jpeg","pdf","postscript"),
                                      cex=1.2,
                                      na.action = "na.roughfix")


###################################################
### code chunk number 102: Ex4 Optimal Threshold Table
###################################################
opt.thresh.a <- read.table( paste(MODELfn.a,"_pred_optthresholds.csv",sep=""),
                            header=TRUE,
                            sep=",",
                            stringsAsFactors=FALSE)
opt.thresh.a[,-1]<-signif(opt.thresh.a[,-1],2)

opt.thresh.b <- read.table( paste(MODELfn.b,"_pred_optthresholds.csv",sep=""),
                            header=TRUE,
                            sep=",",
                            stringsAsFactors=FALSE)
opt.thresh.b[,-1]<-signif(opt.thresh.b[,-1],2)


###################################################
### code chunk number 103: Ex4a Optimal Threshold Table Show
###################################################
opt.thresh.a


###################################################
### code chunk number 104: Ex4b Optimal Threshold Table Show
###################################################
opt.thresh.b


###################################################
### code chunk number 105: Ex4CompImp
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex2a,
                       model.obj.2=model.obj.ex4a,
                       model.name.1="RF",
                       model.name.2="SGB",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Presence",
                       device.type="none",
                       cex=0.9)
                      
model.importance.plot( model.obj.1=model.obj.ex2b,
                       model.obj.2=model.obj.ex4b,
                       model.name.1="RF",
                       model.name.2="SGB",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Presence",
                       device.type="none",
                       cex=0.9)
                       
mtext("Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 106: Ex4CompImpFig
###################################################

opar <- par(mfrow=c(2,1),mar=c(3,3,3,3),oma=c(0,0,3,0))

model.importance.plot( model.obj.1=model.obj.ex2a,
                       model.obj.2=model.obj.ex4a,
                       model.name.1="RF",
                       model.name.2="SGB",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Pinyon Presence",
                       device.type="none",
                       cex=0.9)
                      
model.importance.plot( model.obj.1=model.obj.ex2b,
                       model.obj.2=model.obj.ex4b,
                       model.name.1="RF",
                       model.name.2="SGB",
                       sort.by="predList",
                       predList=predList,
                       scale.by="sum",
                       main="Sage Presence",
                       device.type="none",
                       cex=0.9)
                       
mtext("Variable Importance Comparison",side=3,line=0,cex=1.8,outer=TRUE)
par(opar)


###################################################
### code chunk number 107: Ex4 Interaction Plots
###################################################
ELEV<-c(1600,2100)
NCLD<-c(30,80)

for(j in ELEV){
for(k in NCLD){
      pred.means<-list(j,k)
      names(pred.means)<-c("ELEV250","NLCD01_250")
      model.interaction.plot( model.obj.ex4b,
                              x="EVI2005097",
                              y="NDV2005097",
                              main=paste(response.name.b," (ELEV=",j,"  NCLD=",k,")",sep=""),
                              plot.type="persp",
                              device.type="pdf",
                              pred.means=pred.means,
                              theta=65,phi=25,
                              MODELfn=paste(MODELfn.b,"ELEV",j,"NCLD",k,sep="_"),
                              folder=folder)
}}



###################################################
### code chunk number 108: Ex4 Produce Maps
###################################################

model.mapmake( model.obj=model.obj.ex4a,
               folder=folder,           
               MODELfn=MODELfn.a,
               rastLUTfn=rastLUTfn)
            

model.mapmake( model.obj=model.obj.ex4b,
               folder=folder,           
               MODELfn=MODELfn.b,
               rastLUTfn=rastLUTfn)



###################################################
### code chunk number 109: Ex4 define color sequence
###################################################
h=c(  seq(10,30,length.out=10),
	    seq(31,40,length.out=10),
    	seq(41,90,length.out=60),
	    seq(91,100,length.out=10),
	    seq(101,110,length.out=10))
l =c( seq(25,40,length.out=10),
	    seq(40,90,length.out=35),
	    seq(90,90,length.out=10),
	    seq(90,40,length.out=35),
	    seq(40,10,length.out=10))
probpres.ramp <- hcl(h = h, c = 80, l = l)


###################################################
### code chunk number 110: Ex4ProbabilitySurface
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <- raster(paste(MODELfn.a,"_map.img",sep=""))
mapgrid.b <- raster(paste(MODELfn.b,"_map.img",sep=""))

legend.subset<-c(100,80,60,40,20,1)
legend.colors<-probpres.ramp[legend.subset]
legend.label<-c("100%"," 80%"," 60%"," 40%"," 20%","  0%")

image( mapgrid.a, 
       col=probpres.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=c(0,1),
       asp=1,bty="n",main="")
mtext(response.name.a,side=3,line=1,cex=1.2)
image(mapgrid.b, 
      col=probpres.ramp,
      xlab="",ylab="",xaxt="n",yaxt="n",
      zlim=c(0,1),
      asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
        legend=legend.label,
        fill=legend.colors,
        bty="n",
        cex=1.2)
				    
mtext("Probability of Presence",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 111: Ex4ProductionMapsFig
###################################################
opar <- par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(0,0,3,4),xpd=NA)

mapgrid.a <- raster(paste(MODELfn.a,"_map.img",sep=""))
mapgrid.b <- raster(paste(MODELfn.b,"_map.img",sep=""))

legend.subset<-c(100,80,60,40,20,1)
legend.colors<-probpres.ramp[legend.subset]
legend.label<-c("100%"," 80%"," 60%"," 40%"," 20%","  0%")

image( mapgrid.a, 
       col=probpres.ramp,
       xlab="",ylab="",xaxt="n",yaxt="n",
       zlim=c(0,1),
       asp=1,bty="n",main="")
mtext(response.name.a,side=3,line=1,cex=1.2)
image(mapgrid.b, 
      col=probpres.ramp,
      xlab="",ylab="",xaxt="n",yaxt="n",
      zlim=c(0,1),
      asp=1,bty="n",main="")
mtext(response.name.b,side=3,line=1,cex=1.2)

legend( x=xmax(mapgrid.b),y=ymax(mapgrid.b),
        legend=legend.label,
        fill=legend.colors,
        bty="n",
        cex=1.2)
				    
mtext("Probability of Presence",side=3,line=1,cex=1.5,outer=T)
par(opar)


###################################################
### code chunk number 112: Remove Rplots pdf
###################################################
dev.off()
file.remove("Rplots.pdf")
file.remove("Vplots.pdf")



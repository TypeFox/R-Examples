### R code from vignette source 'ExploreModelMap.Rnw'

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
### code chunk number 4: ExSetup Define training and test files
###################################################
qdatafn <- "VModelMapData.csv"
qdata.trainfn <- "VModelMapData_TRAIN.csv"
qdata.testfn  <- "VModelMapData_TEST.csv"


###################################################
### code chunk number 5: ExSetup define folder
###################################################
folder <- getwd()


###################################################
### code chunk number 6: ExSetup split training and test
###################################################
get.test(       proportion.test=0.2,
                qdatafn=qdatafn,
                seed=42,
                folder=folder,
                qdata.trainfn=qdata.trainfn,
                qdata.testfn=qdata.testfn)


###################################################
### code chunk number 7: ExSetup Define predictors
###################################################
predList <- c( "ELEV250",
               "NLCD01_250",
               "EVI2005097",
               "NDV2005097",
               "NIR2005097",
               "RED2005097")
predFactor <- c("NLCD01_250")


###################################################
### code chunk number 8: Ex1 Define Identifier
###################################################
unique.rowname <- "ID"


###################################################
### code chunk number 9: ExSetup update raster LUT
###################################################
rastLUTfn     <- "VModelMapData_LUT.csv"
rastLUTfn     <- read.table( rastLUTfn,
                             header=FALSE,
                             sep=",",
                             stringsAsFactors=FALSE)
rastLUTfn[,1] <- paste(folder,rastLUTfn[,1],sep="/")


###################################################
### code chunk number 10: ExCorr
###################################################

qdata.train <- read.table(file=qdata.trainfn,sep=",",header=TRUE,check.names=FALSE,as.is=TRUE)
correlation.function(	qdata=qdata.train,
                      predList=predList,
                      predFactor=predFactor,
                      MODELpredfn=paste(folder,"Explore",sep="/"),
                      device.type=c("jpeg","pdf","png"),
											cex=1
                      )


###################################################
### code chunk number 11: Ex1a Model Explore Pinyon
###################################################
model.explore( qdata.trainfn=qdata.trainfn,
               folder=folder,		
               predList=predList,
               predFactor=predFactor,
							
               OUTPUTfn="PinyonCover",

               response.name="PINYON",
               response.type="continuous",
	
               unique.rowname=unique.rowname,

               device.type=c("png"),
               #cex=1.2,

               # Raster arguments
               rastLUTfn=rastLUTfn,
               na.value=-9999,

               # colors for continuous predictors
               col.ramp=terrain.colors(101),
							
               #colors for categorical predictors
               col.cat=c("wheat1","springgreen2","darkolivegreen4",
                         "darkolivegreen2","yellow","thistle2",
                         "brown2","brown4")
)


###################################################
### code chunk number 12: MapElevReadGDAL
###################################################

elevfn <- paste(folder,"/VModelMapData_dem_ELEVM_250.img",sep="")
mapgrid <- raster(elevfn)


###################################################
### code chunk number 13: MapElev
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
### code chunk number 14: MapElevFig
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
### code chunk number 15: Ex1a Model Explore Sage Presence
###################################################
model.explore( qdata.trainfn=qdata.trainfn,
               folder=folder,		
               predList=predList,
               predFactor=predFactor,
							
               OUTPUTfn="SagePresence",

               response.name="SAGE",
               response.type="binary",
	
               unique.rowname=unique.rowname,

               device.type=c("png"),
               #cex=1.2,

               # Raster arguments
               rastLUTfn=rastLUTfn,
               na.value=-9999,

               # colors for continuous predictors
               col.ramp=heat.colors(101),
							
               #colors for categorical predictors
               col.cat=c("wheat1","springgreen2","darkolivegreen4",
                         "darkolivegreen2","yellow","thistle2",
                         "brown2","brown4")
)


###################################################
### code chunk number 16: Ex1a Model Explore VEGCAT
###################################################
model.explore( qdata.trainfn=qdata.trainfn,
               folder=folder,		
               predList=predList,
               predFactor=predFactor,
							
               OUTPUTfn="VegCat",

               response.name="VEGCAT",
               response.type="categorical",
	
               unique.rowname=unique.rowname,

               device.type=c("png"),
               #cex=1.2,

               # Raster arguments
               rastLUTfn=rastLUTfn,
               na.value=-9999,

               # colors for continuous predictors
               col.ramp=heat.colors(101),
							
               #colors for categorical predictors
               col.cat=c("wheat1","springgreen2","darkolivegreen4",
                         "darkolivegreen2","yellow","thistle2",
                         "brown2","brown4")
)


###################################################
### code chunk number 17: Remove Rplots pdf
###################################################
dev.off()

file.remove("Vplots.pdf")



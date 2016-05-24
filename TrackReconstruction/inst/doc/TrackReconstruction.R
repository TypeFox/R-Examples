### R code from vignette source 'TrackReconstruction.Snw'

###################################################
### code chunk number 1: TrackReconstruction.Snw:325-329
###################################################
 library(TrackReconstruction)
 betas<-Standardize(1,1,-1,1,1,1,-57.8,68.76,-61.8,64.2,-70.16,58.08,-10.1,9.55,
 -9.75,9.72,-9.91,9.43)
 betas


###################################################
### code chunk number 2: TrackReconstruction.Snw:340-343
###################################################
 data(rawdatagap)
 gaps<-GapFinder(rawdatagap, timediff = 1, timeformat = "%d-%b-%Y %H:%M:%S")
 head(gaps)


###################################################
### code chunk number 3: TrackReconstruction.Snw:354-357
###################################################
 data(gpsdata02)
 gpsformat<-GPStable(gpsdata02)
 head(gpsformat)


###################################################
### code chunk number 4: TrackReconstruction.Snw:369-376
###################################################
 #get declination and inclination data for study area
 decinc<-c(10.228,65.918)
 #the fifth of 10 consecutive data sets from tag data from one individual with start and end 
 #times corresponding to GPS fixes (+ rows=Hz*RmL*0.5 on each end).
 data(rawdata)
 DRoutput<-DeadReckoning(rawdata, betas, decinc, Hz = 16, RmL = 2, DepthHz = 1, SpdCalc=3,
 MaxSpd=3.5)


###################################################
### code chunk number 5: TrackReconstruction.Snw:399-400
###################################################
 Georeferenced<-GeoReference(DRoutput,gpsformat[5:6,])


###################################################
### code chunk number 6: TrackReconstruction.Snw:444-462
###################################################
 require(scatterplot3d)
 require(onion)
 require(RColorBrewer)
 require(lattice)
 library(plotrix)
 library(fields)
 require(rgl)
 library(TrackReconstruction)
 #Import data
 #setwd("G:\\filepath\\gebco_08")
 #bathymetry<- read.table("Gebco1.asc",sep=",",header=TRUE)
 #or get the example data from the package
 data(bathymetry)
 col=gray(0:200/200)
 #format data for graphing
 image.xyz=tapply(bathymetry$Depth, list(bathymetry$Long, bathymetry$Lat), unique)
 #create palatte for depth colors
 Bathymetry.palatte<-colorRampPalette(brewer.pal(9, "Blues"),bias=3)


###################################################
### code chunk number 7: fig1plot
###################################################
 	image.plot(x=as.numeric(dimnames(image.xyz)[[1]]),
 	y=as.numeric(dimnames(image.xyz)[[2]]), 
 	z=image.xyz,
	col=c(rev(Bathymetry.palatte(100)),#gray(0:20/20),
	terrain.colors(100)),
	breaks=round(c(seq(from=min(image.xyz),to=0,length.out=101),
	seq(from=max(image.xyz)/101,to=max(image.xyz),length.out=100))),
	ylab="",
 	xlab=""
	#,smallplot=2 #plots legend off x axis
 	)


###################################################
### code chunk number 8: TrackReconstruction.Snw:498-504
###################################################
 data(georef1min01)
 limits<-GraphLimits(georef1min01)
 	Sminlat=limits$miny
 	Smaxlat=limits$maxy
 	Sminlong=limits$minx
 	Smaxlong=limits$maxx


###################################################
### code chunk number 9: fig2plot
###################################################
 scatterplot3d(georef1min01$Longitude,georef1min01$Latitude,(georef1min01$Depth*-1),
 	color="black",#ifelse(georef1min01$SunTimes==1,"red","black"), shades night
 		#and day if you have the data
 	type="l",
 	lwd=1,
 	#pch=".",
 	highlight.3d=F,
 	angle=55,
 	xlim=c(Sminlong,Smaxlong),
 	ylim=c(Sminlat,Smaxlat),
 	zlim=c(0,-80),
 	zlab="Depth",
 	ylab="Latitude",
 	xlab="Longitude",
 	#x.ticklabs=round(seq(from=Sminlong,to=Smaxlong, by=(Smaxlong-Sminlong)/4),digits=2),
 	#y.ticklabs=round(seq(from=Sminlat,to=Smaxlat, by=(Smaxlat-Sminlat)/4),digits=2),
 	#z.ticklabs=c(-80,-60,-40,-20,0),
 	cex.lab=1,
 	cex.axis=1,
 	cex.symbols=1,
 	#lab=c(3, 4),
 	lab.z=5
 )


###################################################
### code chunk number 10: TrackReconstruction.Snw:537-538
###################################################
 scatterplot3d(georef1min01$Longitude,georef1min01$Latitude,(georef1min01$Depth*-1),
 	color="black",#ifelse(georef1min01$SunTimes==1,"red","black"), shades night
 		#and day if you have the data
 	type="l",
 	lwd=1,
 	#pch=".",
 	highlight.3d=F,
 	angle=55,
 	xlim=c(Sminlong,Smaxlong),
 	ylim=c(Sminlat,Smaxlat),
 	zlim=c(0,-80),
 	zlab="Depth",
 	ylab="Latitude",
 	xlab="Longitude",
 	#x.ticklabs=round(seq(from=Sminlong,to=Smaxlong, by=(Smaxlong-Sminlong)/4),digits=2),
 	#y.ticklabs=round(seq(from=Sminlat,to=Smaxlat, by=(Smaxlat-Sminlat)/4),digits=2),
 	#z.ticklabs=c(-80,-60,-40,-20,0),
 	cex.lab=1,
 	cex.axis=1,
 	cex.symbols=1,
 	#lab=c(3, 4),
 	lab.z=5
 )


###################################################
### code chunk number 11: TrackReconstruction.Snw:586-589
###################################################
 DateData<-seq(ISOdatetime(2009,07,14,00,00,00, tz="GMT"),ISOdatetime(2009,07,28,00,00,00,
 	tz="GMT"), by="min")
 head(DateData)



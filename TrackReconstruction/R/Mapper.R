Mapper <-function(inFile, gpsFile, bathyFile, ExpFact=500, minlat=51, maxlat=60, minlong=-177, maxlong=-163, Title="Main")
{
	#This part creates a low resolution image
	#image.plot(bathyFile,
	#	col=c(rev(Bathymetry.palatte(200)),terrain.colors(100)),#gray(0:20/20),
	#	breaks=round(c(seq(from=min(bathyFile),to=0,length.out=201),seq(from=max(bathyFile)/101,to=max(bathyFile),length.out=100))),
	#	smallplot=2 #plots legend off x axis
	#	)
	
	#find limits of new graph
	limits<-GraphLimits(inFile)
	Sminlat=limits$miny
	Smaxlat=limits$maxy
	Sminlong=limits$minx
	Smaxlong=limits$maxx
	
	#Create smoothed image data using functions from fields package
	vec<-as.vector(bathyFile)
	matt=matrix(vec,nrow=((maxlong-minlong)*120),ncol=((maxlat-minlat)*120))
	lat=seq(((Sminlat-minlat)*120),((Smaxlat-minlat)*120),,ExpFact)#((Smaxlat-Sminlat)*120*10))
	long=seq(((Sminlong-minlong)*120),((Smaxlong-minlong)*120),,ExpFact)#((Smaxlong-Sminlong)*120*10))
	loc<-make.surface.grid(list(long,lat))
	objj<-list(x=1:((maxlong-minlong)*120),y=1:((maxlat-minlat)*120),z=matt)
	interp.surface(objj, loc)-> look
	
	#Create palatte for depth colors
	Bathymetry.palatte<-colorRampPalette(brewer.pal(9, "Blues"),bias=3)
	
	#plot the smoothed image	
	image.plot(as.surface(loc, look),
		col=c(rev(Bathymetry.palatte(200)),terrain.colors(100)),#gray(0:20/20),
		breaks=round(c(seq(from=min(bathyFile),to=0,length.out=201),seq(from=max(bathyFile)/101,to=max(bathyFile),length.out=100))),
		main=Title,
		#smallplot=10,#plots legend off x axis
		#lab.breaks=c("","","","","","","")
		xaxt="n",
		yaxt="n",
		xlab="Longitude",
		ylab="Latitude"
	)
	points((inFile$Longitude-minlong)*120,(inFile$Latitude-minlat)*120,
		#col=ifelse(inFile$SunTimes==1,"red","black"),
		pch=".",
		cex=1
	)
	points((gpsFile$Longitude-minlong)*120,(gpsFile$Latitude-minlat)*120,
		col="red",
		pch="*",
		cex=1
	)	
	axis(1,at=seq(from=((Sminlong-minlong)*120),to=((Smaxlong-minlong)*120), by=(((Smaxlong-minlong)*120)-((Sminlong-minlong)*120))/5),labels=round(seq(from=Sminlong,to=Smaxlong, by=(Smaxlong-Sminlong)/5),digits=1))
	axis(2,at=seq(from=((Sminlat-minlat)*120),to=((Smaxlat-minlat)*120), by=(((Smaxlat-minlat)*120)-((Sminlat-minlat)*120))/5),labels=round(seq(from=Sminlat,to=Smaxlat, by=(Smaxlat-Sminlat)/5),digits=1))
	#dev.off()
}

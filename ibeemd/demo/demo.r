library(rgdal)

# polygon data
mydata <- system.file("extdata/simu.shp", package = "ibeemd")

layer <- basename(mydata)
layer <- substr(layer, 1, nchar(layer)-4)
mydataDf <- readOGR(dsn=mydata, layer=layer)
#spplot(mydataDf)

rslt <- iBEEMD(
		spPolysDf = mydataDf, 
		valueField = "value",  
		nMaxIMF = 10, 
		tolSift = 0.05,
		neemd = 500,
		wnsd = 0.05,
		fmodel = "thinplate",
		fig = TRUE)

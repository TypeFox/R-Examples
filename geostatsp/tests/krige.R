# land type, categorical variable
library("geostatsp")
data("swissRain")
swissRain$lograin = log(swissRain$rain)
swissRain$elevation = extract(swissAltitude, swissRain)
swissAltitude[1:50,1:50] = NA

swissRaster = raster(extent(swissBorder), ncols=20, nrows=20, 
		crs=swissRain@proj4string)	


swissRain$land = raster::extract(swissLandType, swissRain)
# get rid of land types with few observations
landTable = table(swissRain$land)
landTable = as.numeric(names(landTable)[landTable > 5])
swissRain2 = swissRain [swissRain$land %in% landTable, ]

swissFit3 = likfitLgm(
    data=swissRain2[1:60,], 
		formula=lograin~ elevation + factor(land),
		param=c(range=46500, nugget=0.05,shape=1,  
				anisoAngleDegrees=35, anisoRatio=12),
		paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"),
		parscale = c(range=5000,nugget=0.01, 
				anisoRatio=1,anisoAngleDegrees=5)
)

swissKrige3 = krigeLgm(
    data=swissRain2[1:60,], 
		formula = swissFit3$model$formula,
		param=swissFit3$param, 
		covariates = list(elevation = swissAltitude,land=swissLandType),
		grid = swissRaster, expPred=TRUE)

pdf("krige3.pdf")
plot(swissKrige3[["predict"]])	
plot(swissBorder, add=TRUE)
dev.off()

if(interactive()  | Sys.info()['user'] =='patrick') {

# now change land to a factor
landTypes = swissLandType@data@attributes[[1]]
# remove land types with no observations on them
landTypes=landTypes[landTypes[,1] %in% unique(swissRain2$land),]

swissRain2$landFac = factor(swissRain2$land, 
		levels=landTypes[,1],
		labels=landTypes[,2])

swissFit4 = likfitLgm(data=swissRain2[1:60,], 
		formula=lograin~ elevation + landFac,
		param=c(range=46500, nugget=0.05,shape=1,  
				anisoAngleDegrees=35, anisoRatio=12),
		paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"),
		parscale = c(range=5000,nugget=0.01, 
				anisoRatio=1,anisoAngleDegrees=5)
)
swissKrige4 = krigeLgm(data=swissRain2[1:60,], formula = swissFit4$model$formula,
		param=swissFit4$param, 
		covariates = list(elevation = swissAltitude,landFac=swissLandType),
		grid = swissRaster,expPred=TRUE )




pdf("krige4.pdf")
plot(swissKrige4[["predict"]])	
plot(swissBorder, add=TRUE)
dev.off()

  

swissRain2$landFac2 = as.character(swissRain2$landFac)

swissFit5= likfitLgm(lograin~ elevation + factor(landFac2),
		data=swissRain2, 
		param=c(range=46500, nugget=0.05,shape=1,  
				anisoAngleDegrees=35, anisoRatio=12),
		paramToEstimate = c("range","nugget", 
				"anisoAngleDegrees", "anisoRatio"),
		parscale = c(range=5000,nugget=0.01, 
				anisoRatio=1,anisoAngleDegrees=5)
)


swissKrige5 = krigeLgm(data=swissRain2[1:60,], 
		formula = swissFit5$model$formula,
		param=swissFit5$param, 
		covariates = list(elevation = swissAltitude,landFac2=swissLandType),
		grid = swissRaster,expPred=TRUE)
pdf("krige5.pdf")
plot(swissKrige5[["predict"]])	
plot(swissBorder, add=TRUE)
dev.off()

}
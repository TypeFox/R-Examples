library("geostatsp")

model <- c(var=5, range=20,shape=0.5)

if (requireNamespace("RandomFields", quietly = TRUE)) { 
	myraster = raster(nrows=50,ncols=50,xmn=100,ymn=100,xmx=110,ymx=110, 
			crs="+init=epsg:2081")
} else {
  myraster = raster(nrows=20,ncols=20,xmn=100,ymn=100,xmx=110,ymx=110, 
		crs="+init=epsg:2081")
}

set.seed(0)
simu = RFsimulate(rbind(a=model, b=model+0.1), 
		x=myraster, n=3
		)

    set.seed(0)
    simu2 = RFsimulate(rbind(a=model, b=model+0.1), 
        x=as(myraster,"SpatialPixels"),
        n=3
    )
    

par(mfrow=c(length(names(simu2)),2))

for(D in 1:length(names(simu2))) {
			plot(simu[[D]])
			plot(raster(simu2,layer=D))
}


if(interactive()  | Sys.info()['user'] =='patrick') {
  simu2 = RFsimulate(rbind(a=model, b=model+0.1), 
      x=as(myraster,"SpatialPoints")[
          sample(ncell(myraster), 12)
          ,]
  )
  
  simu2 = RFsimulate(rbind(a=model, b=model+0.1), 
      x=as(myraster,"SpatialGrid")
  )
  
  for(Dn in c(1,3)) {
    set.seed(0) 
    simu <- RFsimulate(model, x=myraster, n=Dn)
    set.seed(0) 
    simu2 <- RFsimulate(model, x=as(myraster,"SpatialPixels"), n=Dn)
    
    print(proj4string(simu))
    print(proj4string(simu2))
    
    par(mfrow=c(nlayers(simu),2))
    for(D in 1:nlayers(simu)) {
      plot(simu[[D]])
      plot(raster(simu2,layer=D))
    }
  }
  
  
  
data("swissRain")
swissRain$sqrtrain = sqrt(swissRain$rain)

# estimate parameters


# isotropic
	swissRes =  lgm(data=swissRain, 
			grid=20, formula="sqrtrain",
			covariates=swissAltitude,   
			shape=1, fixShape=TRUE,
			aniso=FALSE, fixNugget=FALSE,
			nuggetInPrediction=FALSE
	)
	
	
 # anisotropic
swissRes =  lgm("sqrtrain",
		swissRain, grid=20, 
		covariates=swissAltitude,   
		shape=1, fixShape=TRUE,
		aniso=TRUE, fixNugget=FALSE,
		nuggetInPrediction=FALSE
)
 


# uncoinditional simulation
swissSim = RFsimulate(
		model=swissRes$param,
		x=swissRes$predict,
		n=3
)


# simulate from the random effect conditional on
#   the observed data

swissSim = RFsimulate(
		model=swissRes$param,
		data=swissRes$data[,'resid'],
		x=swissRes$predict,
		err.model=swissRes$param["nugget"],
		n=3
)

# plot the simulated random effect
plot(swissSim[[1]])
plot(swissBorder, add=TRUE)

# now with multiple parameter sets 
swissSim = RFsimulate(model=
				rbind(
						swissRes$param,
						swissRes$param*0.99),
		data=swissRes$data[,'resid'],
		x=swissRes$predict,
		err.model=c(1, 0.99)*swissRes$param["nugget"],
		n=3
)
# plot the simulated random effect
plot(swissSim[[1]])
plot(swissBorder, add=TRUE)

# and multiple simulations
# now with multiple datasets 
swissSim = RFsimulate(model=
				rbind(
						swissRes$param,	
						0.99*swissRes$param,
						1.01*swissRes$param),
		data=swissRes$data[,rep('resid',3)],
		err.model=c(1, 0.99, 1.01)*swissRes$param["nugget"],
		x=swissRes$predict,
		n=3
)
# plot the simulated random effect
plot(swissSim[[1]])
plot(swissBorder, add=TRUE)




# create a small raster of elevation data
swissAltSmall = aggregate(swissAltitude,fact=5)
swissAltSmall = resample(swissAltSmall, swissSim)

# calculate the fixed effects portion of the rainfall process
rainMean = swissRes$param["(Intercept)"] +
		swissRes$param[ "CHE_alt" ] * swissAltSmall

# define a function to identify the location of maximum rainfall	
maxRainLocation = function(x, ...) {
	rain =  (rainMean + x)^2
	 which.max(rain)
}


swissLocation = raster::cellStats(swissSim,   maxRainLocation)
swissLocation = xyFromCell(swissSim, swissLocation)
plot(swissRes$predict[["predict"]])
plot(swissBorder, add=TRUE)
points(swissLocation)
}
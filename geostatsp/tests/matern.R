options(CBoundsCheck=TRUE)
library("geostatsp")

param = c(range=1, shape=1.5,	anisoRatio=2, anisoAngleDegrees=-25)

matern(c(0, 0.001, 100000), param=param)

#x=c(0, 0.001, 100000);param=c(param, variance=1)

#resultFull = .C("matern", as.double(x), as.integer(length(x)),
#		as.double(param["range"]), as.double(param["shape"]),
#		as.double(param["variance"]))


# example with raster
myraster = raster(nrows=40,ncols=60,xmn=-3,xmx=3,ymn=-2,ymx=2)

# plot correlation of each cell with the origin
myMatern = matern(myraster, y=c(0,0), param=param)
as.matrix(myMatern)[1:3,1:3]



bob = function(x) {
thepar = attributes(x)$param
pdf(tempfile("matern", tmpdir=".", fileext=".pdf"))
plot(x, main=
				paste(
				paste(names(thepar), thepar, sep="="),
				collapse=", "),cex.main=0.5
)
dev.off()
}

bob(myMatern)


bob(matern(myraster, y=c(1,-0.5), 
				param =	c(range=1, shape=1.5,	anisoRatio=2, anisoAngleDegrees=25)			
						)
)

bob(matern(myraster,y= c(0,0), 
				param =	c(range=1, shape=25.1,	anisoRatio=2, anisoAngleDegrees=-25)			
		)
)

bob(matern(myraster,y= c(0,0), 
				param =	c(range=0, shape=1.5,	anisoRatio=2, anisoAngleDegrees=-25)			
		)
)

bob(matern(myraster, y=c(0,0), 
				param =	c(range=100000, shape=1.5,	anisoRatio=2, anisoAngleDegrees=-25)			
		)
)


				
				

# correlation matrix for all cells with each other
myraster = raster(nrows=4,ncols=6,xmn=-3,xmx=3,ymn=-2,ymx=2)
myMatern = matern(myraster, param=c(range=0, shape=2))

dim(myMatern)
myMatern[1:3,1:3]



param = c(range=0.06, shape=1.5,	anisoRatio=2, anisoAngleDegrees=-25)
mypoints = SpatialPointsDataFrame(cbind(runif(10), runif(10)),data=data.frame(id=1:10))
matern(mypoints, param=param)

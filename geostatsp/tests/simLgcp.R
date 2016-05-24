library('geostatsp')
mymodel = c(mean=-1.5, variance=1, 
				range=2, shape=2)

myraster = raster(nrows=15,ncols=15,xmn=0,xmx=10,ymn=0,ymx=10)

# some covariates, deliberately with a different resolution than myraster
covA = covB = myoffset = raster(extent(myraster), 10, 10)
values(covA) = as.vector(matrix(1:10, 10, 10))
values(covB) = as.vector(matrix(1:10, 10, 10, byrow=TRUE))
values(myoffset) = round(seq(-1, 1, len=ncell(myoffset)))

myCovariate = list(a=covA, b=covB, offsetFooBar = myoffset)

set.seed(0)
myLgcp=simLgcp(mymodel, myCovariate, 
    betas=c(a=-0.1, b=0.25), 
	offset='offsetFooBar',
	rasterTemplate=myraster)

if(requireNamespace("INLA", quietly=TRUE)) {
res = lgcp(data=myLgcp$events, 
		formula = ~ a + b + offset(offsetFooBar),
		grid=squareRaster(myoffset, 15), 
		covariates=myCovariate,
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41)),
		control.mode=list(theta=c(0.022, -0.5),restart=TRUE)
)

res$parameters$summary[,c(1,3,5)]

lgcpRoc =  spatialRoc(
    res, 
	rr=log(10), 
	truth=myLgcp, 
	random=FALSE)

head(lgcpRoc)

plot(lgcpRoc[,'onemspec'] , 
	lgcpRoc[,'sens'], 
	type='l', 
	xlim=c(0,1), ylim=c(0,1),
	ylab='sensitivity', xlab='1-specificity'
)

}

 





	
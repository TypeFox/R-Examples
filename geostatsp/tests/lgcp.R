havePackages = c(
    'INLA' = requireNamespace('INLA', quietly=TRUE)
)

print(havePackages)

# as in example
require('geostatsp')
myPoints = SpatialPoints(cbind(rbeta(100,2,2), rbeta(100,3,4)))
myPoints@bbox = cbind(c(0,0), c(1,1))

mycov = raster(matrix(rbinom(100, 1, 0.5), 10, 10), 0, 1, 0, 1)
names(mycov)="x1"

if(all(havePackages)) {
  
res = lgcp(data=myPoints, grid=20, covariates=mycov,
		formula=~factor(x1),
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

	
# intercept only

res = lgcp(data=myPoints, grid=20, covariates=mycov,
		formula=~1,
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

if(interactive()  | Sys.info()['user'] =='patrick') {

# some missing values

temp = values(mycov)
temp[1:4] = NA
values(mycov) = temp

res = lgcp(data=myPoints, grid=20, covariates=mycov,
		formula=~factor(x1),
		priorCI=list(sd=c(0.9, 1.1), range=c(0.4, 0.41))
)

plot(res$raster[["predict.exp"]])
plot(myPoints,add=TRUE,col="#0000FF30",cex=0.5)

data('murder')
data('torontoPop')
myCov = list(
    pop=torontoPdens,
    inc = log(torontoIncome)
)

formula = ~ inc + offset(pop, log=TRUE)

resL=lgcp(formula, data=murder, 
    grid=squareRaster(murder, 30),
    covariates=myCov,
    border=torontoBorder)

resO = lgcp( ~ inc + pop, 
    data=murder, 
    grid=squareRaster(murder, 30),
    covariates=list(inc=myCov$inc, pop=log(myCov$pop)),
    border=torontoBorder)

rbind(resL$param$summary, resO$param$summary)



}	
}


# check spdfToBrock

if(require('diseasemapping', quiet=TRUE)){
	
	
	data('kentucky')
	
	popList = list(
			'2002' = kentucky[,c('M.50', 'M.55', 'M.60')],
			'2005' = kentucky[,c('F.50', 'F.55', 'F.60')]
			)
	for(D in names(popList))
		names(popList[[D]]) = paste('expected_', seq(as.numeric(D)-1, len=3), sep='')
	
	popBrick = spdfToBrick(
			x=popList,
			template=squareRaster(kentucky, 10),
    	logSumExpected=FALSE
			)
	sum(popList[['2002']]$expected_2001, na.rm=TRUE)		
	sum(values(popBrick[['expected_2001']]), na.rm=TRUE)*prod(res(popBrick))
	
	popBrick2 = spdfToBrick(
			x=popList,
			template=squareRaster(kentucky, 10),
    	logSumExpected=TRUE
	)
	
	sum(unlist(lapply(popList, function(qq) apply(qq@data, 2, sum, na.rm=TRUE))))
	sum(exp(values(popBrick2)), na.rm=TRUE)*prod(res(popBrick2))
			
}
	